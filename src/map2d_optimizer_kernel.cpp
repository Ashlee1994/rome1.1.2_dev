/***************************************************************************
 *
 * Authors: "Bevin Brett"
 * Intel Corporation
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include "map2d_optimizer_kernel.h"

#include <atomic>

static const bool poolKernels = true;
static const bool poolKernels_uop  = poolKernels & true;
static const bool poolKernels_bp   = poolKernels & true;
static const bool poolKernels_gasd = poolKernels & true;

// #define CHECK_NANS_ETC
//		used for debugging

static bool use3x3 =
#ifdef MY_MACHINE_HAS_AVX512
	true;
#else
	false;
#endif

namespace Map2dOptimizer_Kernel {
	void setKernel3x3(bool to) {
#ifdef MY_MACHINE_HAS_AVX512
		use3x3 = to;
#else
		if (to) std::cerr << "Map2dOptimizer_Kernel::setKernel3x3 machine does not support AVX512" << std::endl;
#endif
	}
};

// Compiler command lines:
// Linux:
//	g++  -o ML2D_Microkernel ML2D_Microkernel.cpp  -std=c++0x -m64 -mavx|-mavx2 -lrt -g -fopenmp
//	icpc -o ML2D_Microkernel ML2D_Microkernel.cpp  -std=c++0x -m64 -mavx|-mavx2 -lrt -g -fopenmp
//		Note:  -mavx is for Sandybridge/Ivybridge ; -mavx2 is for Haswell / Broadwell
// 
// Windows - from the 64-bit compilation environment
//	cl ML2D_Microkernel.cpp -Zi
//
#include <assert.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>

#ifdef _WIN32

#include <intrin.h> 
#define POPCNT64   __popcnt64

#else

#include <x86intrin.h> 
#define POPCNT64   _popcnt64

#endif

#if !defined(CHECK_NANS_ETC)
static void checkNotNegOrNan(double value) { }
static float checkNotNanOrInf(float value) { return value; }
#else
static void checkNotNegOrNan(double value) {
	if (value < 0 || isnan(value)) {
		std::cerr << "Dist::add value is bad " 
			<< " value:"  << value
			<< std::endl;
		std::cout << "Dist::add value is bad " 
			<< " value:"  << value
			<< std::endl;
		exit(1);
	}
}

static float checkNotNanOrInf(float value) {
	if (isinf(value)) {
		if (value < 0) return -std::numeric_limits<float>::max();
		else		   return +std::numeric_limits<float>::max();
	}
	if (isinf(value) || isnan(value)) {
		std::cerr << "Capture value is bad " 
			<< " value:"  << value
			<< std::endl;
		std::cout << "Capture value is bad " 
			<< " value:"  << value
			<< std::endl;
		exit(1);
	}
	return value;
}
#endif

namespace Map2dOptimizer_Kernel {
using namespace MapOptimizer_base_new;

typedef __int64			 S64;
typedef unsigned __int64 U64;

static const int maxTranslationIndexsCapacity = 5;

template<typename T>
void initVectorWithNULLs(std::vector<T*>& v) {
	for (auto i = v.begin(); i != v.end(); i++) (*i) = NULL;
}

template<typename T>
void assignLhsIsSameOrNULL(T* & lhs, T* rhs) {
	assert(!lhs || (lhs == rhs));
	lhs = rhs;
}

template<typename T>
void assignLhsIsSameOrZero(T & lhs, T rhs) {
	assert(!lhs || (lhs == rhs));
	lhs = rhs;
}

#define mallocAlignedDoubleNotZeroed(size) mallocDoubles(size)
#define mallocAlignedDoubleZeroed(size)    mallocZeroedDoubles(size)
#define mallocAlignedFloatNotZeroed(size)  mallocCacheAligned(float,size)
#define mallocAlignedFloatZeroed(size)     mallocZeroedFloats(size)

inline bool bitScanForward(U64& whereFound, U64 mask) {
#ifdef _WIN32
	unsigned long wf;
	auto result = _BitScanForward64(&wf, mask);
	whereFound = wf;
	return result;
#elif defined(__APPLE__)
    whereFound = ffsl(mask);
    //least significant bit is position 1
    whereFound--;
    return mask != 0;
#else
//	whereFound = __bsfq(mask);
//	return mask != 0;
    whereFound = ffsl(mask);
    //least significant bit is position 1
    whereFound--;
    return mask != 0;
#endif
}

//============================================================================================================
// To eliminate many locks, some structures have an element per thread
// However tandem execution means two different o/s threads come through here, maybe having the same omp_thread_number each (omp sequential).
// Perhaps also OMP is also capable of reusing idle threads so a thread may require multiple instances of an entity
// So each list is locked with os locks.
//
static std::atomic<size_t> listPerThreadTemplateCells = 0;

template <class Cell>
class ListPerThreadTemplate {
public:
	ListPerThreadTemplate() : lists(omp_get_max_threads()) {}
	~ListPerThreadTemplate() {
		for (auto i = lists.begin(); i != lists.end(); i++) {
			auto& listHead = *i;
			while (auto head = listHead.head) {
				listHead.head = head->next;
				sDelete(head);
			}
		}
	}
	Cell* acquire() {
		auto thread_num = omp_get_thread_num();
		auto& listHead = lists[thread_num];
		ScopedAcquire scopedAcquire(listHead, __FILE__, __LINE__);
		for (auto p = listHead.head; p; p = p->next) {
			if (!p->inUse) { p->inUse = true; return &p->cell; }
		}
		listPerThreadTemplateCells++;
		std::cerr << "kernel ListPerThreadTemplate making a new ListElement, #" << listPerThreadTemplateCells << std::endl;
		auto p = ListElement::make(listHead);
		p->inUse = true;
		return &p->cell;
	}
	template <class Key>
	void release(Key key) {
		auto thread_num = omp_get_thread_num();
		auto& listHead = lists[thread_num];
		ScopedAcquire scopedAcquire(listHead, __FILE__, __LINE__);
		for (auto p = listHead.head; p; p = p->next) {
			if (!p->cell.matches(key)) continue;
			assert(p->inUse);
			p->inUse = false;
			return;
		}
		assert(false);
	}
private:
	struct ListElement;
	struct ListHead : public Lock {
		ListElement *	head;
		ListHead() : head(nullptr) {}
	};
	struct ListElement {
		ListElement *	next;
		bool			inUse;
		Cell			cell;
		ListElement(ListHead & listHead) : next(listHead.head), inUse(false) { listHead.head = this; }
		static ListElement* make(ListHead & listHead) {
#include "./util_heap_undefs.h"
			return sNewA(ListElement, (listHead));
#include "./util_heap_defs.h"
		}
	};
	std::vector<ListHead> lists;
};


//============================================================================================================
RotationState::RS RotationState::getRotationState(
	int   iclass,
	int   rotation) const {
	RotationState::RS result;
	auto& crs = currentRotationStates[iclass*numberOfRotations + rotation];
	return crs.rs;
}

void RotationState::setRotationState(
	int   iclass,
	int   rotation,
	RS    to,
	void* mainAddr,
	void* addAddr) {

    auto& crs = currentRotationStates[iclass*numberOfRotations + rotation];

#ifdef NDEBUG
#define CHECK(EXP) EXP
#else
#define CHECK(EXP) check(EXP,#EXP)
	auto check = [&](bool passed, const char* condition) {
		if (passed) return;
		std::cerr << "RotationState::setRotationState failed "  << condition << std::endl;
		std::cout << "RotationState::setRotationState failed "  << condition << std::endl;
		std::cerr << "crs.rs   = " << crs.rs << std::endl;
		std::cerr << "to       = " << to     << std::endl;
		std::cerr << "mainAddr = " << (void*)mainAddr << std::endl;
		std::cerr << "addAddr  = " << (void*)addAddr  << std::endl;
		EXIT_ABNORMALLY;
	};
#endif

	switch (to) {
    //case RS_unbuffered:			// Constructed here and never comes back here
    //	crs.mainAddr = 0;
    //	crs.addAddr  = 0;
    //	break;
    case RS_zeroedBuffer:			// the shifted weighted numbers are being accumulated
    	CHECK(crs.rs == RS_unbuffered);
    	CHECK(!crs.mainAddr);
    	CHECK(!crs.addAddr );
    	crs.mainAddr = mainAddr;
    	crs.addAddr  = addAddr;
    	break;
    case RS_captured:				// the shifted weighted numbers are being put in the accumulator by the BackProjection_Kernel
    	CHECK(crs.rs == RS_zeroedBuffer 
			|| crs.rs == RS_captured 
			|| crs.rs == RS_addPending);
    	CHECK(!mainAddr || crs.mainAddr == mainAddr);
    	CHECK(crs.addAddr  == addAddr);
    	break;
    case RS_addPending:				// the accumulator is captured by a kernel
    	CHECK(crs.rs == RS_captured);
    	CHECK(!mainAddr || crs.mainAddr == mainAddr);
    	CHECK(crs.addAddr == addAddr);
    	break;
    case RS_added:					// the accumulator is released
    	CHECK(crs.rs == RS_zeroedBuffer 
		   || crs.rs == RS_addPending);
    	CHECK(crs.mainAddr == mainAddr);
    	CHECK(crs.addAddr  == addAddr);
    	break;
    default:
    	CHECK(false);
    }
#undef CHECK

    crs.rs = to;
}


enum SchedulerKind            { SchedulerKind_BackProjection,  SchedulerKind_GetAllSquaredDifferences,  SchedulerKind_UpdateOtherParams, SchedulerKind__end};
static const char* toStr(SchedulerKind sk) {
	static const char* strs[] = {"SchedulerKind_BackProjection","SchedulerKind_GetAllSquaredDifferences","SchedulerKind_UpdateOtherParams"};
	return strs[sk];
}

class FloatsForRTPool {
	struct Cell {
		S64		len; 
		float*	v; 
		Cell() : len(0), v(nullptr) {} 
		~Cell() { deallocate(); }
		bool matches(float* v) { return this->v == v; }
		void allocate(U64 len) {
			assert(!v);
			v = mallocAlignedFloatNotZeroed(len);
			this->len = len;
		}
		void deallocate() {
			aFree(v);
			len = 0;
		}
	};
	ListPerThreadTemplate<Cell> cellsPerThread;
public:
	FloatsForRTPool() {}
	float* allocate(S64 len) {
		auto & cell = *cellsPerThread.acquire();
		if (cell.len < len) {
			cell.deallocate();
			cell.allocate(len);
		}
		return cell.v;
	}
	void deallocate(float* v) {
		cellsPerThread.release(v);
	}
};

FloatsForRTPool floatsForRTPools[SchedulerKind__end];


class FloatsForRT {
public:
	const SchedulerKind schedulerKind;
	FloatsForRT(SchedulerKind schedulerKind, S64 numberOfRotations, S64 numberOfTranslations) 
	  : schedulerKind(schedulerKind), numberOfRotations(numberOfRotations), numberOfTranslations(numberOfTranslations), v(floatsForRTPools[schedulerKind].allocate(numberOfRotations*numberOfTranslations)) 
	{
		init();
	}
	~FloatsForRT() {
		floatsForRTPools[schedulerKind].deallocate(v);
	}
	void init() {
	}
	int index(int rotation, int translation) const {
		assert(rotation    < numberOfRotations);
		assert(translation < numberOfTranslations);
		return translation*numberOfRotations + rotation;
	}
	float& operator[](S64 index) {
		assert(0 <= index && index < numberOfRotations*numberOfTranslations);
		return v[index]; 
	}
protected:
	S64 numberOfRotations; S64 numberOfTranslations;
	float* v;	// using std::vector introduced many locks
};

class SparseSampleIndexs {
public:
	SparseSampleIndexs(int numberOfSamples) : numberOfSamples(numberOfSamples) {}
	void init_Mresol_fine(const int* Mresol_fine) {
		int numberOfAcceptedSamples(0);
		for (int n = 0; n < numberOfSamples; n++) {
			if (Mresol_fine[n] <= -1) continue;
			numberOfAcceptedSamples++;
		}
		denseToSparse.resize(numberOfAcceptedSamples);
		numberOfAcceptedSamples = 0;
		for (int n = 0; n < numberOfSamples; n++) {
			if (Mresol_fine[n] < 0) continue;
			denseToSparse[numberOfAcceptedSamples++] = n;
		}
	}
	int numberOfAcceptedSamples() const { return denseToSparse.size(); }
	int sparseIndex(int denseIndex) const { return denseToSparse[denseIndex]; }
protected:
	int numberOfSamples;
	std::vector<int> denseToSparse;
};


class CapturedDataStatsHistogram {
	static const size_t scale = 1;
public:
	CapturedDataStatsHistogram() : buckets(65) { init(); }
	void init() {
		for (size_t i = 0; i < buckets.size(); i++) {
			buckets[i] = 0;
		}
	}
	void record(float opsPerAcq) {
		auto index = std::min(buckets.size()-1, size_t(opsPerAcq * scale));
		buckets[index]++;
	}
	void display() {
		int count = 0;
		for (size_t i = 0; i < buckets.size(); i++) {
			count += buckets[i];
		}
		std::cout << "CapturedDataStatsHistogram{";
		for (size_t i = 0; i < buckets.size(); i++) {
			if (buckets[i] == 0) continue;
			std::cout << " [" << std::setprecision(3) << float(i)/scale;
			if (i == buckets.size() - 1) std::cout << " or gtr";
			std::cout << "]=" 
				<< buckets[i]
				<< "(" << int(float(buckets[i])/count*100) << "%)";
		}
		std::cout << "}" << std::endl;
	}
private:
	std::vector<int> buckets;
};


class ICapturedDataStats {
public:
	virtual void captureMore(S64 moreElts) = 0;
	virtual void add(bool isAlreadyCaptured, int dataLen) = 0;
};


class CapturedDataStats : public ICapturedDataStats {
public:
//#define CAPTURE_STATS
#if !defined(CAPTURE_STATS)
	CapturedDataStats(CapturedDataStatsHistogram& capturedDataStatsHistogram) {}
	void init() {}
	void record() {}
	void fini() {}
	void add(bool isAlreadyCaptured, int dataLen) {}
	void captureMore(S64 elts) {}
	void operate(bool onCaptures, int numberOfOperations) {}
#else
	CapturedDataStats(CapturedDataStatsHistogram& capturedDataStatsHistogram) : capturedDataStatsHistogram(capturedDataStatsHistogram) {
		init();
	}

	void init() {
		uniqueCaptured = 0;
		duplicateCaptured = 0;
		len = 0;
		capturedElts = 0;
		numberOfOperationsOnCaptures = 0;
		numberOfOperationsOnOriginal = 0;
	}

	void record() {
		if (capturedElts == 0 || omp_get_thread_num() != 0) return;
		auto opsPerAcq = float(numberOfOperationsOnCaptures)/uniqueCaptured;
		capturedDataStatsHistogram.record(opsPerAcq);
	}

	void fini() {
		if (capturedElts == 0 || omp_get_thread_num() != 0) return;
		auto opsPerAcq = float(numberOfOperationsOnCaptures)/uniqueCaptured;

		bool badAmortization(numberOfOperationsOnCaptures > 0 && uniqueCaptured*10 >= numberOfOperationsOnCaptures);
		S64 expected = len*uniqueCaptured;
		std::cout << "CapturedDataStats {"
				<< " len:"		<< std::setw(4) << len
				<< " acq:"		<< std::setw(4) << uniqueCaptured 
				<< " dup:"		<< std::setw(4) << duplicateCaptured
				<< " opsOnCap:"	<< std::setw(4) << numberOfOperationsOnCaptures
				<< " opsOnOri:"	<< std::setw(4) << numberOfOperationsOnOriginal
				<< " ops/acq:"	<< std::setw(4) << opsPerAcq
				<< " capElts:"	<< std::setw(4) << capturedElts;
		if (capturedElts != expected) std::cout << " out of an expected:" << expected;
		std::cout << "}" << std::endl;
		if (badAmortization)          std::cout << "    ******************** Note: Not enough to amortize capture over";
		if (capturedElts > expected)  std::cout << "    ******************** Note: Too many captured";
		std::cout << std::endl;
	}

	~CapturedDataStats() {
	}

	void add(bool isAlreadyCaptured, int dataLen) {
		assert(!isAlreadyCaptured || uniqueCaptured>0);
		(isAlreadyCaptured ? duplicateCaptured : uniqueCaptured)++;
		assignLhsIsSameOrZero(len, dataLen);
	}

	void captureMore(S64 elts) {
		assert(len > 0);
		this->capturedElts += elts;
		if (capturedElts > len*uniqueCaptured) 
			std::cout << "******************** Note: Too many captured";
	}

	void operate(bool onCaptures, int numberOfOperations) {
		(onCaptures ? this->numberOfOperationsOnCaptures : this->numberOfOperationsOnOriginal) += numberOfOperations;
	}

private:
	int uniqueCaptured;
	int duplicateCaptured;
	int len;
	int capturedElts;
	int numberOfOperationsOnCaptures;
	int numberOfOperationsOnOriginal;
	CapturedDataStatsHistogram& capturedDataStatsHistogram;
#endif
};


static std::atomic<int> capturedDataNextIndex = 1;
static int const interestingCapturedDataIndex = 0;

class CapturedData {
	S64	         	almostConst_len;
	float*			almostConst_ptr;
	S64				capturedElts;	// used for debugging check

	ICapturedDataStats* stats;
	const char*     name;
	int				index;
	bool			interesting;
	enum Mode {Mode_undef, Mode_forSummingInto, Mode_const, Mode_constSparse, Mode_sqrt, Mode_MultipliedByWeightsSqrtForIS};
	Mode			mode;
	bool			zeroed;
	S64				capturedLen;	// supports deferring the capture until certain needed and until about to be used
	double*			originalNonConst;
	const double*	originalConst;

	// Info for the various modes
	const SparseSampleIndexs*	sparseSampleIndexs;
	CapturedData*				weightsSqrt;

	// Debug options
	static const bool deferCapture		= true;
	static const bool testDeferCapture	= deferCapture && false;

public:
	CapturedData() : almostConst_len(0), almostConst_ptr(NULL), capturedElts(0) {
		init(); 
		if (testDeferCapture) {
			index = capturedDataNextIndex++; 
		}
	}
	void init() {
		stats              = NULL;
		name               = "unnamed";
		index              = 0;
		interesting        = false;
		mode               = Mode_undef;
		zeroed			   = false;
		capturedElts	   = 0;
		capturedLen        = 0;
		originalNonConst   = NULL;
		originalConst      = NULL;
		sparseSampleIndexs = NULL;
		weightsSqrt		   = NULL;
	}
	~CapturedData() { 
		if (!almostConst_ptr) return;
		TUNING_SCOPE_STEP(CapturedData_dtr)
		aFree(almostConst_ptr); 
		almostConst_ptr = NULL; 
	}
	
	void setInteresting() { interesting = true; }
	bool isInteresting () { return interesting; }

	double*			ptrToSumInto      () const { return originalNonConst; }
	double const*   ptrToOriginalConst() const { return originalConst;    }

	bool isCapturedForSummingInto() const { return !!originalNonConst; }

	void captureForSummingInto(
		ICapturedDataStats* stats,
		const char* name,
		S64		dataLen,
		double* data) {
		this->name = name;
		setCaptureModeAndLen(stats, Mode_forSummingInto, dataLen);
		assignLhsIsSameOrNULL(originalNonConst, data);
		captureOrDeferUpto();
	}

	void capture(
		ICapturedDataStats* stats,
		const char*   name,
		S64			  dataLen,
		const double* data)
	{
		this->name = name;
		setCaptureModeAndLen(stats, Mode_const, dataLen);
		assignLhsIsSameOrNULL(originalConst, data);
		captureOrDeferUpto();
	}

	// The original formula was
	//			diffReal = x0_real - x1_real
	//			diffImag = x0_imag - x1_imag
	//
	//			(diffReal*diffReal + diffImag*diffImag) * weight
	//
	// But sometimes this goes out of range, with a zero weight, resulting in inf * 0 which is undefined
	// so I have rewritten this by eliminating the * weight by instead multiplying the x0_real etc. by the sqrt of the weight
	//
	void captureSqrt(
		ICapturedDataStats* stats,
		const char*   name,
		S64			  dataLen,
		const double* data)
	{
		this->name = name;
		setCaptureModeAndLen(stats, Mode_sqrt, dataLen);
		assignLhsIsSameOrNULL(originalConst, data);
		captureOrDeferUpto();
	}

	void captureMultipliedByWeightsSqrtForIS(
		ICapturedDataStats* stats,
		const char*   name,
		S64			  dataLen,
		const double* data,
		CapturedData& weightsSqrt)
	{
		this->name = name;
		this->weightsSqrt = &weightsSqrt;
		setCaptureModeAndLen(stats, Mode_MultipliedByWeightsSqrtForIS, dataLen);
		assignLhsIsSameOrNULL(originalConst, data);
		captureOrDeferUpto();
	}

	void captureSparse(
		ICapturedDataStats* stats,
		const char*   name,
		S64			  dataLen,
		const double* data,
		SparseSampleIndexs const & sparseSampleIndexs) {

		this->name = name;
		this->sparseSampleIndexs = &sparseSampleIndexs;
		setCaptureModeAndLen(stats, Mode_constSparse, dataLen);
		assignLhsIsSameOrNULL(originalConst, data);
		captureOrDeferUpto();
	}

	void sumIntoOriginal() {
		// The original is locked by the caller
		assert(!!originalNonConst);
		for (S64 i = 0; i < almostConst_len; i++) {
			originalNonConst[i] += double(almostConst_ptr[i]);
		}
	}

	float* ptrValidUpto(S64 newCapturedLen) {
		if (capturedLen < newCapturedLen) captureUpto(newCapturedLen);
		if (isInteresting()) {
			print(std::cout);
		}
		return almostConst_ptr;
	}

private:
	void setCaptureModeAndLen(ICapturedDataStats* stats, Mode mode, S64 dataLen) {
		bool isAlreadyCaptured = (this->mode != Mode_undef);
		assert(!isAlreadyCaptured || this->mode == mode);
		this->mode = mode;
		this->stats = stats;
		stats->add(isAlreadyCaptured, dataLen);
		assignLhsIsSameOrZero(almostConst_len, dataLen);
	}

	void __declspec(noinline) captureOrDeferUpto() {
		if      (!deferCapture)		 captureUpto(almostConst_len);
		else if (testDeferCapture) { captureUpto(almostConst_len); capturedLen = 0; }
	}

	void __declspec(noinline) captureUpto(S64 newCapturedLen) {
		TUNING_SCOPE_STEP_BEGIN(kernel::captureUpto)

		assert(0 <= newCapturedLen && newCapturedLen <= almostConst_len);
		newCapturedLen = std::max(newCapturedLen, capturedLen);

		if (almostConst_ptr == NULL) {
			TUNING_SCOPE_STEP(captureUpto_malloc)
			almostConst_ptr = mallocAlignedFloatNotZeroed(almostConst_len);
			if (testDeferCapture && (index == interestingCapturedDataIndex)) {
				std::cerr << "interestingCapturedDataIndex being captured" << std::endl;
			}
		} else if (testDeferCapture) {
			testCaptured(newCapturedLen);
			capturedLen = newCapturedLen;
			return;
		}

		auto moreElts = newCapturedLen - capturedLen;
		assert(capturedElts == capturedLen);
		assert(mode != Mode_undef);
		if (mode != Mode_forSummingInto) {
			stats->captureMore(moreElts);
		}
		capturedElts += moreElts;
		assert(capturedElts <= almostConst_len);

		switch (mode) {
		case Mode_forSummingInto: {
			if (!zeroed) {
				zeroed = true;
				for (S64 i = 0; i < almostConst_len; i++) {
					almostConst_ptr[i] = 0.0;
				}
			}
		} break;
		case Mode_const: {
			#pragma ivdep
			for (S64 i = capturedLen; i < newCapturedLen; i++) 
				almostConst_ptr[i] = checkNotNanOrInf(float(originalConst[i]));
		}	break;
		case Mode_constSparse: {
			#pragma ivdep
			for (S64 i = capturedLen; i < newCapturedLen; i++) 
				almostConst_ptr[i] = checkNotNanOrInf(float(originalConst[(*sparseSampleIndexs).sparseIndex(i)]));
		}	break;
		case Mode_sqrt: {
			#pragma ivdep
			for (S64 i = capturedLen; i < newCapturedLen; i++) 
				almostConst_ptr[i] = sqrtf(checkNotNanOrInf(float(originalConst[i])));
		}	break;
		case Mode_MultipliedByWeightsSqrtForIS: {
			weightsSqrt->captureUpto(newCapturedLen);
			#pragma ivdep
			for (S64 i = capturedLen; i < newCapturedLen; i++) 
				almostConst_ptr[i] = checkNotNanOrInf(float(originalConst[i])*weightsSqrt->almostConst_ptr[i]);
		}	break;
		default: assert(false);
		}	// switch

		capturedLen = newCapturedLen;
		TUNING_SCOPE_STEP_END
	}

	void testCaptured(S64 newCapturedLen) {
		if (!testDeferCapture) return;

		// The data was captured when the ptr was set
		// It SHOULD be the same as the data that would be captured here - check that it is

		auto check = [&](int i, float rhs) {
			float lhs = almostConst_ptr[i];
			if (lhs != rhs) {
				std::cerr << "Initial Captured data[" << i << "]:" << lhs << " != Deferred Captured data " << rhs << " for name " << name << " index:" << index << std::endl;
				std::cout << "Initial Captured data[" << i << "]:" << lhs << " != Deferred Captured data " << rhs << " for name " << name << " index:" << index << std::endl;
				EXIT_ABNORMALLY;
			}
		};

		switch (mode) {
		case Mode_forSummingInto: {
			for (S64 i = capturedLen; i < newCapturedLen; i++) check(i, 0.0f);
				// [0..capturedLen-1] may have already been added to
		} break;
		case Mode_const: {
			for (S64 i = 0; i < newCapturedLen; i++) check(i,float(originalConst[i]));
		}	break;
		case Mode_sqrt: {
			for (S64 i = 0; i < newCapturedLen; i++) check(i,sqrtf(checkNotNanOrInf(float(originalConst[i]))));
		}	break;
		case Mode_MultipliedByWeightsSqrtForIS: {
			if (weightsSqrt->capturedLen < newCapturedLen) {
				std::cerr << "weightsSqrt has only been captured up to " << capturedLen << " instead of " << newCapturedLen << name << std::endl;
				std::cout << "weightsSqrt has only been captured up to " << capturedLen << " instead of " << newCapturedLen << name << std::endl;
				EXIT_ABNORMALLY;
			}
			weightsSqrt->testCaptured(newCapturedLen);
			for (S64 i = 0; i < newCapturedLen; i++) check(i,checkNotNanOrInf(float(originalConst[i]*weightsSqrt->almostConst_ptr[i])));
		}	break;
		default: assert(false);
		}	// switch
	}

	void print(std::ostream & os) {
		os << "capture of " << name << ", showing sums, not individual values" << std::endl;
		double sum = 0.0;
		S64 showSumOf = 1;
		for (S64 i = 0; i < capturedLen; i++) {
			sum += almostConst_ptr[i];
			if (i+1 != showSumOf) continue;
			std::cout << i << ":" << sum << " ";
			showSumOf *= 2;
		}
		std::cout << std::endl;
	}
};


struct SingletonStep {
	int rotation;
	int translation;
	void print(const char* why) const {
		std::cout << why 
			<< " rotation:"    << rotation
			<< " translation:" << translation
			<< std::endl;
	}
};

struct MergedStep {
	int rotationBegin;
	int rotationLen;
	int translationIndexsLen;
	int translationIndexs[maxTranslationIndexsCapacity];

	void print(const char* why) const {
		std::cout << why << " rotations:";
		for (int i = 0; i < rotationLen; i++) std::cout << rotationBegin+i << " ";
		std::cout << "translations:";
		for (int i = 0; i < translationIndexsLen; i++) std::cout << translationIndexs[i] << " ";
		std::cout << std::endl;
	}
};

class StepsPool {
	// The Plan's allocating and deallocating steps showed up very high in the Locks and Waits
	// and this gets rid of that...
	struct Cell {
		U64				stepsCapacity;
		SingletonStep*	singletonSteps;
		MergedStep*		mergedSteps;
		Cell() : stepsCapacity(0), singletonSteps(NULL), mergedSteps(NULL) {}
		~Cell() { deallocate(); }
		bool matches(SingletonStep*	singletonSteps) { return this->singletonSteps == singletonSteps; }
		void allocate(U64 stepsCapacity) {
			assert(!singletonSteps);
			singletonSteps = vNew(SingletonStep,stepsCapacity);
			mergedSteps    = vNew(MergedStep   ,stepsCapacity);
			this->stepsCapacity = stepsCapacity;
		}
		void deallocate() {
			vDelete(mergedSteps);
			vDelete(singletonSteps);
			stepsCapacity = 0;
		}
	};
	ListPerThreadTemplate<Cell> cellsPerThread;
public:
	StepsPool() {}
	~StepsPool() {}
	void allocate  (U64 stepsCapacity, SingletonStep*& singletonSteps, MergedStep*& mergedSteps) {
		auto & cell = *cellsPerThread.acquire();
		if (cell.stepsCapacity < stepsCapacity) {
			cell.deallocate();
			cell.allocate(stepsCapacity);
		}
		singletonSteps = cell.singletonSteps;
		mergedSteps    = cell.mergedSteps;
	}
	void deallocate(SingletonStep*& singletonSteps, MergedStep*& mergedSteps) {
		cellsPerThread.release(singletonSteps);
		singletonSteps = nullptr;
		mergedSteps    = nullptr;
	}
};
static StepsPool stepsPool;


class MustDoPool {
	// This is part of getting rid of almost all cross-thread locks
	struct Cell { 
		int		numberOfRotations;
		int		numberOfTranslations;
		S64		len; 
		int*	rotationUseCount; 
		int*	translationUseCount;
		U64*	bits;
		Cell() : numberOfRotations(0), numberOfTranslations(0), len(0), rotationUseCount(NULL), translationUseCount(NULL), bits(0) {} 
		~Cell() { deallocate(); }
		bool matches(int* rotationUseCount) { return this->rotationUseCount == rotationUseCount; }
		void allocate(int numberOfRotations, int numberOfTranslations, S64 len) {
			assert(!rotationUseCount); assert(!translationUseCount); assert(!bits);
			rotationUseCount    = vNew(int,numberOfRotations);		this->numberOfRotations    = numberOfRotations;
			translationUseCount = vNew(int,numberOfTranslations);	this->numberOfTranslations = numberOfTranslations;
			bits                = vNew(U64,len);					this->len				   = len;
		}
		void deallocate() {
			vDelete(rotationUseCount)   ; numberOfRotations    = 0;
			vDelete(translationUseCount); numberOfTranslations = 0;
			vDelete(bits)               ; len					 = 0;
		}
	};
	ListPerThreadTemplate<Cell> cellsPerThread;
public:
	MustDoPool() {}
	void allocate(
		int* & rotationUseCount,
		int* & translationUseCount,
		U64* & bits,
		const int numberOfRotations, 
		const int numberOfTranslations,
		const S64 len) {
		auto & cell = *cellsPerThread.acquire();
		if (cell.numberOfRotations < numberOfRotations || cell.numberOfTranslations < numberOfTranslations || cell.len < len) {
			cell.deallocate();
			cell.allocate(numberOfRotations, numberOfTranslations, len);
		}
		rotationUseCount    = cell.rotationUseCount;
		translationUseCount = cell.translationUseCount;
		bits				= cell.bits;
	}
	void deallocate(int* & rotationUseCount, int* & translationUseCount, U64* & bits) {
		cellsPerThread.release(rotationUseCount);
		rotationUseCount = nullptr;
		translationUseCount = nullptr;
		bits = nullptr;
	}
};
MustDoPool mustDoPools[SchedulerKind__end];


class Scheduler {
public:

	// The data being operated on
	//
	const SchedulerKind schedulerKind;
	const int translationIndexsCapacity;

	const S64 numberOfRotations;
	const S64 numberOfTranslations;
	const S64 numberOfSamples;

	CapturedDataStatsHistogram capturedDataStatsHistogram;
	CapturedDataStats capturedDataStats;

	static bool showPlan(int iter) {
		if (1) return false;
		if (omp_get_thread_num() != 0) return false;
		static int prevIter   = -1;
		static U64 nextSample = 1;
		static U64 count;
		if (prevIter != iter) {
			prevIter = iter;
			nextSample = 1;
			count      = 0;
		}
		if (++count != nextSample) return false;
		nextSample = (nextSample<100) ? nextSample+1 : nextSample*2; 
		std::cout << "Showing plans for " << count << " kernels so far during iter " << iter << ", next will be at " << nextSample << std::endl;
		return true;
	}
	bool _showPlan;

	Scheduler(
		SchedulerKind schedulerKind,
		const int translationIndexsCapacity,
		const int numberOfRotations,
		const int numberOfTranslations,
		const int numberOfSamples) 
	  :
		schedulerKind(schedulerKind),
	    translationIndexsCapacity(translationIndexsCapacity),
		numberOfRotations   (numberOfRotations), 
		numberOfTranslations(numberOfTranslations), 
		numberOfSamples(numberOfSamples),
		capturedDataStats(capturedDataStatsHistogram),
		mustDo(schedulerKind, numberOfRotations,numberOfTranslations),
		_showPlan(false)
	{
		init(0);
	}
	void init(int iter) {
		_showPlan = showPlan(iter);
		capturedDataStats.init();
		mustDo.init();
	}
	void fini() {
		mustDo.fini();
		capturedDataStats.fini();
	}
	void append(
		int rotation,
		int translation) {
		mustDo.append(rotation, translation);
	}

	class Plan {
		// Plans doing all the samples of all the rotations and translations for one tile.
		// Executes one block at a time.
	public:
		static const S64 maxRotStep   = 64;
		static const S64 maxTransStep = 64;

		Scheduler & parent;
		const S64 tileRotationBegin;    const S64 tileRotationEnd;	  const S64 tileRotationStep;
		const S64 tileTranslationBegin; const S64 tileTranslationEnd; const S64 tileTranslationStep;
															          const S64 sampleStep;
		Plan(
			Scheduler & parent,
			S64 tileRotationBegin,    S64 tileRotationEnd, 
			S64 tileTranslationBegin, S64 tileTranslationEnd, 
			S64 sampleStep,
			bool showPlan)
		  : parent(parent),
			tileRotationBegin   (tileRotationBegin),    tileRotationEnd   (tileRotationEnd),    tileRotationStep   (tileRotationEnd-tileRotationBegin),
			tileTranslationBegin(tileTranslationBegin), tileTranslationEnd(tileTranslationEnd), tileTranslationStep(tileTranslationEnd-tileTranslationBegin),
																			                    sampleStep         (sampleStep),
			translationIndexsCapacity(parent.translationIndexsCapacity),
			showPlan(showPlan)
		{
			if (tileRotationStep > 64 || tileTranslationStep > 64) {
				std::cerr << "GetAllSquaredDifferences_Kernel_implementation::Plan has a limit of 64 on tileRotationStep:" << tileRotationStep << " and tileTranslationStep:" << tileTranslationStep << std::endl;
				std::cout << "GetAllSquaredDifferences_Kernel_implementation::Plan has a limit of 64 on tileRotationStep:" << tileRotationStep << " and tileTranslationStep:" << tileTranslationStep << std::endl;
				assert(!"too many tileRotationStep or tileTranslationStep");
				EXIT_ABNORMALLY;
			}
			stepsCapacity = tileRotationStep*tileTranslationStep + 1;
			stepsPool.allocate(stepsCapacity, singletonSteps, mergedSteps);
			prepare();
		}

		~Plan() {
			stepsPool.deallocate(singletonSteps, mergedSteps);
		}

	private:
		const bool showPlan;
		const int translationIndexsCapacity;

		U64	  stepsCapacity;

		S64				mergedStepsLen;
		MergedStep*		mergedSteps;

		S64				singletonStepsLen;
		SingletonStep*	singletonSteps;

		void prepare() {
			prepareSingletons();
			prepareAnyRemaining();
		}

		void prepareSingletons() {
			// Finds all the cases where neither the r nor the t is used again
			// so that they should be done together using different optimizations
			//		No benefit to tiling, except the weights
			//		No benefit to converting to float ahead of time
			//		No benefit to sqrt'ing the weight
			//
			singletonStepsLen = 0;
		}

		void prepareAnyRemaining() {
			const bool logging = false;

			// Build the rTransPending table that shows for each rotation, which translations have to be done
			//
			if (logging) {
				std::cout << "Scanning"
				    << " tileRotationBegin:" << tileRotationBegin << "..tileRotationEnd:" << tileRotationEnd 
				    << " tileTranslationBegin:" << tileTranslationBegin << "..tileTranslationEnd:" << tileTranslationEnd 
					<< std::endl;
			}

			static_assert(sizeof(U64)*8 == maxRotStep,   "Using a bitmask for the rots");
			static_assert(sizeof(U64)*8 == maxTransStep, "Using a bitmask for the trans");
			U64 rTransPending[maxRotStep];
			U64 rTransPendingCount;

			parent.mustDo.getSlice(
				rTransPending,
				rTransPendingCount,
				tileRotationBegin,
				tileRotationEnd,
				tileTranslationBegin,
				tileTranslationEnd);

			S64 rTransPendingFirstMaybeNonZero = 0;

			// If we try to keep the whole tile in L2, so "inCache" refers to in L1 cache
			// We further tile into the L1 cache and hence set the age by the number of distinct rotations that can be held in it
			//
			// L1 cache is 32KB
			//
			const U64 floatsNeededForOneITransOrIRotOrMinvsigma2 = sampleStep*sizeof(float);							// typically 256*4 = 1K
			const U64 maxAge = 
				( 24000												// L1 size we are willing to use						// 24K
				- 3*floatsNeededForOneITransOrIRotOrMinvsigma2		// Needed for the rot real and imag, and the Minvsigma2	// -3K
				) / (2*floatsNeededForOneITransOrIRotOrMinvsigma2)	// Needed for the trans real and image					// /2K
				*2/3;												// But don't fill the cache								// So will be about 6
			if (logging) std::cout << "maxAge:" << maxAge << std::endl;
			
			U64 clock(100);

			typedef U64 InCacheSince[maxTransStep];
			auto addToCache = [&](InCacheSince & inCacheSince, U64 & inCache, U64 maskToAdd, const char* cache) {
				while (maskToAdd) {
					U64 x;
					bitScanForward(x, maskToAdd);
					auto m = U64(1) << x;
					maskToAdd ^= m;
					inCacheSince[x] = clock;
					inCache |= m;
					if (logging) std::cout << "Added " << cache << ":" << x << " at clock:" << clock << std::endl;
				}
			};
			auto updateInCache = [&](InCacheSince const & inCacheSince, U64 & inCache, const char* cache) {
				auto remaining = inCache;
				while (remaining) {
					U64 x;
					bitScanForward(x, remaining);
					auto m = U64(1) << x;
					remaining ^= m;
					if (inCacheSince[x] + maxAge >= clock) return;
					inCache ^= m;
					if (logging) std::cout << "Aged out " << cache << ":" << x << ", entered at clock:" << inCacheSince[x] << std::endl;
				}
			};

			InCacheSince transInCacheSince;
			InCacheSince rotInCacheSince;
			U64 transInCache(0);
			U64 rotInCache(0);

			auto addToTransInCache = [&](U64 maskToAdd) {
				addToCache(transInCacheSince, transInCache, maskToAdd, "t");
			};

			auto addToRotInCache = [&](U64 maskToAdd) {
				addToCache(rotInCacheSince, rotInCache, maskToAdd, "r");
			};

			auto updateTransInCache = [&]() {
				updateInCache(transInCacheSince, transInCache, "t");
			};

			auto updateRotInCache = [&]() {
				updateInCache(rotInCacheSince, rotInCache, "r");
			};

			auto appendMergedSteps = [&]() {
				if (mergedStepsLen > 0) {
					assert(mergedStepsLen-1 < stepsCapacity);
					auto & doneTranslationIndexsLen = mergedSteps[mergedStepsLen - 1].translationIndexsLen;
					auto & doneRotationLen          = mergedSteps[mergedStepsLen - 1].rotationLen;
					if (doneTranslationIndexsLen == 0) return;
					if (doneRotationLen != 1) {
						assert(doneTranslationIndexsLen == 3);
					}
				}
				assert(mergedStepsLen < stepsCapacity);
				mergedSteps[mergedStepsLen].rotationBegin        = 0;
				mergedSteps[mergedStepsLen].rotationLen          = 0;
				mergedSteps[mergedStepsLen].translationIndexsLen = 0;
				mergedStepsLen++;
			};

			// Start the first merge
			mergedStepsLen    = 0;
			singletonStepsLen = 0;
			appendMergedSteps();

			// Decide on the next one to do, and append it to the mergedSteps, until all done
			//
			for (;;) {
				if (logging) std::cout << "rTransPendingCount:" << rTransPendingCount << std::endl;

				if (rTransPendingCount == 0) break;

				// Find a pending that has both or either of its needs in the cache
				// Prefer ones with a common irot to what we are currently merging
				//		Try the current r for a t in the cache
				//		Try the current r for any t	- the unrolling is VERY valuable
				//		Pick the r that has the most t's in the cache
				//		Pick the r that has the most t's, so that we get a good chance at merging
				//
				char const* why = "no reason";
				U64 t;	U64 tMask = 0;
				S64 r;	S64 rLen = 1;
				
				{
					assert(mergedStepsLen < stepsCapacity);
					auto const mergedStep			= mergedSteps[mergedStepsLen-1];
					auto const translationIndexsLen = mergedStep.translationIndexsLen;

					if (translationIndexsLen > 0) {
						r = mergedStep.rotationBegin - tileRotationBegin;
						assert(0 <= r && r < tileRotationStep);

						why = "r same as last, t in cache";
						auto tInCache = rTransPending[r] & transInCache;
						if (bitScanForward(t, tInCache)) {
							goto DoIt;
						}

						why = "r same as last, t pending";
						if (bitScanForward(t, rTransPending[r])) {
							goto DoIt;
						}
					} else if (use3x3 && parent.translationIndexsCapacity >= 3) {
						// If {r1, r2, .. rR} and {t1, t2, .. tT} require all (ri, tj) to be done
						// it can be done by reading w, R real and imag, T real and imag and combining 
						// so the formula for the reads per pair done is   (1+2R+2T)/(RT)
						// which is symmetric so has a minimum when R==T.  
						// The above code aims for 1x5 which gives 13/5  = 2.6		memory reads per cell
						// 2x2 would be                             9/4  = 2.25
						// 3x3 would be                            13/9  = 1.4, so may be worth trying for, even though the nine combos must be explicitly coded in the execute
						// 4x4 would be                            17/16 = 1.1, but sure to be compute bound

						S64 bestPopCount = 0;
						S64 bestI  = 0;
						U64 bestMT = 0;

						for (S64 i = rTransPendingFirstMaybeNonZero; 
								i < std::min(tileRotationStep-2, rTransPendingFirstMaybeNonZero + 10);	// Don't search too far
								i++) {
							auto mt  = rTransPending[i+0];
								 mt &= rTransPending[i+1];
								 mt &= rTransPending[i+2];
							S64 popCount = POPCNT64(mt);
							if (popCount < 3) continue;
							auto mtc = mt & transInCache;
							auto mrc = ((U64(1)<<(i+2))-(U64(1)<<i)) & rotInCache;
							popCount += POPCNT64(mtc) + POPCNT64(mrc);
							// std::cout << "FOUND A 3x3 with " << POPCNT64(mtc) << " trans in cache, " << POPCNT64(mrc) << " rot in cache, totalling " << popCount << std::endl;
							if (bestPopCount >= popCount) continue;
							bestPopCount = popCount;
							bestI = i;
							bestMT = mt;
						}

						if (bestMT != 0) {
							
							auto bestToKeep   = bestMT &  transInCache;
							auto othersToKeep = bestMT & ~transInCache;
							bestMT = 0;
							for (int i = 0; i < 3; i++) {
								U64 x;
								if (!bitScanForward(x, bestToKeep)) {
									bestToKeep = othersToKeep;
									bitScanForward(x, bestToKeep);
								}
								U64 keep = U64(1) << x;
								bestMT     ^= keep;
								bestToKeep ^= keep;
							}

							r     = bestI;
							rLen  = 3;
							tMask = bestMT;
							bitScanForward(t, tMask);

							goto DoIt;
						}
					}

					{
						why = "r has the most t pending in the cache";
						S64 bestPopCount = 0;
						for (S64 i = rTransPendingFirstMaybeNonZero; 
								i < std::min(tileRotationStep, rTransPendingFirstMaybeNonZero + 10);	// Don't search too far
								i++) {
				    		auto m = rTransPending[i] & transInCache;
				    		S64 popCount = POPCNT64(m);
				    		if (bestPopCount < popCount) {
				    			bestPopCount = popCount;
				    			r = i;
				    			bitScanForward(t, m);
				    		}
						}
				    
						if (bestPopCount) goto DoIt;
				    
						// std::cout << "DEBUG" << std::endl;
				    
						why = "r has the most t pending";
						for (S64 i = rTransPendingFirstMaybeNonZero; i < tileRotationStep; i++) {
				    		auto m = rTransPending[i];
				    		S64 popCount = POPCNT64(m);
				    		if (bestPopCount < popCount) {
				    			bestPopCount = popCount;
				    			r = i;
				    			bitScanForward(t, m);
				    		}
						}
						assert(bestPopCount > 0);
					}
				}

			DoIt:
				// Do r,t, but merge ones that have a common irot
				// TODO - a similar game with a common itrans
				//
				assert(0 <= t && t      <  tileTranslationStep);
				assert(0 <= r && r+rLen <= tileRotationStep);

				auto rMask = (U64(1)<<U64(r+rLen)) - (U64(1)<<U64(r));

				if (tMask == 0) tMask = U64(1L)<<U64(t);

				for (int rx = r; rx < r + rLen; rx++) {
					assert((rTransPending[rx] & tMask) == tMask);
					rTransPending[rx]   ^= tMask;
					rTransPendingCount -= POPCNT64(tMask);

					if (rx == rTransPendingFirstMaybeNonZero && rTransPending[rx] == 0 && rTransPendingCount > 0) {
						while (rTransPending[rTransPendingFirstMaybeNonZero] == 0) {
							rTransPendingFirstMaybeNonZero++;
						}
						assert(rTransPendingFirstMaybeNonZero < tileRotationStep);
					}
				}

				// Age items in the cache
				//
				auto oldTransInCache = transInCache;
				auto oldRotInCache   = rotInCache;
				transInCache |= tMask;
				rotInCache   |= rMask;

				addToTransInCache(tMask);
				addToRotInCache  (rMask);

				auto oldClock = clock;
				if (oldTransInCache != transInCache) clock += POPCNT64(tMask);
				if (oldRotInCache   != rotInCache)   clock += POPCNT64(rMask);

				if (oldClock != clock) {
					if (logging) std::cout << "Clock:" << clock << std::endl;
					updateTransInCache();
					updateRotInCache();
				}

				// Do this one
				//
				if (rLen == 1) {

					int rotation    = r + tileRotationBegin;
					int translation = t + tileTranslationBegin;

					// There are some singletons that can be done more efficiently using the original data
					// because there is not enough work to justify capturing them
					//
					if (parent.mustDo.rotationUseCount(rotation)       <= 1			// BEVIN - WHAT IS THE RIGHT NUMBER?
					||  parent.mustDo.translationUseCount(translation) <= 1) {

						parent.capturedDataStats.operate(false, rLen*1);

						// Shift these ones to the singletonSteps
						//
						assert(singletonStepsLen < stepsCapacity);
						singletonSteps[singletonStepsLen].rotation    = rotation;
						singletonSteps[singletonStepsLen].translation = translation;
						singletonStepsLen++;

					} else {

						parent.capturedDataStats.operate(true, rLen*1);

						if (mergedSteps[mergedStepsLen-1].rotationBegin != rotation) {
							appendMergedSteps();
							mergedSteps[mergedStepsLen-1].rotationBegin = rotation;
						}
						mergedSteps[mergedStepsLen-1].rotationLen = 1;

						if (logging) std::cout << (mergedSteps[mergedStepsLen-1].translationIndexsLen ? "Merging" : "Starting")
								<< " r:" << r << " t:" << t << " because " << why << std::endl;

						auto & mergedStep = mergedSteps[mergedStepsLen-1];
						assert(mergedStep.translationIndexsLen < translationIndexsCapacity);
						mergedStep.translationIndexs[mergedStep.translationIndexsLen++] = translation;

						if (mergedStep.translationIndexsLen == translationIndexsCapacity) appendMergedSteps();
					}

				} else {
					auto & step = mergedSteps[mergedStepsLen-1];
					
					assert(use3x3);
					assert(step.translationIndexsLen == 0);
					assert(rLen == 3);
					assert(POPCNT64(tMask) == 3);
					parent.capturedDataStats.operate(true, 9);

					assert(3 <= parent.translationIndexsCapacity);
					step.rotationBegin = r + tileRotationBegin;
					step.rotationLen   = 3;

					while (tMask) {
						U64 x;
						bitScanForward(x, tMask);
						assert(step.translationIndexsLen < translationIndexsCapacity);
						step.translationIndexs[step.translationIndexsLen++] = x + tileTranslationBegin;
						tMask ^= U64(1) << x;
					}

					appendMergedSteps();
				}
			}

			appendMergedSteps();

			if (showPlan) {
				auto dashes = [&](){ std::cout << "    "; for (int t = 0; t < tileTranslationStep; t++) { std::cout << '_'; }  std::cout << std::endl; };
				std::cout << "    Horizontal is " << tileTranslationStep << " translations, vertical is " << tileRotationStep << " rotations" << std::endl;
				U64 toDo(0);					// BEVIN Does not account for singletons
				int  pendingBlankLines(0);		// Strips off leading and trailing blank lines
				bool nonBlankOutputLines(false);
				for (S64 rotation = tileRotationBegin; rotation < tileRotationEnd; rotation++) {
					int pendingSpaces(4);
					bool nonBlankOutput(false);
					for (S64 translation = tileTranslationBegin; translation < tileTranslationEnd; translation++) {
						auto isSet = parent.mustDo.isSet(rotation, translation);
						if (!isSet) pendingSpaces++;
						else {
							if (toDo == 0) dashes();
							for (;pendingBlankLines>0; pendingBlankLines--) std::cout << std::endl;
							for (;pendingSpaces>0; pendingSpaces--) std::cout << ' ';
							std::cout << '*';
							nonBlankOutput = true;
							toDo++;
						}
					}
					if (nonBlankOutput) {
						std::cout << std::endl;
						nonBlankOutputLines = true;
					} else if (nonBlankOutputLines) {
						pendingBlankLines++;
					}
				}
				if (!toDo) {
					std::cout << "    Nothing to do" << std::endl;
				} else {
					std::cout << "    Number of stars: " << toDo << std::endl;
					dashes();
				}

				U64 done(0);

				U64 counts[maxTranslationIndexsCapacity+1];
				for (int i = 0; i < translationIndexsCapacity+1; i++) counts[i] = 0;
				for (int i = 0; i < mergedStepsLen; i++) {
					auto& step = mergedSteps[i];
					auto x = step.translationIndexsLen;
					assert(x <= translationIndexsCapacity);
					counts[x]++;
					done += step.rotationLen*step.translationIndexsLen;
				}
				std::cout << "translationIndexsLen histogram";
				for (int i = 1; i < translationIndexsCapacity+1; i++) {
					if (counts[i]) std::cout << "  i:" << i << " " << counts[i];
				}
				std::cout << std::endl;

				done += singletonStepsLen;

				if (toDo != done) {
					std::cerr << "toDo != done by steps" << std::endl;
					for (U64 i = 0; i < mergedStepsLen; i++) {
						std::cout << i; mergedSteps[i].print(" is");
					}
				}
			}
		}

	public:
		// The plan is executed for each brick
		//
		void executeMergedSteps(int sLo, int sHi) {
			for (int i = 0; i < mergedStepsLen; i++) {
				parent.executeMergedStep(mergedSteps[i], *this, sLo, sHi);
			}
		}
		void executeSingletons() {
			if (singletonStepsLen > 0) parent.executeSingletons(singletonSteps, singletonStepsLen);
		}
	};  // class Plan

	virtual void executeSingletons(SingletonStep const * steps, int numberOfSteps) = 0;
	virtual void executeMergedStep(MergedStep const & step, Plan& plan, int sLo, int sHi) = 0;

	void compute(const S64 rotationStepProposed, const S64 translationStepProposed, const S64 sampleStep) {
		TUNING_SCOPE_STEP_BEGIN(kernel::compute)
		const S64 rotationStep    = std::min(rotationStepProposed,    Plan::maxRotStep);
		const S64 translationStep = std::min(translationStepProposed, Plan::maxTransStep);

		if (_showPlan) { 
			std::cout << "===========================================================================" << std::endl;
			std::cout << toStr(schedulerKind) << "::compute" << std::endl;
		}

		for (int translationBegin = 0; translationBegin < numberOfTranslations; translationBegin += translationStep) {
			auto translationEnd = std::min(numberOfTranslations, translationBegin + translationStep);

			for (int rotationBegin = 0; rotationBegin < numberOfRotations; rotationBegin += rotationStep) {
				auto rotationEnd = std::min(numberOfRotations, rotationBegin + rotationStep);

				auto sEndBeforePlan = numberOfAcceptedSamples();

				if (_showPlan) {	// Must be before Plan constructed
					std::cout << "Scheduler.plan for "
						<< " r:" << rotationBegin    << ".." << rotationEnd 
						<< " t:" << translationBegin << ".." << translationEnd 
						<< " s:0.." << sEndBeforePlan
						<< std::endl;
				}

				Plan plan(*this, rotationBegin, rotationEnd, translationBegin, translationEnd, sampleStep, _showPlan);

				auto sEnd = numberOfAcceptedSamples();
				if (sEndBeforePlan != sEnd) { 
					std::cerr << "*** Plan changes sEnd ***" << std::endl;
					EXIT_ABNORMALLY; 
				}

				for (int sLo = 0; sLo < sEnd;      sLo += sampleStep) {
					auto sHi = std::min<int>(sEnd, sLo +  sampleStep);
					plan.executeMergedSteps(sLo, sHi);
		        }

				plan.executeSingletons();
			}
		}
		if (_showPlan) { 
			std::cout << "Scheduler.capturedDataStatsHistogram : ";
			capturedDataStats.record();
			capturedDataStatsHistogram.display(); 
			std::cout << std::endl;
			capturedDataStatsHistogram.init();
		}
		TUNING_SCOPE_STEP_END
	}

	virtual int numberOfAcceptedSamples() const { return numberOfSamples; }

private:
	class MustDo {
		static S64 computeLen(const int numberOfRotations, const int numberOfTranslations) {
			S64 needed = numberOfRotations*numberOfTranslations;
			needed = (needed+63)/64;
			return needed;
		}
	public:
		SchedulerKind const schedulerKind;
		MustDo(SchedulerKind schedulerKind, const int numberOfRotations, const int numberOfTranslations)
		  : schedulerKind(schedulerKind),
		    _numberOfRotations(numberOfRotations), _numberOfTranslations(numberOfTranslations), len(computeLen(numberOfRotations, numberOfTranslations)),
		    _prevRotation(-1), _prevTranslation(0), _rotationsWithTranslations(0), _appends(0), _skips(0) { 
			mustDoPools[schedulerKind].allocate(_rotationUseCount, _translationUseCount, bits, numberOfRotations, numberOfTranslations, len);
			init();
		}
		void init() {
			for (S64 i = 0; i < _numberOfRotations;    i++) _rotationUseCount   [i] = 0;
			for (S64 i = 0; i < _numberOfTranslations; i++) _translationUseCount[i] = 0;
			for (S64 i = 0; i < len;				   i++) bits                [i] = 0;
		}
		void fini() {}
		~MustDo() {
			if (false && _rotationsWithTranslations > 0) {
				std::cerr << "~MustDo appends:" << _appends 
					<< " skips:" << _skips 
					<< " skips/_rotationsWithTranslations:" << std::setprecision(2) << float(_skips)/float(_rotationsWithTranslations)
					<< std::endl;
			}
			mustDoPools[schedulerKind].deallocate(_rotationUseCount, _translationUseCount, bits);
		}
		void append(int rotation, int translation) {
			S64 index; U64 mask; U64 maskOffset;
			indexAndMask(index, mask, maskOffset, rotation, translation);
			assert(!(bits[index]&mask));
			bits[index] |= mask;
			_rotationUseCount   [rotation   ]++;
			_translationUseCount[translation]++;
			_appends++;
			if (_prevRotation != rotation) {
				_rotationsWithTranslations++;
				_prevRotation = rotation;
			} else if (_prevTranslation+1 != translation) {
				_skips++;
			}
			_prevTranslation = translation; 
		}
		bool isSet(int rotation, int translation) const {
			S64 index; U64 mask; U64 maskOffset;
			indexAndMask(index, mask, maskOffset, rotation, translation);
			return (bits[index]&mask) != 0;
		}
		int rotationUseCount   (int rotation)    const { assert(rotation    < _numberOfRotations   ); return _rotationUseCount   [rotation];    }
		int translationUseCount(int translation) const { assert(translation < _numberOfTranslations); return _translationUseCount[translation]; }

		void getSlice(
			U64 (& rTransPending)[Plan::maxRotStep],
			U64 & rTransPendingCount,
			S64 tileRotationBegin,
			S64 tileRotationEnd,
			S64 tileTranslationBegin,
			S64 tileTranslationEnd) {

			// The simple way of implementing this operation is surprisingly expensive
			//
			static const bool testing = false;
			if (testing) {
				rTransPendingCount = 0;
				for (S64 r = 0; r < Plan::maxRotStep; r++) rTransPending[r] = 0;

				for (S64 rotation = tileRotationBegin; rotation < tileRotationEnd; rotation++) {
					for (S64 translation = tileTranslationBegin; translation < tileTranslationEnd; translation++) {
						if (!isSet(rotation, translation)) continue;
						auto m = U64(1) << (translation-tileTranslationBegin);
						rTransPending[rotation-tileRotationBegin] |= m;
						rTransPendingCount++;
					}
				}
			}

			// Given that the operation is a bit substring copy, and there is at most maxTransStep bits, it is this...
			//
		DoAgain:
			U64 new_rTransPendingCount = 0;
			static_assert(Plan::maxTransStep <= 64, "Too large maxTransStep");
			S64 numberOfTranslations = tileTranslationEnd - tileTranslationBegin;
			assert(numberOfTranslations <= 64);
			for (S64 rotation = tileRotationBegin; rotation < tileRotationEnd; rotation++) {
				S64 index;
				U64 mask;
				U64 maskOffset;
				indexAndMask(index, mask, maskOffset, rotation, tileTranslationBegin);
				// Get these bits out, whether or not they span the boundary 
				//						index+1																			index
				//	[............................................................####][######################################..........................]
				//																^maskOffset+numberOfTranslations			^maskOffset
				//
				// Do the two cases separately
				//
				U64 m0 = bits[index+0];

				S64 endOffset = numberOfTranslations + maskOffset;
				if (endOffset <= 64) {
					// Doesn't span boundary
					m0        <<= 64 - endOffset;
					maskOffset += 64 - endOffset;
				} else {
					// Spans boundary
					U64 m1 = bits[index+1] << (128 - endOffset);
					m0        >>= endOffset - 64;
					m0 |= m1;
					maskOffset -= endOffset - 64;
				}
				m0 >>= maskOffset;
				if (!testing) {
					rTransPending[rotation - tileRotationBegin] = m0;
				} else {
					if (rTransPending[rotation - tileRotationBegin] != m0) {
						std::cerr << "Wrong mask" << std::endl;
						goto DoAgain;
					}
				}
				new_rTransPendingCount += POPCNT64(m0);
			}
			if (!testing) {
				rTransPendingCount = new_rTransPendingCount;
			} else {
				if (rTransPendingCount != new_rTransPendingCount) {
					std::cerr << "Wrong new_rTransPendingCount" << std::endl;
					goto DoAgain;
				}
			}
		}

	private:
		int  _numberOfRotations;
		int  _numberOfTranslations;
		int  _prevRotation, _prevTranslation, _rotationsWithTranslations, _appends, _skips;
		int* _rotationUseCount;
		int* _translationUseCount;
		S64  len;
		U64* bits;
		void indexAndMask(S64 & index, U64 & mask, U64 & maskOffset, int rotation, int translation) const {
			S64 x = rotation*_numberOfTranslations + translation;
			assert(0 <= x);
			index = x/64;
			assert(index < len);
			maskOffset = (x % 64);
			mask = U64(1) << maskOffset;
		}  
	} mustDo;
public:
	bool mustDoIsSet(int rotation, int translation) const { return mustDo.isSet(rotation, translation); }
};  // class Scheduler

const S64 Scheduler::Plan::maxRotStep;
const S64 Scheduler::Plan::maxTransStep;
	// Linux needs these

class BackProjection_Kernel_implementation : public BackProjection_Kernel, Scheduler {
public:
	typedef BackProjection_Kernel_implementation Kernel;

	static BackProjection_Kernel_implementation* make(
		const int translationIndexsCapacity,
		const int numberOfRotations,
		const int numberOfTranslations,
		const int numberOfSamples) {
#include "./util_heap_undefs.h"
		return sNewA(BackProjection_Kernel_implementation, (
			translationIndexsCapacity,
			numberOfRotations,
			numberOfTranslations,
			numberOfSamples));
#include "./util_heap_defs.h"
	}

	BackProjection_Kernel_implementation(
		const int translationIndexsCapacity,
		const int numberOfRotations,
		const int numberOfTranslations,
		const int numberOfSamples) 
	  :
	    Scheduler(
			SchedulerKind_BackProjection,
#ifdef MY_MACHINE_HAS_AVX512
			translationIndexsCapacity,
#else
			std::min(4, translationIndexsCapacity),
				// the 5x1 uses too many vector registers
#endif
			numberOfRotations, numberOfTranslations, numberOfSamples),
		rotationState   (NULL),
		iclass          (-1),
		rotations_real  (numberOfRotations),
		rotations_imag  (numberOfRotations),
		rotations_weight(numberOfRotations),
		translations_real(numberOfTranslations),
		translations_imag(numberOfTranslations),
		weights(SchedulerKind_BackProjection, numberOfRotations,numberOfTranslations)
	{
	}

	bool suitable(
		const int numberOfRotations,
		const int numberOfTranslations,
		const int numberOfSamples) {

		return true
			&& this->numberOfRotations	  == numberOfRotations
			&& this->numberOfTranslations == numberOfTranslations
			&& this->numberOfSamples	  == numberOfSamples;
	}

	void init(
		int iter,
		RotationState& rotationState,
		const int iclass) {

		Scheduler::init(iter);

		this->rotationState = &rotationState;
		this->iclass        = iclass;

		for (int i = 0; i < numberOfRotations; i++) {
			rotations_real  [i].init();
			rotations_imag  [i].init();
			rotations_weight[i].init();
		}
		for (int i = 0; i < numberOfTranslations; i++) {
			translations_real  [i].init();
			translations_imag  [i].init();
		}
		ctf.init();
		weights.init();
	}

	virtual void release();

	virtual ~BackProjection_Kernel_implementation() {}

	RotationState*   		  rotationState;
	int						  iclass;
	std::vector<CapturedData> rotations_real;
	std::vector<CapturedData> rotations_imag;
	std::vector<CapturedData> rotations_weight;
	std::vector<CapturedData> translations_real;
	std::vector<CapturedData> translations_imag;
	CapturedData              ctf;
	FloatsForRT				  weights;

	virtual void append(
		int rotation,
		int translation,
		      double* rotation_real,          double* rotation_imag,		// where to add back into, in a critical region
		const double* translation_real, const double* translation_imag,
		      double* rotation_weight,        const double* ctf,
			  double  weight) {
		Scheduler::append(rotation,translation);
		rotations_real   [rotation   ].captureForSummingInto(&capturedDataStats, "bpk_rotations_real",    numberOfSamples, rotation_real);
		rotations_imag   [rotation   ].captureForSummingInto(&capturedDataStats, "bpk_rotations_imag",    numberOfSamples, rotation_imag);
		rotations_weight [rotation   ].captureForSummingInto(&capturedDataStats, "bpk_rotations_weight",  numberOfSamples, rotation_weight);
		translations_real[translation].capture              (&capturedDataStats, "bpk_translations_real", numberOfSamples, translation_real);
		translations_imag[translation].capture              (&capturedDataStats, "bpk_translations_imag", numberOfSamples, translation_imag);
		this->ctf                     .capture              (&capturedDataStats, "bpk_ctf",               numberOfSamples, ctf);

		weights[weights.index(rotation, translation)] = weight;

		rotationState->setRotationState(iclass, rotation, rotationState->RS_captured, 0, rotation_real);
	}

	virtual void executeSingletons(SingletonStep const * steps, int numberOfSteps) {
		auto ctf = this->ctf.ptrValidUpto(numberOfSamples);

		for (int step = 0; step < numberOfSteps;) {
			static const int MaxNumberOfStepsExecuted = 3;
			int			     numberOfStepsExecuted    = std::min(MaxNumberOfStepsExecuted, numberOfSteps - step);

#define HEADSTEP(SD) \
			const auto rotation##SD    = steps[step+SD].rotation;										\
			const auto translation##SD = steps[step+SD].translation;									\
			const auto r_real##SD      = rotations_real  [rotation##SD].ptrValidUpto(numberOfSamples);	\
			const auto r_imag##SD      = rotations_imag  [rotation##SD].ptrValidUpto(numberOfSamples);	\
   				  auto rw##SD          = rotations_weight[rotation##SD].ptrValidUpto(numberOfSamples);	\
			assert(!!r_real##SD);																		\
			assert(!!r_imag##SD);																		\
			assert(U64(r_real##SD) % 64 == 0);															\
			assert(U64(r_imag##SD) % 64 == 0);															\
			auto t_real##SD        = translations_real[translation##SD].ptrValidUpto(numberOfSamples);	\
			auto t_imag##SD        = translations_imag[translation##SD].ptrValidUpto(numberOfSamples);	\
			assert(!!t_real##SD);																		\
			assert(!!t_imag##SD);																		\
			assert(U64(t_real##SD) % 64 == 0);															\
			assert(U64(t_imag##SD) % 64 == 0);															\
			auto w##SD = weights[weights.index(rotation##SD, translation##SD)];							// end of macro
			
		#define BODYSTEP(SD) \
			{   r_real##SD[n] += t_real##SD[n] * w##SD;													\
				r_imag##SD[n] += t_imag##SD[n] * w##SD;													\
				rw##SD    [n] += ctf       [n] * w##SD;													\
			}																							// end of macro

			switch (numberOfStepsExecuted) {
			case 1: {
				HEADSTEP(0)
				#pragma vector aligned
    		    #pragma ivdep
				for (int n = 0; n < numberOfSamples; n++) {
					BODYSTEP(0)
				}
			}	break;
			case 2: {
				HEADSTEP(0) HEADSTEP(1)
				#pragma vector aligned
    		    #pragma ivdep
				for (int n = 0; n < numberOfSamples; n++) {
					BODYSTEP(0) BODYSTEP(1) 
				}
			}	break;
			case 3: {
				HEADSTEP(0) HEADSTEP(1) HEADSTEP(2)
				#pragma vector aligned
    		    #pragma ivdep
				for (int n = 0; n < numberOfSamples; n++) {
					BODYSTEP(0) BODYSTEP(1) BODYSTEP(2)
				}
			}	break;
			default:
				std::cerr << "BackProjection_Kernel::executeSingletons: numberOfStepsExecuted too large" << std::endl;
				EXIT_ABNORMALLY;
			}

#undef BODYSTEP
#undef HEADSTEP

			step += numberOfStepsExecuted;
		}
	}

	virtual void executeMergedStep(MergedStep const & step, Plan& plan, int sLo, int sHi) {
		auto  translationIndexsLen = step.translationIndexsLen;
		auto  rotationBegin        = step.rotationBegin;
		auto  rotationLen		   = step.rotationLen;
		auto& translationIndexs    = step.translationIndexs;
		
		if (translationIndexsLen == 0) return;

		const bool logging = false && rotationBegin==0 && translationIndexs[0]==0;

		auto r_real0 = rotations_real  [rotationBegin + 0].ptrValidUpto(sHi);
		auto r_imag0 = rotations_imag  [rotationBegin + 0].ptrValidUpto(sHi); 
   	    auto rw0	 = rotations_weight[rotationBegin + 0].ptrValidUpto(sHi);
		auto ctf     = this->ctf						  .ptrValidUpto(sHi);
		assert(U64(r_real0) % 64 == 0); 
		assert(U64(r_imag0) % 64 == 0); 
		assert(U64(rw0    ) % 64 == 0); 

		#define HEADSTEP_SHARED(TD) \
			auto t_real##TD = translations_real[translationIndexs[TD]].ptrValidUpto(sHi);		\
			auto t_imag##TD = translations_imag[translationIndexs[TD]].ptrValidUpto(sHi);		\
			assert(!!t_real##TD);																\
			assert(!!t_imag##TD);																\
			assert(U64(t_real##TD) % 64 == 0); 													\
			assert(U64(t_imag##TD) % 64 == 0); 													// end of macro

		#define HEADSTEP_R(RD, TD) \
			auto w##RD##TD = weights[weights.index(rotationBegin+RD, translationIndexs[TD])];

		#define HEADSTEP(TD) \
			HEADSTEP_SHARED(TD) HEADSTEP_R(0, TD)
				
		#define BODYSTEP_R(RD, TD) \
			{   r_real##RD[n] += t_real##TD[n] * w##RD##TD;						\
				r_imag##RD[n] += t_imag##TD[n] * w##RD##TD;						\
				rw##RD    [n] += ctf       [n] * w##RD##TD;						\
			}																	// end of macro

		#define BODYSTEP(TD) \
			BODYSTEP_R(0,TD)

		if (rotationLen == 3) {
			// 3x3 is done as simply three repeats of a 1x3
			// It is the only case known how to do
			//
			assert(rotationLen  == 3);
			assert(translationIndexsLen == 3);

			auto r_real1 = rotations_real  [rotationBegin + 1].ptrValidUpto(sHi);
			auto r_imag1 = rotations_imag  [rotationBegin + 1].ptrValidUpto(sHi); 
			auto r_real2 = rotations_real  [rotationBegin + 2].ptrValidUpto(sHi);
			auto r_imag2 = rotations_imag  [rotationBegin + 2].ptrValidUpto(sHi); 
   			auto rw1	 = rotations_weight[rotationBegin + 1].ptrValidUpto(sHi); 
   			auto rw2	 = rotations_weight[rotationBegin + 2].ptrValidUpto(sHi); 
			assert(U64(r_real1) % 64 == 0); 
			assert(U64(r_imag1) % 64 == 0); 
			assert(U64(r_real2) % 64 == 0); 
			assert(U64(r_imag2) % 64 == 0); 
			assert(U64(rw1    ) % 64 == 0); 
			assert(U64(rw2    ) % 64 == 0); 

			HEADSTEP_SHARED(0) HEADSTEP_R(0,0)  HEADSTEP_R(1,0)  HEADSTEP_R(2,0)
			HEADSTEP_SHARED(1) HEADSTEP_R(0,1)  HEADSTEP_R(1,1)  HEADSTEP_R(2,1)
			HEADSTEP_SHARED(2) HEADSTEP_R(0,2)  HEADSTEP_R(1,2)  HEADSTEP_R(2,2)

			#pragma vector aligned
    	    #pragma ivdep
			for (int n = sLo; n < sHi; n++) {
				// BEVIN - TODO - see if there is enough registers to do this
				// The compiler did not generate good code for the macros
				// BODYSTEP_R(0,0)  BODYSTEP_R(1,0)  BODYSTEP_R(2,0)
				// BODYSTEP_R(0,1)  BODYSTEP_R(1,1)  BODYSTEP_R(2,1)
				// BODYSTEP_R(0,2)  BODYSTEP_R(1,2)  BODYSTEP_R(2,2)
				r_real0[n] += t_real0[n] * w00 + t_real1[n] * w01 + t_real2[n] * w02; 
				r_imag0[n] += t_imag0[n] * w00 + t_imag1[n] * w01 + t_imag2[n] * w02; 
				    rw0[n] += ctf    [n] * w00 + ctf    [n] * w01 + ctf    [n] * w02;

				r_real1[n] += t_real1[n] * w10 + t_real1[n] * w11 + t_real2[n] * w12; 
				r_imag1[n] += t_imag1[n] * w10 + t_imag1[n] * w11 + t_imag2[n] * w12; 
				    rw1[n] += ctf    [n] * w10 + ctf    [n] * w11 + ctf    [n] * w12;

				r_real2[n] += t_real2[n] * w20 + t_real1[n] * w21 + t_real2[n] * w22; 
				r_imag2[n] += t_imag2[n] * w20 + t_imag1[n] * w21 + t_imag2[n] * w22; 
				    rw2[n] += ctf    [n] * w20 + ctf    [n] * w21 + ctf    [n] * w22;
			}

			return;
		}

		assert(rotationLen == 1);
		static_assert(maxTranslationIndexsCapacity <= 5, "Need more or fewer cases");
		switch (translationIndexsLen) {
		case 5: {
#ifndef MY_MACHINE_HAS_AVX512
			std::cerr << "BackProjection_Kernel_implementation using case 5 but not enough vector registers" << std::endl;
			std::cout << "BackProjection_Kernel_implementation using case 5 but not enough vector registers" << std::endl;
			EXIT_ABNORMALLY;
#endif
			HEADSTEP(0)
			HEADSTEP(1)
			HEADSTEP(2)
			HEADSTEP(3)
			HEADSTEP(4)
			#pragma vector aligned
    	    #pragma ivdep
			for (int n = sLo; n < sHi; n++) {
				if (false) {
					BODYSTEP(0)
					BODYSTEP(1)
					BODYSTEP(2)
					BODYSTEP(3)
					BODYSTEP(4)
				} else {
					r_real0[n] += t_real0[n] * w00 
					           +  t_real1[n] * w01 
					           +  t_real2[n] * w02 
					           +  t_real3[n] * w03
							   +  t_real4[n] * w04; 
					r_imag0[n] += t_imag0[n] * w00 
					           +  t_imag1[n] * w01 
					           +  t_imag2[n] * w02 
					           +  t_imag3[n] * w03 
					           +  t_imag4[n] * w04; 
					rw0[n]     += ctf    [n] * w00	
					           +  ctf    [n] * w01	
					           +  ctf    [n] * w02	
					           +  ctf    [n] * w03	
					           +  ctf    [n] * w04;	
				}
			}
		}; break;
		case 4: {
			HEADSTEP(0)
			HEADSTEP(1)
			HEADSTEP(2)
			HEADSTEP(3)
			#pragma vector aligned
    	    #pragma ivdep
			for (int n = sLo; n < sHi; n++) {
				if (false) {
					BODYSTEP(0)
					BODYSTEP(1)
					BODYSTEP(2)
					BODYSTEP(3)
				} else {
					// The compiler was generating unneccessary fetches and stores with the macros version, 
					// so manually expand and reorder it to get marginally better code
					//
					// BODYSTEP_R(0,0)
					// BODYSTEP_R(0,1)
					// BODYSTEP_R(0,2)
					// BODYSTEP_R(0,3)
					// {   r_real0[n] += t_real0[n] * w00; r_imag0[n] += t_imag0[n] * w00; rw0[n] += ctf[n] * w00; }
					// {   r_real0[n] += t_real1[n] * w01; r_imag0[n] += t_imag1[n] * w01; rw0[n] += ctf[n] * w01; }
					// {   r_real0[n] += t_real2[n] * w02; r_imag0[n] += t_imag2[n] * w02; rw0[n] += ctf[n] * w02; }
					// {   r_real0[n] += t_real3[n] * w03; r_imag0[n] += t_imag3[n] * w03; rw0[n] += ctf[n] * w03; }
					r_real0[n] += t_real0[n] * w00 + t_real1[n] * w01 + t_real2[n] * w02 + t_real3[n] * w03; 
					r_imag0[n] += t_imag0[n] * w00 + t_imag1[n] * w01 + t_imag2[n] * w02 + t_imag3[n] * w03; 
					rw0[n]     += ctf    [n] * w00 + ctf    [n] * w01 + ctf    [n] * w02 + ctf    [n] * w03;
				}
			}
		}; break;
		case 3: {
			HEADSTEP(0)
			HEADSTEP(1)
			HEADSTEP(2)
			#pragma vector aligned
    	    #pragma ivdep
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP(0)
				BODYSTEP(1)
				BODYSTEP(2)
			}
		}; break;
		case 2: {
			HEADSTEP(0)
			HEADSTEP(1)
			#pragma vector aligned
    	    #pragma ivdep
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP(0)
				BODYSTEP(1)
			}
		}; break;
		case 1: {
			HEADSTEP(0)
			#pragma vector aligned
    	    #pragma ivdep
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP(0)
			}
		}; break;
		default: assert(false);
		}	// switch

#undef BODYSTEP
#undef BODYSTEP_R
#undef HEADSTEP
#undef HEADSTEP_R
#undef HEADSTEP_SHARED

		// This can do the tail, if sLo updated :-) which it isn't
		//
		if (false) {
		    for (int r = 0; r < step.rotationLen; r++) {
		    	auto rotation = step.rotationBegin + r;
		    	for (int ti = 0; ti < step.translationIndexsLen; ti++) {
		    		auto translation = step.translationIndexs[ti];
		    
		    	    auto rotation_real		= rotations_real   [rotation   ].ptrValidUpto(sHi); 
		    	    auto rotation_imag		= rotations_imag   [rotation   ].ptrValidUpto(sHi); 
		    	    auto rotation_weight	= rotations_weight [rotation   ].ptrValidUpto(sHi);
		    	    auto translation_real	= translations_real[translation].ptrValidUpto(sHi); 
		    	    auto translation_imag	= translations_imag[translation].ptrValidUpto(sHi); 
		    	    auto weight             = weights[weights.index(rotation, translation)];
		    
		    		if (0) {
		    			if (rotation == 0 && translation == 0) {
		    				std::cout << "BackProjection_Kernel_implementation::execute: w:" << weight << std::endl;
		    				for (int n = sLo; n < sHi; n++) {
		    					std::cout << "n:" << n 
		    						<< " tr:" << translation_real[n]
		    						<< " ti:" << translation_imag[n]
		    						<< " ct:" << ctf[n]
		    						<< std::endl;
		    				}
		    			}
		    		}
		    
		    	    #pragma ivdep
		    	    for (int n = sLo; n < sHi; n++)
		    	    {
		    	    	rotation_real  [n] += translation_real[n] * weight;
		    	    	rotation_imag  [n] += translation_imag[n] * weight;
		    	    	// now Fweight stores sum of all w
		    	    	// Note that CTF needs to be squared in Fweight, weightxinvsigma2 already contained one copy
		    	    	rotation_weight[n] += ctf[n] * weight;
		    	    }
		    	}
			}
		}
	}

	virtual void compute() {
		// numbers found by experimention
		compute(48, 64, 256);
	}

	virtual void compute(const S64 rotationStepProposed, const S64 translationStepProposed, const S64 sampleStep) {
		Scheduler::compute(rotationStepProposed, translationStepProposed, sampleStep);
	}

	virtual void releaseCapturedOutputs() {
		TUNING_SCOPE_STEP(bpk_releaseCapturedOutputs)
		// This is safe because the original is thread private
		for (int rotation = 0; rotation < numberOfRotations; rotation++) {
			auto& rr = rotations_real[rotation];
			if (rr.isCapturedForSummingInto()) {
				rotationState->setRotationState(iclass, rotation, rotationState->RS_addPending, 0, rr.ptrToSumInto());
				rotations_real  [rotation].sumIntoOriginal();	// Note: the original is the accumulator, not the exp_Fref...
				rotations_imag  [rotation].sumIntoOriginal();	// Note: the original is the accumulator, not the exp_Fref...
				rotations_weight[rotation].sumIntoOriginal();	// Note: the original is the accumulator, not the exp_Fref...
			}
		}
	}
};

class BpkPool {
	struct Cell {
		BackProjection_Kernel_implementation* v;
		Cell() : v(nullptr) {}
		~Cell() { deallocate(); }
		bool matches(BackProjection_Kernel_implementation* v) { return this->v == v; }
		void allocate(BackProjection_Kernel_implementation* v) {
			assert(!this->v);
			this->v = v;
		}
		void deallocate() {
			sDelete(v);
		}
	};
	ListPerThreadTemplate<Cell> cellsPerThread;
public:
	BpkPool() {}
	BackProjection_Kernel_implementation* allocate(
		int iter,
		RotationState& rotationState,
		int iclass,
		int numberOfTranslations,
		int numberOfSamples) {
		auto & cell = *cellsPerThread.acquire();
		if (!cell.v || !cell.v->suitable(rotationState.numberOfRotations, numberOfTranslations, numberOfSamples)) { 
			cell.deallocate();
			cell.allocate(BackProjection_Kernel_implementation::make(maxTranslationIndexsCapacity, rotationState.numberOfRotations, numberOfTranslations, numberOfSamples));
		}
		cell.v->init(iter, rotationState, iclass);
		return cell.v;
	}
	void deallocate(BackProjection_Kernel_implementation* v) {
		cellsPerThread.release(v);
	}
} bpk_pool;


BackProjection_Kernel* BackProjection_Kernel::acquire(
	int iter,
	RotationState& rotationState,
	int iclass,
	int numberOfTranslations,
	int numberOfSamples) 
{
	return bpk_pool.allocate(iter, rotationState, iclass, numberOfTranslations, numberOfSamples);
}

void BackProjection_Kernel_implementation::release() {
	Scheduler::fini();
	this->rotationState = NULL;
	this->iclass        = -1;
	bpk_pool.deallocate(this);
}

class GetAllSquaredDifferences_Kernel_implementation : public GetAllSquaredDifferences_Kernel, Scheduler {
public:
	typedef GetAllSquaredDifferences_Kernel_implementation Kernel;

	static GetAllSquaredDifferences_Kernel_implementation* make(
		const int translationIndexsCapacity,
		const int numberOfRotations,
		const int numberOfTranslations,
		const int numberOfSamples) 
	{
#include "./util_heap_undefs.h"
		return sNewA(GetAllSquaredDifferences_Kernel_implementation, (
			translationIndexsCapacity,
			numberOfRotations,
			numberOfTranslations,
			numberOfSamples));
#include "./util_heap_defs.h"
	}

	GetAllSquaredDifferences_Kernel_implementation(
		const int translationIndexsCapacity,
		const int numberOfRotations,
		const int numberOfTranslations,
		const int numberOfSamples) 
	  : Scheduler(SchedulerKind_GetAllSquaredDifferences, translationIndexsCapacity, numberOfRotations, numberOfTranslations, numberOfSamples),
		rotations_real	 (numberOfRotations),
		rotations_imag	 (numberOfRotations),
		translations_real(numberOfTranslations),
		translations_imag(numberOfTranslations),
		dist			 (*this)
	{
	}

	virtual void release();

	bool suitable(
		const int numberOfRotations,
		const int numberOfTranslations,
		const int numberOfSamples) {

		return true
			&& this->numberOfRotations	  == numberOfRotations
			&& this->numberOfTranslations == numberOfTranslations
			&& this->numberOfSamples	  == numberOfSamples;
	}

	void init(int iter, const double* weightsDouble)
	{
		Scheduler::init(iter);
		for (int i = 0; i < numberOfRotations; i++) {
			rotations_real[i].init();
			rotations_imag[i].init();
		}
		for (int i = 0; i < numberOfTranslations; i++) {
			translations_real[i].init();
			translations_imag[i].init();
		}
		dist.init();
		weightsSqrt.init();
		appended.resize(0);
		weightsSqrt.captureSqrt(&capturedDataStats, "weightsSqrt", numberOfSamples, weightsDouble);
	}

	~GetAllSquaredDifferences_Kernel_implementation() {
	}

	CapturedData weightsSqrt;
	std::vector<CapturedData> translations_real, translations_imag, rotations_real, rotations_imag;

	struct Appended {
		int rotation;
		int translation;
		int rot_trans_over;
		Appended() {}
		Appended(
			int rotation,
			int translation,
			int rot_trans_over) : rotation(rotation), translation(translation), rot_trans_over(rot_trans_over) {}
	};
	std::vector<Appended> appended;

	virtual void append(
		int rotation,
		int translation,
		int rot_trans_over,
		const double* rotation_real,    const double* rotation_imag,
		const double* translation_real, const double* translation_imag,
		bool interesting) {
		Scheduler::append(rotation,translation);
		appended.push_back(Appended(rotation,translation,rot_trans_over));
		if (interesting) {
			rotations_real   [   rotation].setInteresting();
			rotations_imag   [   rotation].setInteresting();
			translations_real[translation].setInteresting();
			translations_imag[translation].setInteresting();
		}
		rotations_real   [   rotation].captureMultipliedByWeightsSqrtForIS(&capturedDataStats, "gasd_rotations_real",    numberOfSamples, rotation_real   , weightsSqrt);
		rotations_imag   [   rotation].captureMultipliedByWeightsSqrtForIS(&capturedDataStats, "gasd_rotations_imag",    numberOfSamples, rotation_imag   , weightsSqrt);
		translations_real[translation].captureMultipliedByWeightsSqrtForIS(&capturedDataStats, "gasd_translations_real", numberOfSamples, translation_real, weightsSqrt);
		translations_imag[translation].captureMultipliedByWeightsSqrtForIS(&capturedDataStats, "gasd_translations_imag", numberOfSamples, translation_imag, weightsSqrt);
	}


	class Dist {
		Kernel & parent;
		const S64 vLen;
		double*   v;
		bool	  zeroed;
	public:
		Dist(Kernel & parent) 
		  : parent(parent),
			vLen(parent.numberOfRotations*parent.numberOfTranslations),
			v(mallocAlignedDoubleNotZeroed(vLen))
		{
			init();
		}
		void init() {
			zeroed = false;
		}
		~Dist() {
			aFree(v);
		}
		void zero() {
			for (S64 i = 0; i < vLen; i++) v[i] = 0;
			zeroed = true;
		}
		void zero(S64 r, S64 t) {
			v[index(r,t)] = 0.0;
		}
		const double& operator()(S64 r, S64 t) { return v[index(r,t)]; }

		void add(S64 r, S64 t, double value) {
			checkNotNegOrNan(value);
			v[index(r,t)] += value;
		}

		S64 index(S64 r, S64 t) {
			assert(0 <= r); assert(r < parent.numberOfRotations);
			assert(0 <= t); assert(t < parent.numberOfTranslations);
			return r*parent.numberOfTranslations + t;
		};
	} dist;

	virtual void executeSingletons(SingletonStep const * steps, int numberOfSteps) {
		// std::cerr << "executeSingletons numberOfSteps:"<< numberOfSteps << std::endl;

		for (int step = 0; step < numberOfSteps;) {
			static const int MaxNumberOfStepsExecuted = 3;
			int			     numberOfStepsExecuted    = std::min(MaxNumberOfStepsExecuted, numberOfSteps - step);

			const bool logging = false;

			if (logging) steps[step].print("Executing");

#define HEADSTEP(SD) \
			const auto rotation##SD    = steps[step+SD].rotation;										\
			const auto translation##SD = steps[step+SD].translation;									\
			auto r_real##SD = rotations_real[rotation##SD].ptrToOriginalConst();						\
			auto r_imag##SD = rotations_imag[rotation##SD].ptrToOriginalConst();						\
			assert(!!r_real##SD);																		\
			assert(!!r_imag##SD);																		\
			auto t_real##SD = translations_real[translation##SD].ptrToOriginalConst();					\
			auto t_imag##SD = translations_imag[translation##SD].ptrToOriginalConst();					\
			assert(!!t_real##SD);																		\
			assert(!!t_imag##SD);																		\
			double dist##SD = 0;																		// end of macro
			
			const double* w = this->weightsSqrt.ptrToOriginalConst();

		#define BODYSTEP(SD) \
			{	double diff_real = r_real##SD[n] - t_real##SD[n];										\
				double diff_imag = r_imag##SD[n] - t_imag##SD[n];										\
				double dist_addend = (diff_real * diff_real + diff_imag * diff_imag) * w[n];			\
				checkNotNegOrNan(dist_addend);															\
				dist##SD += dist_addend;																\
			}																							// end of macro

		#define TAILSTEP(SD) \
			dist.add(rotation##SD, translation##SD, dist##SD);											// end of macro

			switch (numberOfStepsExecuted) {
			case 3: {
				HEADSTEP(0)
				HEADSTEP(1)
				HEADSTEP(2)
				for (int n = 0; n < numberOfSamples; n++) {
					BODYSTEP(0)
					BODYSTEP(1)
					BODYSTEP(2)
				}
				TAILSTEP(0)
				TAILSTEP(1)
				TAILSTEP(2)
			}; break;
			case 2: {
				HEADSTEP(0)
				HEADSTEP(1)
				for (int n = 0; n < numberOfSamples; n++) {
					BODYSTEP(0)
					BODYSTEP(1)
				}
				TAILSTEP(0)
				TAILSTEP(1)
			}; break;
			case 1: {
				HEADSTEP(0)
				for (int n = 0; n < numberOfSamples; n++) {
					BODYSTEP(0)
				}
				TAILSTEP(0)
			}; break;
			default:
				std::cerr << "GetAllSquaredDifferences_Kernel::executeSingletons: numberOfStepsExecuted too large" << std::endl;
				EXIT_ABNORMALLY;
			}

#undef TAILSTEP
#undef BODYSTEP
#undef HEADSTEP

			step += numberOfStepsExecuted;
		}
	}


	virtual void executeMergedStep(MergedStep const & step, Plan& plan, int sLo, int sHi) {
		auto  translationIndexsLen = step.translationIndexsLen;
		auto  rotationBegin        = step.rotationBegin;
		auto  rotationLen		   = step.rotationLen;
		auto& translationIndexs    = step.translationIndexs;
		
		if (translationIndexsLen == 0) return;

		const bool logging = false && rotationBegin==0 && translationIndexs[0]==0;

		if (logging) step.print("Executing");

		float* r_real0 = rotations_real[rotationBegin + 0].ptrValidUpto(sHi);
		float* r_imag0 = rotations_imag[rotationBegin + 0].ptrValidUpto(sHi); 
		assert(!!r_real0);
		assert(!!r_imag0);
		assert(U64(r_real0) % 64 == 0); 
		assert(U64(r_imag0) % 64 == 0); 

		#define HEADSTEP_SHARED(TD) \
			float* t_real##TD = translations_real[translationIndexs[TD]].ptrValidUpto(sHi);	\
			float* t_imag##TD = translations_imag[translationIndexs[TD]].ptrValidUpto(sHi);	\
			assert(!!t_real##TD);															\
			assert(!!t_imag##TD);															\
			assert(U64(t_real##TD) % 64 == 0);												\
			assert(U64(t_imag##TD) % 64 == 0);												\
			// end of macro

		#define HEADSTEP_R(RD, TD) \
			double dist##RD##TD = 0;

		#define HEADSTEP(TD) \
			HEADSTEP_SHARED(TD) HEADSTEP_R(0, TD)
				
		#define BODYSTEP_R(RD, TD) \
			{	float diff_real = r_real##RD[n] - t_real##TD[n];							\
				float diff_imag = r_imag##RD[n] - t_imag##TD[n];							\
				float dist_addend = (diff_real * diff_real + diff_imag * diff_imag);		\
				checkNotNegOrNan(dist_addend);												\
				dist##RD##TD += dist_addend;												\
			}																				// end of macro

		#define BODYSTEP(TD) \
			BODYSTEP_R(0,TD)

		#define TAILSTEP_R(RD,TD) \
			dist.add(rotationBegin + RD, translationIndexs[TD], dist##RD##TD);

		#define TAILSTEP(TD) \
			TAILSTEP_R(0,TD)

		if (rotationLen == 3) {
			// 3x3 is done as simply three repeats of a 1x3
			// It is the only case known how to do
			//
			assert(rotationLen  == 3);
			assert(translationIndexsLen == 3);

			float* r_real1 = rotations_real[rotationBegin + 1].ptrValidUpto(sHi);
			float* r_imag1 = rotations_imag[rotationBegin + 1].ptrValidUpto(sHi); 
			float* r_real2 = rotations_real[rotationBegin + 2].ptrValidUpto(sHi);
			float* r_imag2 = rotations_imag[rotationBegin + 2].ptrValidUpto(sHi); 
			assert(U64(r_real1) % 64 == 0); 
			assert(U64(r_imag1) % 64 == 0); 
			assert(U64(r_real2) % 64 == 0); 
			assert(U64(r_imag2) % 64 == 0); 

			HEADSTEP_SHARED(0) HEADSTEP_R(0,0)  HEADSTEP_R(1,0)  HEADSTEP_R(2,0)
			HEADSTEP_SHARED(1) HEADSTEP_R(0,1)  HEADSTEP_R(1,1)  HEADSTEP_R(2,1)
			HEADSTEP_SHARED(2) HEADSTEP_R(0,2)  HEADSTEP_R(1,2)  HEADSTEP_R(2,2)

			#pragma vector aligned
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP_R(0,0)  BODYSTEP_R(1,0)  BODYSTEP_R(2,0)
				BODYSTEP_R(0,1)  BODYSTEP_R(1,1)  BODYSTEP_R(2,1)
				BODYSTEP_R(0,2)  BODYSTEP_R(1,2)  BODYSTEP_R(2,2)
			}

			TAILSTEP_R(0,0)	TAILSTEP_R(1,0)	TAILSTEP_R(2,0)
			TAILSTEP_R(0,1)	TAILSTEP_R(1,1)	TAILSTEP_R(2,1)
			TAILSTEP_R(0,2)	TAILSTEP_R(1,2)	TAILSTEP_R(2,2)

			return;
		}

		assert(rotationLen == 1);
		static_assert(maxTranslationIndexsCapacity <= 5, "Need more or fewer cases");
		switch (translationIndexsLen) {
		case 5: {
			HEADSTEP(0)
			HEADSTEP(1)
			HEADSTEP(2)
			HEADSTEP(3)
			HEADSTEP(4)

			#pragma vector aligned
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP(0)
				BODYSTEP(1)
				BODYSTEP(2)
				BODYSTEP(3)
				BODYSTEP(4)
			}

			TAILSTEP(0)
			TAILSTEP(1)
			TAILSTEP(2)
			TAILSTEP(3)
			TAILSTEP(4)
		}; break;
		case 4: {
			HEADSTEP(0)
			HEADSTEP(1)
			HEADSTEP(2)
			HEADSTEP(3)
			#pragma vector aligned
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP(0)
				BODYSTEP(1)
				BODYSTEP(2)
				BODYSTEP(3)
			}
			TAILSTEP(0)
			TAILSTEP(1)
			TAILSTEP(2)
			TAILSTEP(3)
		}; break;
		case 3: {
			HEADSTEP(0)
			HEADSTEP(1)
			HEADSTEP(2)
			#pragma vector aligned
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP(0)
				BODYSTEP(1)
				BODYSTEP(2)
			}
			TAILSTEP(0)
			TAILSTEP(1)
			TAILSTEP(2)
		}; break;
		case 2: {
			HEADSTEP(0)
			HEADSTEP(1)
			#pragma vector aligned
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP(0)
				BODYSTEP(1)
			}
			TAILSTEP(0)
			TAILSTEP(1)
		}; break;
		case 1: {
			HEADSTEP(0)
			#pragma vector aligned
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP(0)
			}
			TAILSTEP(0)
		}; break;
		default: assert(false);
		}	// switch

#undef TAILSTEP
#undef TAILSTEP_R
#undef BODYSTEP
#undef BODYSTEP_R
#undef HEADSTEP
#undef HEADSTEP_R
#undef HEADSTEP_SHARED

		if (logging) {
			double partial_diff200;
			getRT(partial_diff200, 0,0);
			std::cout << "Diff2:" << partial_diff200 << std::endl;
			std::cout << std::endl;
		}
	}

	virtual void compute() {
		// numbers found by experimention
		compute(64, 64, 512);
	}

	virtual void compute(const S64 rotationStepProposed, const S64 translationStepProposed, const S64 sampleStep) {
		dist.zero();
		Scheduler::compute(rotationStepProposed, translationStepProposed, sampleStep);
	}

	virtual void recompute(int rotation, int translation) {
		const S64 sampleStep = 256;
		for (int i = 0; i < 2; i++) {
			Plan plan(*this, rotation, rotation+1, translation, translation+1, sampleStep, false);
			auto const & d = dist(rotation, translation);
			dist.zero(rotation, translation);
			for (int sLo = 0; sLo < numberOfSamples; sLo += sampleStep) {
				auto sHi = std::min(numberOfSamples, sLo +  sampleStep);
				plan.executeMergedSteps(sLo, sHi);
			}
			std::cout << "Recomputed as kernel_diff2:" << d << std::endl;
		}
	}

	virtual bool getRT(
		double&	diff2,
		int		rotation,
		int		translation) {
		if (!mustDoIsSet(rotation, translation)) { diff2 = -666.777; return false; }
		diff2 = dist(rotation, translation);
		return true;
	}
	virtual bool getI(
		double&	diff2,
		int&    rot_trans_over,
		int		i) {
		if (i >= appended.size()) { diff2 = -666.777; rot_trans_over = 0; return false; }
		auto& a = appended[i];
		diff2 = dist(a.rotation, a.translation);
		rot_trans_over = a.rot_trans_over;
		return true;
	}
};


class GasdkPool {
	struct Cell {
		GetAllSquaredDifferences_Kernel_implementation* v;
		Cell() : v(nullptr) {}
		~Cell() { deallocate(); }
		bool matches(GetAllSquaredDifferences_Kernel_implementation * v) { return this->v == v; }
		void allocate(GetAllSquaredDifferences_Kernel_implementation* v) {
			assert(!this->v);
			this->v = v;
		}
		void deallocate() {
			sDelete(v);
		}
	};
	ListPerThreadTemplate<Cell> cellsPerThread;
public:
	GasdkPool() {}
	GetAllSquaredDifferences_Kernel_implementation* allocate(
		int iter,
		int numberOfRotations,
		int numberOfTranslations,
		int numberOfSamples,
		const double* weights) {
		auto & cell = *cellsPerThread.acquire();
		if (!cell.v || !cell.v->suitable(numberOfRotations, numberOfTranslations, numberOfSamples)) { 
			cell.deallocate();
			auto v = GetAllSquaredDifferences_Kernel_implementation::make(maxTranslationIndexsCapacity, numberOfRotations, numberOfTranslations, numberOfSamples);
			cell.allocate(v);
		}
		cell.v->init(iter, weights);
		return cell.v;
	}
	void deallocate(GetAllSquaredDifferences_Kernel_implementation* v) {
		cellsPerThread.release(v);
	}
} gasdk_pool;


GetAllSquaredDifferences_Kernel* GetAllSquaredDifferences_Kernel::acquire(
	int iter,
	int numberOfRotations,
	int numberOfTranslations,
	int numberOfSamples,
	const double* weights) 
{
	return gasdk_pool.allocate(iter, numberOfRotations, numberOfTranslations, numberOfSamples, weights);
}

void GetAllSquaredDifferences_Kernel_implementation::release() {
	Scheduler::fini();
	gasdk_pool.deallocate(this);
}


class UpdateOtherParams_Kernel_implementation : public UpdateOtherParams_Kernel, Scheduler {
public:
	typedef UpdateOtherParams_Kernel_implementation Kernel;

	static UpdateOtherParams_Kernel_implementation* make(
		const int translationIndexsCapacity,
		const int numberOfRotations,
		const int numberOfTranslations,
		const int numberOfSamples)
	{
#include "./util_heap_undefs.h"
		return sNewA(UpdateOtherParams_Kernel_implementation, (
			translationIndexsCapacity,
			numberOfRotations,
			numberOfTranslations,
			numberOfSamples));
#include "./util_heap_defs.h"
	}

	UpdateOtherParams_Kernel_implementation(
		const int translationIndexsCapacity,
		const int numberOfRotations,
		const int numberOfTranslations,
		const int numberOfSamples)
	  : Scheduler(SchedulerKind_UpdateOtherParams,
#ifdef MY_MACHINE_HAS_AVX512
			translationIndexsCapacity,
#else
			std::min(4, translationIndexsCapacity),
				// the 5x1 uses too many vector registers
#endif
			numberOfRotations, numberOfTranslations, numberOfSamples),
	    sparseSampleIndexs  (numberOfSamples),
		rotations_real	    (numberOfRotations),
		rotations_imag	    (numberOfRotations),
		translations_real   (numberOfTranslations),
		translations_imag   (numberOfTranslations),
		weights 			(SchedulerKind_UpdateOtherParams, numberOfRotations,numberOfTranslations),
		dist				(*this),
		sparse_wsum_sigma2_noise(NULL),
		dense_wsum_sigma2_noise_available(0),
		dense_wsum_sigma2_noise_double(NULL),
		dense_wsum_sigma2_noise_float (NULL)
	{
	}

	bool suitable(
		const int numberOfRotations,
		const int numberOfTranslations,
		const int numberOfSamples) {

		return true
			&& this->numberOfRotations	  == numberOfRotations
			&& this->numberOfTranslations == numberOfTranslations
			&& this->numberOfSamples	  == numberOfSamples;
	}

	void init(
		int iter,
		const int* Mresol_fine,
		double* wsum_sigma2_noise) {

		Scheduler::init(iter);

		for (int i = 0; i < numberOfRotations; i++) {
			rotations_real[i].init();
			rotations_imag[i].init();
		}
		for (int i = 0; i < numberOfTranslations; i++) {
			translations_real[i].init();
			translations_imag[i].init();
		}
		weights.init();
		appended.resize(0);
		dist.init();

		sparse_wsum_sigma2_noise = wsum_sigma2_noise;
	
		sparseSampleIndexs.init_Mresol_fine(Mresol_fine);

		auto needed = sparseSampleIndexs.numberOfAcceptedSamples();

		if (needed > dense_wsum_sigma2_noise_available) {
			needed += needed/3;		// grow exponentially
			needed = std::max(needed, 128);
			dense_wsum_sigma2_noise_double = mallocAlignedDoubleNotZeroed(needed);
			dense_wsum_sigma2_noise_float  = mallocAlignedFloatNotZeroed (needed);
			dense_wsum_sigma2_noise_available = needed;
		}

	}

	virtual void release();

	~UpdateOtherParams_Kernel_implementation() {
		aFree(dense_wsum_sigma2_noise_double);
		aFree(dense_wsum_sigma2_noise_float);
	}

	SparseSampleIndexs sparseSampleIndexs;
	std::vector<CapturedData> translations_real, translations_imag, rotations_real, rotations_imag;
	FloatsForRT weights;
	double* sparse_wsum_sigma2_noise;

	int     dense_wsum_sigma2_noise_available;
	float*  dense_wsum_sigma2_noise_float;
	double* dense_wsum_sigma2_noise_double;

	virtual int numberOfAcceptedSamples() const { 
		return sparseSampleIndexs.numberOfAcceptedSamples(); 
	}

	struct Appended {
		int rotation;
		int translation;
		Appended() {}
		Appended(
			int rotation,
			int translation) : rotation(rotation), translation(translation) {}
	};
	std::vector<Appended> appended;

	virtual void append(
		int rotation,
		int translation,
		const double* rotation_real,    const double* rotation_imag,
		const double* translation_real, const double* translation_imag, 
		double weight,
		bool interesting) {
		Scheduler::append(rotation,translation);
		appended.push_back(Appended(rotation,translation));
		weights[weights.index(rotation, translation)] = checkNotNanOrInf(weight);
		if (interesting) {
			rotations_real   [   rotation].setInteresting();
			rotations_imag   [   rotation].setInteresting();
			translations_real[translation].setInteresting();
			translations_imag[translation].setInteresting();
		}
		rotations_real   [   rotation].captureSparse(&capturedDataStats, "uop_rotations_real",    numberOfSamples, rotation_real   , sparseSampleIndexs);
		rotations_imag   [   rotation].captureSparse(&capturedDataStats, "uop_rotations_imag",    numberOfSamples, rotation_imag   , sparseSampleIndexs);
		translations_real[translation].captureSparse(&capturedDataStats, "uop_translations_real", numberOfSamples, translation_real, sparseSampleIndexs);
		translations_imag[translation].captureSparse(&capturedDataStats, "uop_translations_imag", numberOfSamples, translation_imag, sparseSampleIndexs);
	}

	class Dist {
		Kernel & parent;
		const S64 vLen;
		double*   v;
		bool	  zeroed;
	public:
		Dist(Kernel & parent) 
		  : parent(parent),
			vLen(parent.numberOfRotations*parent.numberOfTranslations),
			v(mallocAlignedDoubleNotZeroed(vLen))
		{
			init();
		}
		void init() {
			zeroed = false;
		}
		~Dist() {
			aFree(v);
		}
		void zero() {
			for (S64 i = 0; i < vLen; i++) v[i] = 0;
			zeroed = true;
		}
		void zero(S64 r, S64 t) {
			v[index(r,t)] = 0.0;
		}
		const double& operator()(S64 r, S64 t) { return v[index(r,t)]; }

		void add(S64 r, S64 t, double value) {
			checkNotNegOrNan(value);
			v[index(r,t)] += value;
		}

		S64 index(S64 r, S64 t) {
			assert(0 <= r); assert(r < parent.numberOfRotations);
			assert(0 <= t); assert(t < parent.numberOfTranslations);
			return r*parent.numberOfTranslations + t;
		};
	} dist;

	virtual void executeSingletons(SingletonStep const * steps, int numberOfSteps) {
		
		typedef void numberOfSamples;
			// Poison this because the of sparseness

		// std::cerr << "executeSingletons numberOfSteps:"<< numberOfSteps << std::endl;

		const int numberOfAcceptedSamples = sparseSampleIndexs.numberOfAcceptedSamples();

		for (int step = 0; step < numberOfSteps;) {
			static const int MaxNumberOfStepsExecuted = 3;
			int			     numberOfStepsExecuted    = std::min(MaxNumberOfStepsExecuted, numberOfSteps - step);

			const bool logging = false;

			if (logging) steps[step].print("Executing");

		#define HEADSTEP(SD) \
			const auto rotation##SD    = steps[step+SD].rotation;										\
			const auto translation##SD = steps[step+SD].translation;									\
			auto r_real##SD = rotations_real[rotation##SD].ptrToOriginalConst();						\
			auto r_imag##SD = rotations_imag[rotation##SD].ptrToOriginalConst();						\
			assert(!!r_real##SD);																		\
			assert(!!r_imag##SD);																		\
			auto t_real##SD = translations_real[translation##SD].ptrToOriginalConst();					\
			auto t_imag##SD = translations_imag[translation##SD].ptrToOriginalConst();					\
			assert(!!t_real##SD);																		\
			assert(!!t_imag##SD);																		\
			double weight##SD = weights[weights.index(rotation##SD,translation##SD)];                   \
			double dist##SD = 0;																		// end of macro
		
		#define BODYSTEP_BEFORE double add_wsum_sigma2_noise = 0.0;

		#define BODYSTEP(SD) \
			{	double diff_real = r_real##SD[n] - t_real##SD[n];										\
				double diff_imag = r_imag##SD[n] - t_imag##SD[n];										\
				double dist_addend = (diff_real * diff_real + diff_imag * diff_imag) * weight##SD;		\
				checkNotNegOrNan(dist_addend);															\
				dist##SD += dist_addend;																\
				add_wsum_sigma2_noise += dist_addend;													\
			}																							// end of macro

		#define BODYSTEP_AFTER  dense_wsum_sigma2_noise_double[n] += add_wsum_sigma2_noise;				// end of macro

		#define TAILSTEP(SD) \
			dist.add(rotation##SD, translation##SD, dist##SD);											// end of macro

			switch (numberOfStepsExecuted) {
			case 3: {
				HEADSTEP(0)
				HEADSTEP(1)
				HEADSTEP(2)
    		    #pragma ivdep
				for (int n = 0; n < numberOfAcceptedSamples; n++) {
					BODYSTEP_BEFORE
					BODYSTEP(0)
					BODYSTEP(1)
					BODYSTEP(2)
					BODYSTEP_AFTER
				}
				TAILSTEP(0)
				TAILSTEP(1)
				TAILSTEP(2)
			}; break;
			case 2: {
				HEADSTEP(0)
				HEADSTEP(1)
    		    #pragma ivdep
				for (int n = 0; n < numberOfAcceptedSamples; n++) {
					BODYSTEP_BEFORE
					BODYSTEP(0)
					BODYSTEP(1)
					BODYSTEP_AFTER
				}
				TAILSTEP(0)
				TAILSTEP(1)
			}; break;
			case 1: {
				HEADSTEP(0)
    		    #pragma ivdep
				for (int n = 0; n < numberOfAcceptedSamples; n++) {
					BODYSTEP_BEFORE
					BODYSTEP(0)
					BODYSTEP_AFTER
				}
				TAILSTEP(0)
			}; break;
			default:
				std::cerr << "UpdateOtherParams_Kernel::executeSingletons: numberOfStepsExecuted too large" << std::endl;
				EXIT_ABNORMALLY;
			}

#undef TAILSTEP
#undef BODYSTEP_AFTER
#undef BODYSTEP
#undef BODYSTEP_BEFORE
#undef HEADSTEP

			step += numberOfStepsExecuted;
		}
	}


	virtual void executeMergedStep(MergedStep const & step, Plan& plan, int sLo, int sHi) {
		auto  translationIndexsLen = step.translationIndexsLen;
		auto  rotationBegin        = step.rotationBegin;
		auto  rotationLen		   = step.rotationLen;
		auto& translationIndexs    = step.translationIndexs;
		
		if (translationIndexsLen == 0) return;

		const bool logging = false && rotationBegin==0 && translationIndexs[0]==0;

		if (logging) step.print("Executing");

		float* r_real0 = rotations_real[rotationBegin + 0].ptrValidUpto(sHi);
		float* r_imag0 = rotations_imag[rotationBegin + 0].ptrValidUpto(sHi); 
		assert(!!r_real0);
		assert(!!r_imag0);
		assert(U64(r_real0) % 64 == 0); 
		assert(U64(r_imag0) % 64 == 0); 

		#define HEADSTEP_SHARED(TD) \
			float* t_real##TD = translations_real[translationIndexs[TD]].ptrValidUpto(sHi);	\
			float* t_imag##TD = translations_imag[translationIndexs[TD]].ptrValidUpto(sHi);	\
			assert(!!t_real##TD);															\
			assert(!!t_imag##TD);															\
			assert(U64(t_real##TD) % 64 == 0);												\
			assert(U64(t_imag##TD) % 64 == 0);												\
			// end of macro

		#define HEADSTEP_R(RD, TD) \
			float dist##RD##TD = 0;

		#define HEADSTEP(TD) \
			HEADSTEP_SHARED(TD) HEADSTEP_R(0, TD)
				
		#define BODYSTEP_BEFORE \
			float add_wsum_sigma2_noise = 0;														// end of macro

		#define BODYSTEP_R(RD, TD) \
			{	float diff_real = r_real##RD[n] - t_real##TD[n];									\
				float diff_imag = r_imag##RD[n] - t_imag##TD[n];									\
				float weight##SD = weights[weights.index(rotationBegin+RD,translationIndexs[TD])];	\
				float dist_addend = (diff_real * diff_real + diff_imag * diff_imag) * weight##SD;	\
				checkNotNegOrNan(dist_addend);														\
				dist##RD##TD += dist_addend;														\
				add_wsum_sigma2_noise += dist_addend;												\
			}																						// end of macro

		#define BODYSTEP_AFTER  dense_wsum_sigma2_noise_float[n] += add_wsum_sigma2_noise;			// end of macro

		#define BODYSTEP(TD) \
			BODYSTEP_R(0,TD)

		#define TAILSTEP_R(RD,TD) \
			dist.add(rotationBegin + RD, translationIndexs[TD], dist##RD##TD);

		#define TAILSTEP(TD) \
			TAILSTEP_R(0,TD)

		if (rotationLen == 3) {
			// 3x3 is done as simply three repeats of a 1x3
			// It is the only case known how to do
			//
			assert(rotationLen  == 3);
			assert(translationIndexsLen == 3);

			float* r_real1 = rotations_real[rotationBegin + 1].ptrValidUpto(sHi);
			float* r_imag1 = rotations_imag[rotationBegin + 1].ptrValidUpto(sHi); 
			float* r_real2 = rotations_real[rotationBegin + 2].ptrValidUpto(sHi);
			float* r_imag2 = rotations_imag[rotationBegin + 2].ptrValidUpto(sHi); 
			assert(U64(r_real1) % 64 == 0); 
			assert(U64(r_imag1) % 64 == 0); 
			assert(U64(r_real2) % 64 == 0); 
			assert(U64(r_imag2) % 64 == 0); 

			HEADSTEP_SHARED(0) HEADSTEP_R(0,0)  HEADSTEP_R(1,0)  HEADSTEP_R(2,0)
			HEADSTEP_SHARED(1) HEADSTEP_R(0,1)  HEADSTEP_R(1,1)  HEADSTEP_R(2,1)
			HEADSTEP_SHARED(2) HEADSTEP_R(0,2)  HEADSTEP_R(1,2)  HEADSTEP_R(2,2)

			#pragma vector aligned
   		    #pragma ivdep
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP_BEFORE
				BODYSTEP_R(0,0)  BODYSTEP_R(1,0)  BODYSTEP_R(2,0)
				BODYSTEP_R(0,1)  BODYSTEP_R(1,1)  BODYSTEP_R(2,1)
				BODYSTEP_R(0,2)  BODYSTEP_R(1,2)  BODYSTEP_R(2,2)
				BODYSTEP_AFTER
			}

			TAILSTEP_R(0,0)	TAILSTEP_R(1,0)	TAILSTEP_R(2,0)
			TAILSTEP_R(0,1)	TAILSTEP_R(1,1)	TAILSTEP_R(2,1)
			TAILSTEP_R(0,2)	TAILSTEP_R(1,2)	TAILSTEP_R(2,2)

			return;
		}

		assert(rotationLen == 1);
		static_assert(maxTranslationIndexsCapacity <= 5, "Need more or fewer cases");
		switch (translationIndexsLen) {
		case 5: {
#ifndef MY_MACHINE_HAS_AVX512
			std::cerr << "UpdateOtherParams_Kernel_implementation using case 5 but not enough vector registers" << std::endl;
			std::cout << "UpdateOtherParams_Kernel_implementation using case 5 but not enough vector registers" << std::endl;
			EXIT_ABNORMALLY;
#endif
			HEADSTEP(0)
			HEADSTEP(1)
			HEADSTEP(2)
			HEADSTEP(3)
			HEADSTEP(4)

			#pragma vector aligned
   		    #pragma ivdep
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP_BEFORE
				BODYSTEP(0)
				BODYSTEP(1)
				BODYSTEP(2)
				BODYSTEP(3)
				BODYSTEP(4)
				BODYSTEP_AFTER
			}

			TAILSTEP(0)
			TAILSTEP(1)
			TAILSTEP(2)
			TAILSTEP(3)
			TAILSTEP(4)
		}; break;
		case 4: {
			HEADSTEP(0)
			HEADSTEP(1)
			HEADSTEP(2)
			HEADSTEP(3)
			#pragma vector aligned
   		    #pragma ivdep
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP_BEFORE
				BODYSTEP(0)
				BODYSTEP(1)
				BODYSTEP(2)
				BODYSTEP(3)
				BODYSTEP_AFTER
			}
			TAILSTEP(0)
			TAILSTEP(1)
			TAILSTEP(2)
			TAILSTEP(3)
		}; break;
		case 3: {
			HEADSTEP(0)
			HEADSTEP(1)
			HEADSTEP(2)
			#pragma vector aligned
   		    #pragma ivdep
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP_BEFORE
				BODYSTEP(0)
				BODYSTEP(1)
				BODYSTEP(2)
				BODYSTEP_AFTER
			}
			TAILSTEP(0)
			TAILSTEP(1)
			TAILSTEP(2)
		}; break;
		case 2: {
			HEADSTEP(0)
			HEADSTEP(1)
			#pragma vector aligned
   		    #pragma ivdep
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP_BEFORE
				BODYSTEP(0)
				BODYSTEP(1)
				BODYSTEP_AFTER
			}
			TAILSTEP(0)
			TAILSTEP(1)
		}; break;
		case 1: {
			HEADSTEP(0)
			#pragma vector aligned
   		    #pragma ivdep
			for (int n = sLo; n < sHi; n++) {
				BODYSTEP_BEFORE
				BODYSTEP(0)
				BODYSTEP_AFTER
			}
			TAILSTEP(0)
		}; break;
		default: assert(false);
		}	// switch

#undef TAILSTEP
#undef TAILSTEP_R
#undef BODYSTEP
#undef BODYSTEP_R
#undef BODYSTEP_AFTER
#undef BODYSTEP_BEFORE
#undef BODYSTEP_AFTER
#undef HEADSTEP
#undef HEADSTEP_R
#undef HEADSTEP_SHARED

		if (logging) {
			double partial_diff200;
			getRT(partial_diff200, 0,0);
			std::cout << "Diff2:" << partial_diff200 << std::endl;
			std::cout << std::endl;
		}
	}

	virtual void compute() {
		compute(64, 64, 512+256);
			// numbers found by experimention
	}

	virtual void compute(const S64 rotationStepProposed, const S64 translationStepProposed, const S64 sampleStep) {
		const int numberOfAcceptedSamples = sparseSampleIndexs.numberOfAcceptedSamples();
		for (int n = 0; n < numberOfAcceptedSamples; n++) {
			dense_wsum_sigma2_noise_float [n] = 0.0;
			dense_wsum_sigma2_noise_double[n] = 0.0;
		}
		dist.zero();
		Scheduler::compute(rotationStepProposed, translationStepProposed, sampleStep);
		for (int n = 0; n < numberOfAcceptedSamples; n++) {
			sparse_wsum_sigma2_noise[sparseSampleIndexs.sparseIndex(n)] += 
				dense_wsum_sigma2_noise_double[n] +
				dense_wsum_sigma2_noise_float [n];
		}
	}

	virtual void recompute(int rotation, int translation) {
		const S64 sampleStep = 256;
		for (int i = 0; i < 2; i++) {
			Plan plan(*this, rotation, rotation+1, translation, translation+1, sampleStep, false);
			auto const & d = dist(rotation, translation);
			dist.zero(rotation, translation);
			for (int sLo = 0; sLo < numberOfSamples; sLo += sampleStep) {
				auto sHi = std::min(numberOfSamples, sLo +  sampleStep);
				plan.executeMergedSteps(sLo, sHi);
			}
			std::cout << "Recomputed as kernel_diff2:" << d << std::endl;
		}
	}

	virtual bool getRT(
		double&	diff2,
		int		rotation,
		int		translation) {
		if (!mustDoIsSet(rotation, translation)) { diff2 = -666.777; return false; }
		diff2 = dist(rotation, translation);
		return true;
	}
	virtual bool getI(
		double&	diff2,
		int		i) {
		if (i >= appended.size()) { diff2 = -666.777; return false; }
		auto& a = appended[i];
		diff2 = dist(a.rotation, a.translation);
		return true;
	}
};


class UopkPool {
	struct Cell {
		UpdateOtherParams_Kernel_implementation* v;
		Cell() : v(nullptr) {}
		~Cell() { deallocate(); }
		bool matches(UpdateOtherParams_Kernel_implementation * v) { return this->v == v; }
		void allocate(UpdateOtherParams_Kernel_implementation* v) {
			assert(!this->v);
			this->v = v;
		}
		void deallocate() {
			sDelete(v);
		}
	};
	ListPerThreadTemplate<Cell> cellsPerThread;
public:
	UopkPool() {}
	UpdateOtherParams_Kernel_implementation* allocate(
		int iter,
		int numberOfRotations,
		int numberOfTranslations,
		int numberOfSamples,
		const int* Mresol_fine,
		double* wsum_sigma2_noise,
		int maxRotations,
		int maxTranslations) {
		auto & cell = *cellsPerThread.acquire();
		if (!cell.v || !cell.v->suitable(numberOfRotations, numberOfTranslations, numberOfSamples)) { 
			cell.deallocate();
			cell.allocate(UpdateOtherParams_Kernel_implementation::make(maxTranslationIndexsCapacity, numberOfRotations, numberOfTranslations, numberOfSamples));
		}
		cell.v->init(iter, Mresol_fine, wsum_sigma2_noise);
		return cell.v;
	}
	void deallocate(UpdateOtherParams_Kernel_implementation* v) {
		cellsPerThread.release(v);
	}
} uopk_pool;


UpdateOtherParams_Kernel* UpdateOtherParams_Kernel::acquire(
	int iter,
	int numberOfRotations,
	int numberOfTranslations,
	int numberOfSamples,
	const int* Mresol_fine,
	double* wsum_sigma2_noise,
	int maxRotations,
	int maxTranslations)
{
	return uopk_pool.allocate(iter, numberOfRotations, numberOfTranslations, numberOfSamples, Mresol_fine, wsum_sigma2_noise, maxRotations, maxTranslations);
}

void UpdateOtherParams_Kernel_implementation::release() {
	Scheduler::fini();
	uopk_pool.deallocate(this);
}

}	// namespace



namespace Map2dOptimizer_Kernel {

static S64 hash(S64 a, S64 b, S64 c, S64 d, S64 e) {
	S64 h = ((((a*123+b)*45+c)*67)+d*89)+e;
	return h ^ (h<<31);
}

static bool randomChoice(S64 hash, S64 inverseDensity) {
	return (hash % inverseDensity) == 0;
}

void unitTestCorrectness() {
	std::cout << "Map2dOptimizer_Kernel::unitTestCorrectness begin" << std::endl;
	static const int numberOfRotations		= 80;
	static const int numberOfTranslations	= 80;
	static const int numberOfSamples		= 64*3;

	auto weights = vNew(double,numberOfSamples);
	for (int i = 0; i < numberOfSamples; i++) { weights[i] = float(i); }
	
	std::vector<double*> translations_real(numberOfTranslations);
	std::vector<double*> translations_imag(numberOfTranslations);
	std::vector<double*> rotations_real   (numberOfRotations);
	std::vector<double*> rotations_imag   (numberOfRotations);
	{
		auto makeIons = [&](std::vector<double*>&v) { for (int i = 0; i < v.size(); i++) v[i] = mallocAlignedDoubleZeroed(numberOfSamples+1)+1; };	// deliberately misalign
		makeIons(translations_real);
		makeIons(translations_imag);
		makeIons(rotations_real);
		makeIons(rotations_imag);
	}

	bool errors(false);
	for (int inverseDensity = 1; inverseDensity < numberOfRotations*numberOfTranslations; inverseDensity *= 10) {
	    for (int oneIndex = 0; oneIndex < numberOfSamples; oneIndex = 2*oneIndex+1) {
	    	std::cout << "Map2dOptimizer_Kernel::unitTest oneIndex:" << oneIndex << std::endl;
	    	U64 count(0);
	    	U64 interestingCount(16);
	    	U64 nonZeroDistances(0);
	    
			// Give the pooling a chance to mess up...
			for (int tries = 0; tries < 2; tries++) {

		    	// Need to try all the positions
		    	for (int rotationIsNZ = 0; !errors && rotationIsNZ < 2; rotationIsNZ++) {
				    for (int imagIsNZ = 0; !errors && imagIsNZ < 2; imagIsNZ++) {
	    				for (int eltNZ = 0; eltNZ < (rotationIsNZ ? numberOfRotations : numberOfTranslations); eltNZ++) {
	    
	    				    for (int translationIndexsCapacity = maxTranslationIndexsCapacity; !errors && translationIndexsCapacity >= 1; translationIndexsCapacity--) {
	    	    	    	
								count++;
	    					    if (count == interestingCount) {
	    					    	std::cout << " count:" << count 
							    	#define P(X) << " " << #X << ":" << X
	    					    		P(oneIndex) P(rotationIsNZ) P(imagIsNZ) P(eltNZ) P(translationIndexsCapacity) P(nonZeroDistances)
							    	#undef P
	    					    		<< std::endl; 
	    					    	interestingCount *= 2;
	    					    }
	    					    
	    	    	    	    auto k = GetAllSquaredDifferences_Kernel::acquire(0,numberOfRotations,numberOfTranslations,numberOfSamples,weights);
	    	    	    	    
	    					    double* changed = NULL;
	    	    	    	    switch (rotationIsNZ*2+imagIsNZ) {
	    	    	    	    case 0: changed = translations_real[eltNZ]; break;
	    	    	    	    case 1: changed = translations_imag[eltNZ]; break;
	    	    	    	    case 2: changed = rotations_real   [eltNZ]; break;
	    	    	    	    case 3: changed = rotations_imag   [eltNZ]; break;
	    	    	    	    }
	    	    			    changed[oneIndex] = 1.0;
	    					    
	    	    	    	    for (auto r = 0; !errors && r < numberOfRotations; r++) {
	    	    	    	    	for (auto t = 0; !errors && t < numberOfTranslations; t++) {
							    		if (randomChoice(hash(rotationIsNZ,imagIsNZ,eltNZ,r,t),inverseDensity)) {
	    	    	    	    			k->append(r,t,r*numberOfTranslations+t,rotations_real[r],rotations_imag[r],translations_real[t],translations_imag[t],false);
							    		}
	    	    	    	    	}
	    	    	    	    }
	    	    	    	    k->compute();
	    					    
	    	    	    	    for (auto r = 0; !errors && r < numberOfRotations; r++) {
	    	    	    	    	for (auto t = 0; t < numberOfTranslations; t++) {
	    	    	    	    		double distance;
	    	    	    	    		bool didRT = k->getRT(distance, r, t);
	    	    	    	    		if (!didRT) {
	    	    	    	    			continue;
	    	    	    	    		}
	    	    	    	    		double expectedDistance = weights[oneIndex] * (rotationIsNZ ? (r == eltNZ) : (t == eltNZ) );
	    					    		if (expectedDistance != 0) nonZeroDistances++;
	    	    	    	    		if (!nearEnoughTemplate(distance,expectedDistance)) {
	    					    			// Because the kernel distributes the sqrt of the weight across the factors,
	    					    			// it may not be exactly equal
	    	    	    	    			std::cout << "Wrong distance:" << distance 
	    					    				<< " instead of  expectedDistance:" << expectedDistance 
	    					    				<< " (difference:" << distance - expectedDistance << ")"
	    					    				<< " at r:" << r << " t:" << t << std::endl;
	    					    			k->compute();
	    	    	    	    			k->recompute(r,t);
	    	    	    	    			errors = true;
											EXIT_ABNORMALLY;
	    	    	    	    		}
	    	    	    	    	}
	    	    	    	    }
	    	    	    	    
	    					    changed[oneIndex] = 0.0;
	    					    	// Since the captures are deferred until during the compute, can't restore until after the compute
	    					    	// and since used in the test above, can't undo until here
	    	    	    	    
	    	    	    	    k->release();
							}
	    	    	    }
	    			}
	    	    }
	    	}
	    }
	}
	{
		auto freeIons = [&](std::vector<double*>&v) { 
			for (int i = 0; i < v.size(); i++) { auto p = v[i]-1; aFree(p); }
		};
		freeIons(translations_real);
		freeIons(translations_imag);
		freeIons(rotations_real);
		freeIons(rotations_imag);
	}
	vDelete(weights);
	std::cout << "Map2dOptimizer_Kernel::unitTestCorrectness ended" << std::endl;
}

class FindMinInBrick {
	int xLo; int xHi; int yLo; int yHi; int zLo; int zHi;
	double d[125];
	int computing;
public:
	FindMinInBrick(int xLo, int xHi, int yLo, int yHi, int zLo, int zHi) {
		init(xLo, xHi, yLo, yHi, zLo, zHi);
	}
	void init(int xLo, int xHi, int yLo, int yHi, int zLo, int zHi) {
		this->xLo = xLo; this->xHi=xHi; this->yLo = yLo; this->yHi = yHi; this->zLo = zLo; this->zHi = zHi;
		for (int i = 0; i < 125; i++) d[i] = -1.0;	// not known
		computing = 0;								// next one to examine 
	}
	int x() { int p = computing/ 1%5; return xLo + (xHi - xLo)*p/4; }
	int y() { int p = computing/ 5%5; return yLo + (yHi - yLo)*p/4; }
	int z() { int p = computing/25%5; return zLo + (zHi - zLo)*p/4; }
	void sample(double v) { 
		d[computing] = v; 
		int smallest = 0;
		for (computing = 0; computing < 125; computing++) {
			if (d[computing] < 0) return;
			if (d[computing] > d[smallest]) continue;
			smallest = computing;
		}
		// We have 125 values
		// We have the smallest value, which is going to be the center of a new search cube 
		// We have 27 values within that search cube, that could be used to get a 25% speedup of the search sometime
		computing = smallest;
		int x  = this->x(),   y = this->y(), z = this->z(); 
		int dx = (xHi-xLo)/2, dy = (yHi-yLo)/2, dz = (zHi-zLo)/2;
		init(
			std::max(xLo, x - dx), std::min(xHi, x+dx),
			std::max(yLo, y - dy), std::min(yHi, y+dy),
			std::max(zLo, z - dz), std::min(zHi, z+dz));
	}
};

void unitTestPerformanceWkr(std::ofstream & csv);
void unitTestPerformance() {
	const char* const csvName = "Map2dOptimizer_Kernel_unitTestPerformance.csv";
	std::ofstream csv(csvName);
	csv << "kernelName" << ", " << "rotationStepProposed" << ", " << "translationStepProposed" << ", " << "sampleStep" << ", " << "sumOfElapsed" << std::endl;

	std::cout << "Map2dOptimizer_Kernel::unitTestPerformance begin, writing " << csvName << std::endl;
	unitTestPerformanceWkr(csv);
	unitTestPerformanceWkr(csv);
	unitTestPerformanceWkr(csv);
	std::cout << "Map2dOptimizer_Kernel::unitTestPerformance ended" << std::endl;
}
void unitTestPerformanceWkr(std::ofstream & csv) {
	static const int numberOfRotations		= 100;
	static const int numberOfTranslations	= 100;
	static const int numberOfSamples		= 180*180;

	auto weights = mallocAlignedDoubleNotZeroed(numberOfSamples+1)+1;					// deliberately misalign to catch bugs that assume it is aligned
	auto ctf     = mallocAlignedDoubleNotZeroed(numberOfSamples+1)+1;					// deliberately misalign to catch bugs that assume it is aligned
	for (int i = 0; i < numberOfSamples; i++) { ctf[i] = (weights[i] = float(i % 10000))/10000.0f; }

    const int current_Fsize2 = numberOfSamples*(numberOfSamples/2+1);
    auto Mresol_fine = mallocCacheAligned(int,current_Fsize2+1)+1;				// deliberately misalign to catch bugs that assume it is aligned

	auto tr = mallocAlignedDoubleNotZeroed(numberOfSamples+1)+1;
	auto ti = mallocAlignedDoubleNotZeroed(numberOfSamples+1)+1;
	auto rr = mallocAlignedDoubleNotZeroed(numberOfSamples+1)+1;
	auto ri = mallocAlignedDoubleNotZeroed(numberOfSamples+1)+1;
	for (int i = 0; i < numberOfSamples; i++) {
		tr[i] = ti[i] = rr[i] = ri[i] = double(i);
	}

	for (int kernel = 0; kernel<3; kernel++) {

		const char* kernelName = NULL;
		switch (kernel) {
		case 0:
			kernelName = "GetAllSquaredDifferences_Kernel_implementation";
			break;
		case 1:
			kernelName = "BackProjection_Kernel_implementation";
			break;
		case 2:
			kernelName = "UpdateOtherParams_Kernel_implementation";
			break;
		}
		if (!kernelName) break;
		std::cout << "Map2dOptimizer_Kernel::unitTestPerformance " << kernelName << std::endl;

		bool errors(false);
		size_t count(0);
		for (int inverseDensity = numberOfRotations*numberOfTranslations/8; inverseDensity < numberOfRotations*numberOfTranslations; inverseDensity *= 2) {
			std::cout << "Map2dOptimizer_Kernel::unitTestPerformance inverseDensity:" << inverseDensity << std::endl;

			double minElapsedTime = 999999.0;

			FindMinInBrick findMinInBrick(
				16/8,   64/8,						// Only consider byte-sized
				16/8,   64/8,						// Only consider byte-sized
				128/16, 1024/16);					// Only consider vector registers being full  64 bytes/sizeof(float)=16

			int    timedTrials  = 0;
			double sumOfElapsed = 0;

			for (int trials = 0; trials < 10000; trials++) {
				auto startTime = dtime();

				const S64 rotationStepProposed    = findMinInBrick.x() * 8;
				const S64 translationStepProposed = findMinInBrick.y() * 8;
				const S64 sampleStep              = findMinInBrick.z() * 16;

				RotationState rotationState(1,numberOfRotations);
				KernelBase* k = NULL;
				GetAllSquaredDifferences_Kernel* gsd_k = NULL;
				BackProjection_Kernel*			 bp_k  = NULL;
				UpdateOtherParams_Kernel*		 uop_k = NULL;
				switch (kernel) {
				case 0: {
					auto ki = GetAllSquaredDifferences_Kernel::acquire(0,numberOfRotations,numberOfTranslations,numberOfSamples,weights);
					k = gsd_k = ki;
				}	break;
				case 1: {
					auto ki = BackProjection_Kernel::acquire(0, rotationState, 0, numberOfTranslations, numberOfSamples);
					k = bp_k = ki;
				}	break;
				case 2: {
					auto ki = UpdateOtherParams_Kernel::acquire(0, numberOfRotations,numberOfTranslations,numberOfSamples,Mresol_fine,weights,numberOfRotations,numberOfTranslations);
					k = uop_k = ki;
				}	break;
				}
				    	    
				for (auto r = 0; !errors && r < numberOfRotations; r++) {
					for (auto t = 0; !errors && t < numberOfTranslations; t++) {
						if (!randomChoice(hash(trials,r,t,t*t,r^t),inverseDensity)) continue;
						switch (kernel) {
						case 0:
							gsd_k->append(r,t,r*numberOfTranslations+t,rr,ri,tr,ti,false);
							break;
						case 1:
							bp_k->append(r,t,rr,ri,tr,ti,weights,ctf,2.0);
							break;
						case 2:
							uop_k->append(r,t,rr,ri,tr,ti,2.0,false);
							break;
						}
				    }
				}
				k->compute(rotationStepProposed, translationStepProposed, sampleStep);
				
				switch (kernel) {
				case 0:
					for (auto r = 0; !errors && r < numberOfRotations; r++) {
						for (auto t = 0; t < numberOfTranslations; t++) {
							double distance;
					    	if (gsd_k->getRT(distance, r, t)) count++;
					    }
					}
					gsd_k->release();
					break;
				case 1:
					bp_k->releaseCapturedOutputs();
					count++;
					bp_k->release();
					break;
				case 2:
					for (int i = 0;; i++) {
						double distance;
						if (!uop_k->getI(distance,i)) break;
						count++;
					}
					uop_k->release();
					break;
				}

				auto elapsedTime = dtime() - startTime;

				sumOfElapsed += elapsedTime;
				if (timedTrials++ == 10) {
					minElapsedTime = std::min(minElapsedTime, sumOfElapsed);
					findMinInBrick.sample(sumOfElapsed);
					csv << kernelName << ", " << rotationStepProposed << ", " << translationStepProposed << ", " << sampleStep << ", " << sumOfElapsed << std::endl;
					sumOfElapsed = 0;
					timedTrials  = 0;
				}
			}

			const S64 rotationStepProposed    = findMinInBrick.x() * 8;
			const S64 translationStepProposed = findMinInBrick.y() * 8;
			const S64 sampleStep              = findMinInBrick.z() * 16;
			findMinInBrick.sample(sumOfElapsed);
			std::cout   << "  rp:" << std::setw(3) << rotationStepProposed
						<< "  tp:" << std::setw(3) << translationStepProposed
						<< "  ss:" << std::setw(3) << sampleStep
						<< "  et:" << minElapsedTime
						<< " for " << kernelName << " inverseDensity:" << inverseDensity
						<< std::endl;
		}
		std::cout << "Did " << count << std::endl;
	}
	{ auto p = (ri-1);			 aFree(p); }
	{ auto p = (rr-1);			 aFree(p); }
	{ auto p = (ti-1);			 aFree(p); }
	{ auto p = (tr-1);			 aFree(p); }
	{ auto p = (weights-1);		 aFree(p); }
	{ auto p = (Mresol_fine-1);	 aFree(p); }
	{ auto p = (ctf-1);			 aFree(p); }
}


};	// namespace Map2dOptimizer_Kernel
