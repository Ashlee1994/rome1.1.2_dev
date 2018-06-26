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

#pragma once
#include "./map2d_optimizer.h"


namespace Map2dOptimizer_Kernel {
	void setKernel3x3(bool to);

// For debugging the BackProjection_Kernel
//
// Each class has exp_nr_rot*exp_nr_over_rot rotated images associated with it
//		These are in the exp_Frefs_Rot_{real,imag} and exp_Fweights
//
// For each class and image
//		Start in RS_unbuffered state
//		Move to RS_zeroedBuffer when they are assigned a buffer
//			Move to RS_captureNeeded just before they are appended to the kernel
//				Even if has been given to the kernel already, so delta1 and capturedForSumming lead to delta2
//			Move to RS_capturedForSumming just after the kernel has captured the info
//			The kernel runs
//			The kernel returns
//			Move to RS_addPending as the kernel is deleted
//		Move to RS_unbuffered
//			If RS_addPending then the buffer is summed into its exp_Fref
//
class RotationState {
public:
	enum RS {			
		// This is the phases that BackProjection_Kernel pushes a exp_Frefs_Rot through
		//
		RS_unbuffered,			// Starts here
		RS_zeroedBuffer,		// the shifted weighted numbers might need to be accumulated
		RS_captured,			// the accumulator is captured by a kernel
		RS_addPending,			// the accumulator is released by a kernel but holds the value to be added
		RS_added				// the accumulator has been added
	};

	struct CurrentRotationState {
		RS    rs;
		void* mainAddr;			// in exp_Frefs_Rot
		void* addAddr;			// in the accumulator
		CurrentRotationState() 
		  :	rs(RS_unbuffered), 
			mainAddr(0),
			addAddr(0)
		{}
#if !defined(NDEBUG)
		~CurrentRotationState() {
			assert(rs == RS_unbuffered || rs == RS_added);
		}
#endif
	};

	const int numberOfRotations;

	RotationState(const int numberOfClasses, const int numberOfRotations) 
	  : numberOfRotations(numberOfRotations),
		currentRotationStates(numberOfClasses*numberOfRotations)
	{}

	RS getRotationState(
		int   iclass,
		int   rotation) const;

	void setRotationState(
		int   iclass,
		int   rotation,
		RS    to,
		void* mainAddr,
		void* addAddr);

private:
	std::vector<CurrentRotationState> currentRotationStates;
};


class KernelBase {
public:
	virtual ~KernelBase() {}

	virtual void compute() = 0;
	virtual void compute(const __int64 rotationStepProposed, const __int64 translationStepProposed, const __int64 sampleStep) = 0;
};


class BackProjection_Kernel : public KernelBase {
public:
	virtual void release() = 0;

	static BackProjection_Kernel* acquire(
		int iter,
		RotationState& rotationState,
		int iclass,
		int numberOfTranslations,
		int numberOfSamples);

	virtual void append(
		int rotation,
		int translation,
		      double* rotation_real,          double* rotation_imag,
		const double* translation_real, const double* translation_imag,
		      double* rotation_weight,        const double* ctf,
			  double  weight) = 0;

	virtual void releaseCapturedOutputs() = 0;
protected:
	~BackProjection_Kernel() {}		// call release instead
};


class GetAllSquaredDifferences_Kernel : public KernelBase {
public:
	virtual void release() = 0;

	static GetAllSquaredDifferences_Kernel* acquire(
		int iter,
		int numberOfRotations,
		int numberOfTranslations,
		int numberOfSamples,
		const double* weights);

	virtual void append(
		int rotation,
		int translation,
		int rot_trans_over,
		const double* rotation_real,    const double* rotation_imag,
		const double* translation_real, const double* translation_imag,
		bool interesting) = 0;

	virtual bool getRT(
		double&	distance,
		int		rotation,
		int		translation) = 0;

	virtual bool getI(
		double&	distance,
		int&    rot_trans_over,
		int		i) = 0;

	virtual void recompute(int rotation, int translation) = 0;	// used for debugging
protected:
	~GetAllSquaredDifferences_Kernel() {}		// call release instead
};

class UpdateOtherParams_Kernel : public KernelBase {
public:
	virtual void release() = 0;

	static UpdateOtherParams_Kernel* acquire(
		int iter,
		int numberOfRotations,
		int numberOfTranslations,
		int numberOfSamples,
		const int* Mresol_fine,
		double* wsum_sigma2_noise,
		const int maxRotations,
		const int maxTranslations);

	virtual void append(
		int rotation,
		int translation,
		const double* rotation_real,    const double* rotation_imag,
		const double* translation_real, const double* translation_imag, 
		double weight,
		bool   interesting) = 0;

	virtual bool getRT(
		double&	distance,
		int		rotation,
		int		translation) = 0;

	virtual bool getI(
		double&	distance,
		int		i) = 0;

	virtual void recompute(int rotation, int translation) = 0;	// used for debugging

protected:
	~UpdateOtherParams_Kernel() {}		// call release instead
};
    
void unitTestCorrectness();
void unitTestPerformance();

}
