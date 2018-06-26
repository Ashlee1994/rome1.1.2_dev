/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
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
#define MAP2DOPTIMIZER_NEW_BEING_COMPILED
#include "./map2d_optimizer_new.h"

#include "statusTracer.h"


// Compile-time options
//
namespace Map2dOptimizer_new {

	static const bool do_map         = true;
	static const bool do_zero_mask   = true;
	static const int  padding_factor = 2;

	//#define LOADBALANCEANALYZER
	//#define CHECK_FOR_NANS_ETC
	//#define FAKEMPI
}

//
namespace Map2dOptimizer_new {
    
    using namespace MapOptimizer_base_new;
#include "./tandemExecution_phaseAlgorithms_shared_source_code.h"
    
#define exp_nr_rot exp_nr_rot_NOT_USED
#define exp_nr_psi exp_nr_psi_NOT_USED
    
#define SEP
#define ELT(T, N) T N;
    MAP2DOPTIMIZER_STATIC_SCALARS
#undef ELT
#undef SEP
    
#undef exp_nr_rot
#undef exp_nr_psi
    
    static void ml_new_waypoint(const char* message,bool tag) {
        if (MapOptimizer_base_new::comparing_versions) {
            iterateAlgorithm(NULL)->waypoint(message);
        }
    }
}

// Debugging support
//
namespace Map2dOptimizer_new {

	static auto & emit_test_output          = MapOptimizer_base_new::emit_test_output;
	static auto & emit_test_output_prolific = MapOptimizer_base_new::emit_test_output_prolific;

	// This code only has the optimized versions
	// so that it is more easily enhanced
	//
	static std::ostream& testos() {
		static std::ofstream* ofp = nullptr;
		if (!ofp) {
			std::string ofn = 
				"../../testdata_2d/data8_160/emit_test_output_";		// TODO from command line
			if (emit_test_output_prolific()) ofn += "_prolific";
			ofn += "_opt";
			ofn += "_new";
			ofn += "_stab.txt";
			ofp = ofstreamCheckingCreated::make(ofn);
			*ofp << std::setprecision(4);
			std::cout << "Writing ~~TEST output to " << ofn << std::endl;
		}
		return *ofp;
	}

	static void ml_original_waypoint(const char* message) {
		if (MapOptimizer_base_new::comparing_versions) {
			iterateAlgorithm(NULL)->waypoint(message);
		}
	}
}



// For testing purposes some versions can write out some intermediate results to a file
// that can be compared to look for changes
//
namespace Map2dOptimizer_new {

	// Control the output of text describing intermediate results
	// and useful for comparing implementations
	//
	static const bool emit_test_output_possible = 
		true;												// change this as you desire
	
	bool write_to_all_versions_agree_emit_test_output_possible = 
		(MapOptimizer_base_new::all_versions_agree_emit_test_output_possible &= 
			emit_test_output_possible);
}

namespace Map2dOptimizer_new {
    using namespace MapOptimizer_base_new;
#include "./tandemExecution_phaseAlgorithms_shared_source_code.h"
}

//===============================================================================================
// Context
//
namespace Map2dOptimizer_new {
#include "./tandemExecution_phaseAlgorithms_shared_source_code.h"

	// Sizing
	//		This type is not used yet
	//		I am using it to document the various lengths
	//
	class ArrayLengths {
		int _ori_size;
		int _current_size;
		int _coarse_size;
	public:
		void set_ori_size(int to) { _ori_size = to; }
		int ori_size		  () const { return _ori_size;						}			// side dim of original images
		int pixelsPerImageSide() const { return ori_size();						}
		int pixelsPerImage    () const { return square(pixelsPerImageSide());	}

		void set_current_size(int to) { assert(to <= _ori_size); _current_size = to; }		// side dim of lower resolution version of images
		int current_size      () const { return _current_size; }
		int current_Fsize     () const { return current_size() / 2 + 1;			}
		int current_Fsize2    () const { return current_size()*current_Fsize();	}

		void set_coarse_size  (int to) { assert(to <= _current_size); _coarse_size = to; }	// side dim of even smaller images
		int coarse_size       () const { return _coarse_size; }
		int coarse_Fsize      () const { return coarse_size() / 2 + 1;			}
		int coarse_Fsize2     () const { return coarse_size()*coarse_Fsize();	}

		int ori_Fsize		  () const { return (_ori_size / 2 + 1);					 }	// number of shells
		int numberOfShells    () const { return ori_Fsize();							 }
	};

	// There are lots of arrays of vectors of various...
	// where the vector's length is an expression that is not computed when the object is first made.
	// These macros generate the deferred sizing of the vectors.
	//
#define RUNTIME_SIZED_ARRAY_CLASSES_0 \
    /* CLASSNAME,       BASECLASS,       BASE_COMPONENT,    SIZE_COMPUTABLE,                        LENGTH */                                                       \
    ELT(ImageDouble,    VectorOfDouble,  double,            ori_size != 0,                          ori_size*ori_size                                       )   SEP \
    ELT(ImageFloat,     VectorOfFloat,   float,             ori_size != 0,                          ori_size*ori_size                                       )   SEP \
    ELT(ImageRowDouble, VectorOfDouble,  double,            ori_size != 0,                          ori_size                                                )   SEP \
    ELT(ImageRowFloat,  VectorOfFloat,   float,             ori_size != 0,                          ori_size                                                )   SEP \
    ELT(FimgsData,      VectorOfDouble,  double,            current_size != 0,                      current_size*(current_size/2 + 1)                       )   SEP \
    ELT(PerClassDouble, VectorOfDouble,  double,            nr_classes != 0,                        nr_classes                                              )   SEP \
    ELT(PerImageDouble, VectorOfDouble,  double,            nr_pool != 0,                           nr_pool                                                 )   SEP \
    ELT(PerPsiDouble,   VectorOfDouble,  double,            expCfg.nr_rot_valid(),                  sampling2d.NrPsi(adaptive_oversampling)                 )   SEP \
    ELT(PerTransDouble, VectorOfDouble,  double,            expCfg.nr_trans_valid(),                sampling2d.NrTrans(adaptive_oversampling)               )   SEP \
    ELT(ModelData,      VectorOfDouble,  double,            ori_size != 0,                          ori_size/2+1                                            )   SEP \
    ELT(NrPsiChar,      VectorOfChar,    char,              expCfg.nr_rot_valid(),                  expCfg.nr_rot()                                         )		\
    // end of macro

#define RUNTIME_SIZED_ARRAY_CLASSES_1 \
    ELT(RotsFimgsData,  VectorOfFimgsData, double,          expCfg.nr_psi_adaptive_oversampling_valid(),    expCfg.nr_psi_adaptive_oversampling()			)   SEP \
    ELT(TransFimgsData, VectorOfFimgsData, double,          expCfg.nr_trans_adaptive_oversampling_valid(),  expCfg.nr_trans_adaptive_oversampling()			)		\
    // end of macro

#define RUNTIME_SIZED_ARRAY_CLASSES \
    RUNTIME_SIZED_ARRAY_CLASSES_0 SEP   \
    RUNTIME_SIZED_ARRAY_CLASSES_1       // end of macro

#define SEP
#define ELT(CLASSNAME, BASECLASS, BASE_COMPONENT, SIZE_COMPUTABLE, LENGTH)			\
    class CLASSNAME : public BASECLASS {											\
        inline bool isInitable();													\
        inline int  length();														\
    public:																			\
        CLASSNAME() { if (isInitable()) init(); }									\
        void init() { BASECLASS::init(length()); }									\
        void clear() { BASECLASS::init(0); }										\
        void initWith(double with) { init(); fill(with); }  						\
	    };																				// end of macro
	RUNTIME_SIZED_ARRAY_CLASSES_0
		typedef VectorOfVector<FimgsData> VectorOfFimgsData;
	RUNTIME_SIZED_ARRAY_CLASSES_1
#undef ELT
#undef SEP

}


//===============================================================================================
// Phase 
//
namespace Map2dOptimizer_new {

	void set_iter(int to) {
		assert(to >= 0);
        iter = to;
	}

}


//===============================================================================================
// Phase 
//
namespace Map2dOptimizer_new {

	MetaDataTable metadata;
    
    MLModel mlModel;
    MAPModel mapModel;
    ParticleModel particleModel;
    
	void setupMap2dOptimizer() {
#ifdef DATA_STREAM
        global_data_stream.init(data_stream_out_fn, data_stream_in_fn);
        global_data_stream.foutInt(&nr_iter, 1, "nr_iter", __FILE__, __LINE__);
        global_data_stream.foutInt(&nr_classes, 1, "nr_classes", __FILE__, __LINE__);
        global_data_stream.foutDouble(&pixel_size, 1, "pixel_size", __FILE__, __LINE__);
        global_data_stream.foutInt(&random_seed, 1, "random_seed", __FILE__, __LINE__);
        global_data_stream.foutDouble(&ini_high, 1, "ini_high", __FILE__, __LINE__);
        global_data_stream.foutDouble(&tau2_fudge_factor, 1, "tau2_fudge_factor", __FILE__, __LINE__);
        global_data_stream.foutDouble(&particle_diameter, 1, "particle_diameter", __FILE__, __LINE__);
        global_data_stream.foutInt(&adaptive_oversampling, 1, "adaptive_oversampling", __FILE__, __LINE__);
        global_data_stream.foutInt(&sampling2d.healpix_order, 1, "sampling.healpix_order", __FILE__, __LINE__);
        global_data_stream.foutDouble(&sampling2d.psi_step, 1, "sampling.psi_step", __FILE__, __LINE__);
        global_data_stream.foutDouble(-91, "sampling.limit_tilt", __FILE__, __LINE__);
        global_data_stream.foutDouble(&offset_range, 1, "sampling.offset_range", __FILE__, __LINE__);
        global_data_stream.foutDouble(&offset_step, 1, "sampling.offset_step", __FILE__, __LINE__);
        global_data_stream.foutDouble(0.5, "sampling.perturbation_factor", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
		// assert(iter >= 0);
		if (emit_test_output()) {
			testos() << "~~~~TEST OUTPUT: sigma2_offset " << MapOptimizer_base_new::sigma2_offset << std::endl;
		}
		metadata.readFromStar(star_fn);
        metadata.shuffle();
	}
}


// Context
//
namespace Map2dOptimizer_new {

	// The data for each of the locally stored images.
	// The images are square - equal height and width - but "side" and "size" are too easily confused
	//
	int const & ori_size = MapOptimizer_base_new::ori_size;

    void calspace(int current_size, int set_nr_pool);

	class Image : public ImageDouble {
	public:
		typedef double Elt;
		Image() : _width(ori_size), ImageDouble() { assert(_width*_width == size()); }
		int width() const { return _width; }
		Image(const Image& from) : _width(from._width), ImageDouble() {
			assert(_width*_width == size()); assert(_width == ori_size);
			*this = from;
		}
		void operator=(const Image& from) {
			assert(size() == from.size());
			ImageDouble* p = this;
			*p = from;
		}
	private:
		const int _width;
	};

	static void copy(ImageDouble& lhs, Image const& rhs) {
		auto   lhs_p = lhs.wptrAll();
		auto   rhs_p = rhs.rptrAll();
        ::copy(lhs_p, lhs.size(), rhs_p, rhs.size());
	}


	class Images : public std::vector<Image> {
	public:
		Images() {}
		Images(size_t size) : std::vector<Image>(size) {}
		static Images* make(size_t size) {
#include "./util_heap_undefs.h"
			return sNewA(Images, (size));
#include "./util_heap_defs.h"
		}
	};

	class ListOfImages : public ::ListOfImages {
		std::vector<ImageDouble> & v;
	public:
		ListOfImages(std::vector<ImageDouble> & v) : v(v) { }
		virtual ~ListOfImages() {}
		virtual int nr_images() { return v.size(); }
		virtual int imageSide() { assert(ori_size*ori_size == v[0].size()); return ori_size; }
		virtual const double* rImage(int i) { return v[i].rptrAll(); }
		virtual       double* wImage(int i) { return v[i].wptrAll(); }
		virtual       double* mImage(int i) { return v[i].mptrAll(); }
	};


	static void checkNormalize(Images& images, double particle_diameter, double pixel_size) {
		for (int iimage = 0; iimage < images.size(); iimage++) {
			auto & image = images[iimage];
			::checkNormalize(image.rptrAll(), image.width(), 1, particle_diameter, pixel_size);
		}
	}
}


// Phase 
//
namespace Map2dOptimizer_new {

	void readImages(Images*& imagesPtr, int nr_global_images, int first_local_image, int last_local_image) {

		// Read the local images from the file
		//
		imagesPtr = Images::make(last_local_image - first_local_image + 1);
		auto& images = *imagesPtr;

		const size_t bufLen = ori_size*ori_size;
		auto         buffer = mallocFloats(bufLen);
		for (int i = 0; i < nr_global_images; i++) {

			//TODO read each on the node that needs it, rather than broadcasting to all nodes
			//
			NODE0ONLY{
				MapOptimizer_base_new::readOneImage(bufLen, buffer, metadata[i]);
			}
#if defined(USEMPI)
			MPI::COMM_WORLD.Bcast(buffer,ori_size*ori_size,MPI::FLOAT,0);
#endif

			if (i < first_local_image || last_local_image < i) continue;

			auto & image = images[i - first_local_image];
			copyConverted(image.wptrAll(), image.size(), &buffer[0], bufLen);
		}
		Heap::freeFloats(buffer);
	}

}


// Context
//
namespace Map2dOptimizer_new {
	typedef VectorOfScalar<int> IntPerShell;			//
	typedef VectorOfScalar<int> IntPerCurrentFsize2;	//
	typedef VectorOfScalar<int> IntPerCoarseFsize2;		//
}


// Context
//
namespace Map2dOptimizer_new {

	// Status restored outside the iter loop
	// and so does not include any of the expectation data which is created and used inside the iter loop
	//
	class MyStatusTracer : public StatusTracer {
	public:
		MyStatusTracer() : StatusTracer() {}
		void appendDoubleVec(VectorOfDouble& data, int len, std::string what){
			StatusTracer::appendDoublePtr(data.wptrAll(), len, what);
		}
        void appendDouble(double& data, int len, std::string what){
            StatusTracer::appendDoublePtr(&data, len, what);
        }
        void appendFloatVec(VectorOfFloat& data, int len, std::string what, bool isCompareWithDouble){
            StatusTracer::appendFloatPtr(data.wptrAll(), len, what, isCompareWithDouble);
        }
        void appendFloat(float& data, int len, std::string what, bool isCompareWithDouble){
            StatusTracer::appendFloatPtr(&data, len, what, isCompareWithDouble);
        }
        void appendInt(int& data, int len, std::string what){
            StatusTracer::appendIntPtr(&data, len, what);
        }
		//void appendIntPtr(TBD data, int len, std::string what){
		//	StatusTracer::appendIntPtr();
		//}
	};
	MyStatusTracer statusTracer;

	static void setupStatusTracer();

	static void setupStatusTracer() {
		statusTracer.clear();
        for (int i = 0; i < mapModel.Irefs.size(); i++) {
#if defined(FLOAT_PRECISION)
            statusTracer.appendFloatPtr(mapModel.Irefs[i].wptr(), mapModel.Irefs[i].dimzyx, "mapmodel_Irefs", true);
#else
            statusTracer.appendDoublePtr(mapModel.Irefs[i].wptr(), mapModel.Irefs[i].dimzyx, "mapmodel_Irefs");
#endif
        }
        statusTracer.appendDouble(mapModel.current_resolution, 1, "mapModel_current_resolution");
        //
        for (int igroup = 0; igroup < mlModel.nr_groups; igroup++) {
#if defined(FLOAT_PRECISION)
            statusTracer.appendFloatVec(mlModel.sigma2_noise[igroup], mlModel.ori_Fsize, "mlmodel_sigma2_noise", true);
            statusTracer.appendFloatVec(mlModel.wsum_signal_product_spectra[igroup], mlModel.ori_Fsize, "mlmodel_wsum_signal_product_spectra", true);
            statusTracer.appendFloatVec(mlModel.wsum_reference_power_spectra[igroup], mlModel.ori_Fsize, "mlmodel_wsum_reference_power_spectra", true);
#else
            statusTracer.appendDoubleVec(mlModel.sigma2_noise[igroup],
                                         mlModel.ori_Fsize, "mlmodel_sigma2_noise_"+num2str(igroup,6));
            statusTracer.appendDoubleVec(mlModel.wsum_signal_product_spectra[igroup],
                                         mlModel.ori_Fsize, "mlmodel_wsum_signal_product_spectra_"+num2str(igroup,6));
            statusTracer.appendDoubleVec(mlModel.wsum_reference_power_spectra[igroup],
                                         mlModel.ori_Fsize, "mlmodel_wsum_reference_power_spectra_"+num2str(igroup,6));
#endif
        }
        for (int iclass = 0; iclass < mlModel.nr_classes; iclass++) {
#if defined(FLOAT_PRECISION)
            statusTracer.appendFloatVec(mlModel.tau2_class[iclass], mlModel.ori_Fsize, "mlmodel_tau2_class", true);
            statusTracer.appendFloatVec(mlModel.sigma2_class[iclass], mlModel.ori_Fsize, "mlmodel_sigma2_class", true);
            statusTracer.appendFloatVec(mlModel.data_vs_prior_class[iclass], mlModel.ori_Fsize, "mlmodel_data_vs_prior_class", true);
            statusTracer.appendFloatVec(mlModel.fsc_halves_class[iclass], mlModel.ori_Fsize, "mlmodel_fsc_halves_class", true);
            statusTracer.appendFloatVec(mlModel.pdf_direction[iclass], mlModel.nr_directions, "mlmodel_pdf_direction", true);
#else
            statusTracer.appendDoubleVec(mlModel.tau2_class[iclass],
                                         mlModel.ori_Fsize, "mlmodel_tau2_class_"+num2str(iclass,6));
            statusTracer.appendDoubleVec(mlModel.sigma2_class[iclass],
                                         mlModel.ori_Fsize, "mlmodel_sigma2_class_"+num2str(iclass,6));
            statusTracer.appendDoubleVec(mlModel.data_vs_prior_class[iclass],
                                         mlModel.ori_Fsize, "mlmodel_data_vs_prior_class_"+num2str(iclass,6));
            statusTracer.appendDoubleVec(mlModel.fsc_halves_class[iclass],
                                         mlModel.ori_Fsize, "mlmodel_fsc_halves_class_"+num2str(iclass,6));
            statusTracer.appendDoubleVec(mlModel.pdf_direction[iclass],
                                         mlModel.nr_directions, "mlmodel_pdf_direction_"+num2str(iclass,6));
#endif
        }
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatVec(mlModel.scale_correction, mlModel.nr_groups, "mlmodel_scale_correction", true);
        statusTracer.appendFloatVec(mlModel.prior_offsetx_class, mlModel.nr_classes, "mlmodel_prior_offsetx_class", true);
        statusTracer.appendFloatVec(mlModel.prior_offsety_class, mlModel.nr_classes, "mlmodel_prior_offsety_class", true);
        statusTracer.appendFloatVec(mlModel.pdf_class, mlModel.nr_classes, "mlmodel_pdf_class", true);
        statusTracer.appendFloat(mlModel.avg_norm_correction, 1, "mlmodel_avg_norm_correction", true);
        statusTracer.appendFloat(mlModel.ave_Pmax, 1, "mlmodel_ave_Pmax", true);
        statusTracer.appendFloat(mlModel.sigma2_offset, 1, "mlmodel_sigma2_offset", true);
#else
        statusTracer.appendDoubleVec(mlModel.scale_correction,
                                     mlModel.nr_groups, "mlmodel_scale_correction");
        statusTracer.appendDoubleVec(mlModel.prior_offsetx_class,
                                     mlModel.nr_classes, "mlmodel_prior_offsetx_class");
        statusTracer.appendDoubleVec(mlModel.prior_offsety_class,
                                     mlModel.nr_classes, "mlmodel_prior_offsety_class");
        statusTracer.appendDoubleVec(mlModel.pdf_class,
                                     mlModel.nr_classes, "mlmodel_pdf_class");
        statusTracer.appendDouble(mlModel.avg_norm_correction, 1, "mlmodel_avg_norm_correction");
        statusTracer.appendDouble(mlModel.ave_Pmax, 1, "mlmodel_ave_Pmax");
        statusTracer.appendDouble(mlModel.sigma2_offset, 1, "mlmodel_sigma2_offset");
#endif
        for (int igroup = 0; igroup < mlModel.nr_groups; igroup++) {
#if defined(FLOAT_PRECISION)
            statusTracer.appendFloatVec(mlModel.wsum_sigma2_noise[igroup], mlModel.ori_Fsize, "mlmodel_wsum_sigma2_noise", true);
#else
            statusTracer.appendDoubleVec(mlModel.wsum_sigma2_noise[igroup],
                                         mlModel.ori_Fsize, "mlmodel_wsum_sigma2_noise_"+num2str(igroup,6));
#endif
        }
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatVec(mlModel.wsum_sumw_group, mlModel.nr_groups, "mlmodel_wsum_sumw_group", true);
        statusTracer.appendFloatVec(mlModel.wsum_pdf_class, mlModel.nr_classes, "mlmodel_wsum_pdf_class", true);
        statusTracer.appendFloatVec(mlModel.wsum_prior_offsetx_class, mlModel.nr_classes, "mlmodel_wsum_prior_offsetx_class", true);
        statusTracer.appendFloatVec(mlModel.wsum_prior_offsety_class, mlModel.nr_classes, "mlmodel_wsum_prior_offsety_class", true);
#else
        statusTracer.appendDoubleVec(mlModel.wsum_sumw_group,
                                     mlModel.nr_groups, "mlmodel_wsum_sumw_group");
        statusTracer.appendDoubleVec(mlModel.wsum_pdf_class,
                                     mlModel.nr_classes, "mlmodel_wsum_pdf_class");
        statusTracer.appendDoubleVec(mlModel.wsum_prior_offsetx_class,
                                     mlModel.nr_classes, "mlmodel_wsum_prior_offsetx_class");
        statusTracer.appendDoubleVec(mlModel.wsum_prior_offsety_class,
                                     mlModel.nr_classes, "mlmodel_wsum_prior_offsety_class");
#endif
        for (int iclass = 0; iclass < mlModel.nr_classes; iclass++) {
#if defined(FLOAT_PRECISION)
            statusTracer.appendFloatVec(mlModel.wsum_pdf_direction[iclass], mlModel.nr_directions, "mlmodel_wsum_pdf_direction", true);
#else
            statusTracer.appendDoubleVec(mlModel.wsum_pdf_direction[iclass],
                                         mlModel.nr_directions, "mlmodel_wsum_pdf_direction"+num2str(iclass,6));
#endif
        }
        //
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloat(sampling2d.random_perturbation,1,"sampling3d.random_perturbation",true);
#else
        statusTracer.appendDouble(sampling2d.random_perturbation,1,"sampling3d.random_perturbation");
#endif
        //
        for (int iimage = 0; iimage < nr_global_images; iimage++) {
#if defined(FLOAT_PRECISION)
            statusTracer.appendFloatPtr(&metadata[iimage].NORM, 1, "metadata["+num2str(iimage,6)+"].NORM", true);
            statusTracer.appendFloatPtr(&metadata[iimage].XOFF, 1, "metadata["+num2str(iimage,6)+"].XOFF", true);
            statusTracer.appendFloatPtr(&metadata[iimage].YOFF, 1, "metadata["+num2str(iimage,6)+"].YOFF", true);
            statusTracer.appendFloatPtr(&metadata[iimage].PSI, 1, "metadata["+num2str(iimage,6)+"].PSI", true);
#else
            statusTracer.appendDouble(metadata[iimage].NORM, 1, "metadata["+num2str(iimage,6)+"].NORM");
            statusTracer.appendDouble(metadata[iimage].XOFF, 1, "metadata["+num2str(iimage,6)+"].XOFF");
            statusTracer.appendDouble(metadata[iimage].YOFF, 1, "metadata["+num2str(iimage,6)+"].YOFF");
            statusTracer.appendDouble(metadata[iimage].PSI, 1, "metadata["+num2str(iimage,6)+"].PSI");
#endif
            statusTracer.appendInt(metadata[iimage].CLASS, 1, "metadata["+num2str(iimage,6)+"].CLASS");
        }
	}


	// Expectation tiling support
	//
	class ExpectationConfig {
#define ExpectationConfigMEMBERS \
	ELT(StateMask_adaptive,	   nr_trans_adaptive_oversampling)	SEP \
	ELT(StateMask_adaptive,	   nr_psi_adaptive_oversampling)	SEP \
	ELT(StateMask_rotTran,     nr_rot_EQUAL_psi)				SEP \
	ELT(StateMask_rotTran,     nr_trans)						SEP \
    ELT(StateMask_ipass,       ipass)							SEP \
	ELT(StateMask_ipass,       current_oversampling)			SEP \
	ELT(StateMask_currentSize, current_size)					SEP \
    ELT(StateMask_overRotTran, nr_over_rot)						SEP \
	ELT(StateMask_overRotTran, nr_over_trans)					    \
	// end of macro
	public:
#define SEP
#define ELT(M,N) int N() const { assert(_state & M); return _##N; }  bool N##_valid() const { return (_state & M) != 0; }  int & _##N##_writableDontUseExceptInCompare;
		ExpectationConfigMEMBERS
			int nr_rot() const { return nr_rot_EQUAL_psi(); } bool nr_rot_valid() const { return nr_rot_EQUAL_psi_valid(); }
		int nr_psi() const { return nr_rot_EQUAL_psi(); } bool nr_psi_valid() const { return nr_rot_EQUAL_psi_valid(); }
#undef ELT
#undef SEP
		ExpectationConfig() : _state(State_uninited),
#define SEP ,
#define ELT(M,N) _##N(0), _##N##_writableDontUseExceptInCompare(_##N)
			ExpectationConfigMEMBERS
#undef ELT
#undef SEP
		{
		}
		void init(HealpixSampler& sampling2d) {
			_nr_trans						= sampling2d.NrTrans();
			_nr_rot_EQUAL_psi				= sampling2d.NrPsi();
			_nr_trans_adaptive_oversampling = sampling2d.NrTrans					   (MapOptimizer_base_new::adaptive_oversampling);
			_nr_psi_adaptive_oversampling   = sampling2d.NrPsi						   (MapOptimizer_base_new::adaptive_oversampling);
			_nr_over_rot					= sampling2d.oversamplingFactorOrientations(MapOptimizer_base_new::adaptive_oversampling);
			_nr_over_trans					= sampling2d.oversamplingFactorTranslations(MapOptimizer_base_new::adaptive_oversampling);
			_state |= StateMask_rotTran | StateMask_adaptive | StateMask_overRotTran;
		}
		void setCurrentSize(int current_size) {
			_current_size = current_size;
			_state |= StateMask_currentSize;
		}
		void setIpass(int ipass) {
			_ipass = ipass;
			_current_oversampling = (ipass == 0) ? 0 : MapOptimizer_base_new::adaptive_oversampling;
			_state |= StateMask_ipass;
		}
		void setIpass(HealpixSampler& sampling2d,
			int ipass,
			int current_size,
			int coarse_size) {
			// Use smaller images in the first pass,  larger ones in the second pass
			// Use coarse sampling in the first pass, oversampled on the second pass
			setIpass(ipass);
			setCurrentSize(((ipass == 0) && MapOptimizer_base_new::adaptive_oversampling) ? coarse_size : current_size);
			_nr_over_rot   = sampling2d.oversamplingFactorOrientations(_current_oversampling);
			_nr_over_trans = sampling2d.oversamplingFactorTranslations(_current_oversampling);
			_state |= StateMask_overRotTran;
		}

	protected:
		enum State { 
			State_uninited			=  0, 
			StateMask_rotTran		=  1,
			StateMask_adaptive		=  2, 
			StateMask_currentSize	=  4, 
			StateMask_ipass			=  8, 
			StateMask_overRotTran	= 16};
		int _state;
#define SEP 
#define ELT(M,N) int _##N;
		ExpectationConfigMEMBERS
#undef ELT
#undef SEP
	};


	class ExpectationImages {
#define ExpectationImagesMEMBERS \
		ELT(StateMask_classes, iclass_min)	SEP \
		ELT(StateMask_classes, iclass_max)	SEP \
		ELT(StateMask_images,  first_image)	SEP \
		ELT(StateMask_images,  last_image)	SEP \
		ELT(StateMask_images,  nr_images)		\
		// end of macro
	public:
#define SEP
#define ELT(M,N) int N() const { assert(_state & M); return _##N; }  bool N##_valid() const { return ((_state & M) != 0); }  int & _##N##_writableDontUseExceptInCompare;
		ExpectationImagesMEMBERS
#undef ELT
#undef SEP
			ExpectationImages() : _state(State_uninited),
#define SEP ,
#define ELT(M,N) _##N(0), _##N##_writableDontUseExceptInCompare(_##N)
			ExpectationImagesMEMBERS
#undef ELT
#undef SEP
		{}
		void init(int my_first_image, int my_last_image) {
			_first_image = my_first_image;
			_last_image = my_last_image;
			_nr_images = my_last_image - my_first_image + 1;
			_state |= StateMask_images;
		}
		void setClass(int iclass_min, int iclass_max) {
			_iclass_min = iclass_min;
			_iclass_max = iclass_max;
			_state |= StateMask_classes;
		}
	protected:
		enum State { State_uninited = 0, StateMask_images = 1, StateMask_classes = 2 };
		int _state;
#define SEP 
#define ELT(M,N) int _##N;
		ExpectationImagesMEMBERS
#undef ELT
#undef SEP
	};

	class Exp_FrefsAndWeight {
        enum Need {Need_Readonly, Need_Write1st, Need_Zeroed, Need_ReadWrite, Need_NoAccess};
    public:
        // These are three related globals that are conceptually zeroed early in an iteration, filled in during the iteration, and used throughout
        // However zeroing these is very expensive in the later iterations of the code
        // so instead we defer the zeroing until it is definitely needed, and even share the zeros when possible.
        //
        Exp_FrefsAndWeight() : numberOfSubsections(0), maxSubsectionLength(0), subsectionLength(0), ptrCapacity(0), _ptr_shared_zeroes(NULL) {
            for (int i = 0; i < 3; i++) {
                _ptrs_DONT_USE_UNLESS_ABSOLUTELY_NECESSARY[i] = NULL;
                _stages[i] = NULL;
            }
        }

        void save(Exp_FrefsAndWeight & into) {
            assert(into.numberOfSubsections == 0);
            into.init(numberOfSubsections, maxSubsectionLength);
            copy(into, *this);
        }

        void restore(Exp_FrefsAndWeight & from) {
            copy(*this, from);
            from.fini();
        }

        void init(int numberOfSubsections, int maxSubsectionLength) {
            this->numberOfSubsections = numberOfSubsections;
            this->maxSubsectionLength = maxSubsectionLength;
            this->subsectionLength    = maxSubsectionLength;
            this->ptrCapacity         = numberOfSubsections*maxSubsectionLength;
            for (int i = 0; i < 3; i++) {
                _ptrs_DONT_USE_UNLESS_ABSOLUTELY_NECESSARY[i] = Heap::allocDoubles(ptrCapacity, __FILE__, __LINE__);
                auto s = _stages[i] = vNew(Stage,numberOfSubsections);
                for (int j = 0; j < numberOfSubsections; j++) s[j] = Garbage;
            }
            _ptr_shared_zeroes = Heap::allocZeroedDoubles(maxSubsectionLength, __FILE__, __LINE__);
        }

        void fini() {
            for (int i = 0; i < 3; i++) {
                Heap::freeDoubles(_ptrs_DONT_USE_UNLESS_ABSOLUTELY_NECESSARY[i]);
                _ptrs_DONT_USE_UNLESS_ABSOLUTELY_NECESSARY[i] = NULL;
                vDelete(_stages[i]);
            }
            Heap::freeDoubles(_ptr_shared_zeroes);
            numberOfSubsections = 0;
            maxSubsectionLength = 0;
            subsectionLength    = 0;
        }

        void zero(int lineNumber) {
            TUNING_SCOPE_STEP(serial_zeroing_exp_Frefs_and_Fweight)
            const int ptrSize = numberOfSubsections*subsectionLength;
            assert(ptrCapacity >= ptrSize);
            assert(maxSubsectionLength >= subsectionLength);

            if (false)  // BEVIN
            #pragma omp critical
            {
                std::cerr << "Exp_FrefsAndWeight last used at ";
                for (int i = 0; i < 3; i++) std::cerr << _lastUsed[i] << " ";
                std::cerr << std::endl;

                std::cout << "Exp_FrefsAndWeight zeroing all subsections at " << lineNumber << std::endl;
            }

            if (subsectionLength == 0) subsectionLength = maxSubsectionLength;
                // Cope with the first call to zero() after init() but before any of the ptr calls

            for (int i = 0; i < 3; i++) {
                int alreadyZeroCount = 0;
                for (int j = 0; j < numberOfSubsections; j++) {
                    if (_stages[i][j] == ZeroDeferred || _stages[i][j] == ZeroDone) {
                        alreadyZeroCount++;
                        continue;
                    }
                    _stages[i][j] = ZeroDeferred;
                }
                if (false && alreadyZeroCount > 0) {
                    std::cerr << "Exp_FrefsAndWeight::zero[" << i << "] alreadyZeroCount:" << alreadyZeroCount << " out of " << numberOfSubsections << std::endl;
                }
            }

            subsectionLength = 0;   // This will be increased by the first called of the following
        }

        void zero(int lineNumber, int subsectionIndex, int subsectionLength) {
            for (int i = 0; i < 3; i++) {
                if (_stages[i][subsectionIndex] == ZeroDone) continue;
                if (false && subsectionIndex == 1)
                #pragma omp critical
                {
                    if (false) {    // BEVIN
                        std::cout << "Exp_FrefsAndWeight zeroing subsectionIndex " << subsectionIndex << " at line " << lineNumber << ".  Any subsection last used at line:" << _lastUsed[i] << std::endl;
                    }
                    static std::map<int,int> shown;
                    auto lineNumberLastUsed = _lastUsed[i];
                    if (shown.find(lineNumberLastUsed) == shown.end()) {
                        std::cerr << "Exp_FrefsAndWeight zeroing subsectionIndex " << subsectionIndex << " at line " << lineNumber << ".  Any subsection last used at line:" << _lastUsed[i] << std::endl;
                        shown.insert(std::make_pair(lineNumberLastUsed,0));
                    }
                }
                _stages[i][subsectionIndex] = ZeroDeferred;
				#if !defined(NDEBUG)
                    auto p = ptr(Need_ReadWrite, i, 0, subsectionIndex, subsectionLength);
                    _stages[i][subsectionIndex] = ZeroDeferred;
                    for (int n = 0; n < subsectionLength; n++) {
                        p[n] = -666.999;
                    }
                #endif
            }
        }

        double const * frefs_Rot_real_readonly      (int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_Readonly,  0, lineNumber, subsectionIndex, subsectionLength); }
        double const * frefs_Rot_imag_readonly      (int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_Readonly,  1, lineNumber, subsectionIndex, subsectionLength); }
        double const * fweight_Rot_readonly         (int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_Readonly,  2, lineNumber, subsectionIndex, subsectionLength); }
            // The data is being read but not written
            // So if it is zeros, those zeros will not be changed, and can be shared with other uses

        double*        frefs_Rot_real_write1st      (int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_Write1st,  0, lineNumber, subsectionIndex, subsectionLength); }
        double*        frefs_Rot_imag_write1st      (int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_Write1st,  1, lineNumber, subsectionIndex, subsectionLength); }
        double*        fweight_Rot_write1st         (int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_Write1st,  2, lineNumber, subsectionIndex, subsectionLength); }
            // The data is all being written, so any previous writes - such as zeroing it - can be ignored

        double*        frefs_Rot_real_writableZeroed(int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_Zeroed,    0, lineNumber, subsectionIndex, subsectionLength); }
        double*        frefs_Rot_imag_writableZeroed(int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_Zeroed,    1, lineNumber, subsectionIndex, subsectionLength); }
        double*        fweight_Rot_writableZeroed   (int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_Zeroed,    2, lineNumber, subsectionIndex, subsectionLength); }
            // The data is being updated from its initially zeroed state

        double*        frefs_Rot_real_readWrite     (int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_ReadWrite, 0, lineNumber, subsectionIndex, subsectionLength); }
        double*        frefs_Rot_imag_readWrite     (int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_ReadWrite, 1, lineNumber, subsectionIndex, subsectionLength); }
        double*        fweight_Rot_readWrite        (int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_ReadWrite, 2, lineNumber, subsectionIndex, subsectionLength); }
            // The data is being updated from its initially zeroed state

        double*        frefs_Rot_real_noAccess      (int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_NoAccess, 0, lineNumber, subsectionIndex, subsectionLength); }
        double*        frefs_Rot_imag_noAccess      (int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_NoAccess, 1, lineNumber, subsectionIndex, subsectionLength); }
        double*        fweight_Rot_noAccess         (int lineNumber, int subsectionIndex, int subsectionLength) { return ptr(Need_NoAccess, 2, lineNumber, subsectionIndex, subsectionLength); }
            // The data is not being updated but want the address of where it is writable

        bool           frefs_Rot_real_isZeroed      (int lineNumber, int subsectionIndex, int subsectionLength) { return isZeroed(0, subsectionIndex); }
        bool           frefs_Rot_imag_isZeroed      (int lineNumber, int subsectionIndex, int subsectionLength) { return isZeroed(1, subsectionIndex); }
        bool           fweight_Rot_isZeroed         (int lineNumber, int subsectionIndex, int subsectionLength) { return isZeroed(2, subsectionIndex); }
            // The data is being updated from its initially zeroed state

        int currentNumberOfSubsections() const { return numberOfSubsections; }
        int currentSubsectionLength()    const { return subsectionLength;    }

    private:
        int numberOfSubsections;
        int maxSubsectionLength;
        int ptrCapacity;
        double* _ptr_shared_zeroes;

        int subsectionLength;

        void copy(Exp_FrefsAndWeight & into, Exp_FrefsAndWeight & from) {
#define EQUAL(E) assert(from.E == E);
#define COPY(E)  E = from.E;
            EQUAL(numberOfSubsections)
            EQUAL(maxSubsectionLength)
            EQUAL(ptrCapacity        )
            COPY(subsectionLength    )
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < ptrCapacity;         j++) COPY(_ptrs_DONT_USE_UNLESS_ABSOLUTELY_NECESSARY[i][j]);
                for (int j = 0; j < numberOfSubsections; j++) COPY(_stages[i][j]);
            }
#undef COPY
#undef EQUAL
        }


        enum Stage{Garbage, ZeroDeferred, ZeroDone, Written}* _stages[3];

        bool isZeroed(int i, int subsectionIndex) {
            auto const stage = _stages[i][subsectionIndex];
            return stage == ZeroDeferred || stage == ZeroDone;
        }

        double* ptr(Need need, int i, int lineNumber, int subsectionIndex, int subsectionLength) {

            auto const oldStage = _stages[i][subsectionIndex];

            // Check legit and remember the ongoing subsectionLength
            assert(subsectionIndex  <  this->numberOfSubsections);
            if (this->subsectionLength == 0) {
                assert(subsectionLength <= maxSubsectionLength);
                this->subsectionLength = subsectionLength;
            } else {
                assert(this->subsectionLength == subsectionLength);
            }
            if (oldStage == Garbage
                ) {
                std::cerr << "Exp_FrefsAndWeight garbage accessed" << std::endl;
                EXIT_ABNORMALLY;
            }

            auto result = _ptrs_DONT_USE_UNLESS_ABSOLUTELY_NECESSARY[i] + subsectionIndex*subsectionLength;

            switch (need) {
            case Need_NoAccess:
                break;
            case Need_Readonly:
                if (oldStage == ZeroDone || oldStage == ZeroDeferred) {
                    result = _ptr_shared_zeroes;                // TODO memory protect them readonly
                    if (false) for (int n = 0; n < subsectionLength; n++) if (result[n] != 0.0) { std::cerr << "Exp_FrefsAndWeight zeroes not zero!" << std::endl; EXIT_ABNORMALLY; }
                } else if (oldStage != Written) {
                    std::cerr << "Exp_FrefsAndWeight asked for a readonly non-written section" << std::endl;
                    EXIT_ABNORMALLY;
                }
                break;
            case Need_Write1st:
                // Should not be assuming the zeroing has been done
                // since Need_Write1st is an intent to not read any of the current values
                #if !defined(NDEBUG)
                    for (int n = 0; n < subsectionLength; n++) result[n] = -6.66e88;
                #endif
                break;
            case Need_Zeroed:
                if (oldStage != ZeroDone) {
                    if (oldStage != ZeroDeferred) {
                        std::cerr << "Exp_FrefsAndWeight !ZeroDeferred, is " << _stages[i][subsectionIndex] << std::endl;
                        EXIT_ABNORMALLY;
                    }       // BEVIN
                    for (int n = 0; n < subsectionLength; n++) result[n] = 0.0;
                }
                break;
            case Need_ReadWrite:
                if (oldStage == ZeroDeferred) {
                    for (int n = 0; n < subsectionLength; n++) result[n] = 0.0;
                    _stages[i][subsectionIndex] = ZeroDone;
                }
                break;
            }

            if (need != Need_Readonly && need != Need_NoAccess) _stages[i][subsectionIndex] = Written;      // leave as ZeroDeferred or ZeroDone
            if (lineNumber != 0) _lastUsed[i] = lineNumber;

            if (false && subsectionIndex == 1)
            #pragma omp critical
            {
                std::cout << "Exp_FrefsAndWeight ptr subsectionIndex:" << subsectionIndex
                    << " changing stage from " << oldStage << " to " << _stages[i][subsectionIndex]
                    << " at line " << lineNumber << (lineNumber ? "" : " all subsections being zero()'ed")
                    << std::endl;
            }

            return result;
        }

        int     _lastUsed[3];
    public:
        double* _ptrs_DONT_USE_UNLESS_ABSOLUTELY_NECESSARY[3];  // DONT USE UNLESS ABSOLUTELY NECESSARY - IN THE COMPARE CODE BELOW
    };

	class Exp_Fimgs_shifted {
    public:
        static const bool debugging         = false;
        static const bool debuggingDeferred = debugging || false;   // BEVIN
        static const int  bufferCapacity    = 4;

        class Buffer {
			Exp_Fimgs_shifted* parent;
            const char* const name;
            const int lineNumber;
            const int exp_current_Fsize2;
            int len;
            int ishifts[bufferCapacity];
            int countAppends;
        public:
            Buffer(
				Exp_Fimgs_shifted* parent,
				const char* name, int lineNumber, int exp_current_Fsize2) : parent(parent), name(name), lineNumber(lineNumber), exp_current_Fsize2(exp_current_Fsize2),len(0),countAppends(0) {}
            ~Buffer() {
                static DoSomePerIter doSomePerIter;
                doSomePerIter.note(debuggingDeferred, iter, [&](int count){
                    std::cerr << "Exp_Fimgs_shifted::Buffer::~Buffer count:" << count << " countAppends:" << countAppends << std::endl;
                });
                flush(name, lineNumber);
            }
            void append(const char* why, int lineNumber, int ishift) {
                countAppends++;
                static DoSomePerIter doSomePerIter;
                doSomePerIter.note(debuggingDeferred, iter, [&](int count) {
                    std::cerr << "Exp_Fimgs_shifted::Buffer::append count:" << count << " countAppends:"<< countAppends << std::endl;
                });
                if (len == bufferCapacity) flush(why, lineNumber);
                ishifts[len++] = ishift;
            }
        private:
            void flush(const char* why, int lineNumber);
        };

        Exp_Fimgs_shifted() 
		  : realPtr(NULL), imagPtr(NULL), 
		    max_nr_images(0), max_ishift(0), max_current_Fsize2(0), 
			state(NULL), unusedDeferred(0), unusedWritten(0), usedWritten(0),
            captured_Mctf_invsigma2(NULL),
            deferred(NULL) {}

        bool           efsDeferred          (int ishift, int exp_current_Fsize2) { return efsDeferred(ishift); }
        void           efsUndefer           (const char* why, int lineNumber, int ishift, int exp_current_Fsize2, bool isRead, bool isWrite) {
                                                                                            efsUndeferWkr(why, lineNumber, ishift, exp_current_Fsize2, isRead, isWrite); }
        void           efsUndefer4          (const char* why, int lineNumber, int ishift0, int ishift1, int ishift2, int ishift3,
                                                         int exp_current_Fsize2, bool isRead, bool isWrite) {
                                                                                            efsUndeferWkr4(why, lineNumber, ishift0, ishift1, ishift2, ishift3, exp_current_Fsize2, isRead, isWrite); }

        double const * efsRealConst         (int lineNumber, int ishift, int exp_current_Fsize2) { return realPtr + offsetAndChangeState(lineNumber, ishift, exp_current_Fsize2, true , false); }
        double const * efsImagConst         (int lineNumber, int ishift, int exp_current_Fsize2) { return imagPtr + offsetAndChangeState(lineNumber, ishift, exp_current_Fsize2, true , false); }

        double *       realWriteBeforeRead  (int ishift, int exp_current_Fsize2) { return realPtr + offsetAndChangeState(0, ishift, exp_current_Fsize2, false, true ); }
        double *       imagWriteBeforeRead  (int ishift, int exp_current_Fsize2) { return imagPtr + offsetAndChangeState(0, ishift, exp_current_Fsize2, false, true ); }
        double *       realReadWrite        (int ishift, int exp_current_Fsize2) { return realPtr + offsetAndChangeState(0, ishift, exp_current_Fsize2, true , true ); }
        double *       imagReadWrite        (int ishift, int exp_current_Fsize2) { return imagPtr + offsetAndChangeState(0, ishift, exp_current_Fsize2, true , true ); }

        void malloc(int nr_images, int nr_trans, int current_Fsize2) {
            if (debugging) std::cerr << "Exp_Fimgs_shifted::malloc(nr_images:"<<nr_images<<",nr_trans:"<<nr_trans<<",current_Fsize2:"<<current_Fsize2<<")" << std::endl;
            // if too small, free
            if (nr_images > max_nr_images || nr_images*nr_trans > max_ishift) {
                free();
            }
            // Make big enough
            if (max_nr_images == 0) {
                max_nr_images      = nr_images;
                max_ishift         = nr_images*nr_trans;
                max_current_Fsize2 = current_Fsize2;
                state = vNew(char,max_ishift);
                realPtr                 = Heap::allocDoubles(max_ishift   *max_current_Fsize2, __FILE__, __LINE__);
                imagPtr                 = Heap::allocDoubles(max_ishift   *max_current_Fsize2, __FILE__, __LINE__);
                captured_Mctf_invsigma2 = Heap::allocDoubles(max_nr_images*max_current_Fsize2, __FILE__, __LINE__);
                deferred = vNew(Deferred,max_ishift);
            }

            for (int i = 0; i < max_ishift; i++) { setState(i, State_uninit, true); deferred[i].forget(); }
        }

        void zero(int nr_images, int nr_trans, int current_Fsize2) {
            if (debugging) std::cerr << "Exp_Fimgs_shifted::zero()" << std::endl;
            int ishift = nr_images*nr_trans;
            assert(nr_images      <= max_nr_images);
            assert(ishift         <= max_ishift);
            assert(current_Fsize2 <= max_current_Fsize2);
            countAndResetConsumed();
            memset(realPtr, 0, sizeof(double)*ishift*current_Fsize2);
            memset(imagPtr, 0, sizeof(double)*ishift*current_Fsize2);
            for (int i = 0; i < max_ishift; i++) setState(i, State_zero);
        }

        void doOrDeferShiftImageInFourierTransform(
			int ishift, 
			int expCfg_current_size, int exp_current_Fsize2, 
			SOAComplexReadonly& Fimg_aux, double shiftx, double shifty, bool canDefer) {

            // Once the inputs will survive to the deferral point
            // use code similar to the commented out code to do the shifts of only the ones that are needed!
            //
            if (!canDefer) {
                static DoSomePerIter doSomePerIter;
                doSomePerIter.note(debuggingDeferred, iter, [&](int count){
                    std::cerr << "Exp_Fimgs_shifted::Buffer::doOrDeferShiftImageInFourierTransform !canDefer count:" << count << std::endl;
                });
                doShiftImageInFourierTransform(ishift, expCfg_current_size, exp_current_Fsize2, Fimg_aux, shiftx, shifty, -1);
                return;
            }

            auto & d = deferred[ishift];
            omp_set_lock(&d.lock);

            // Hits     if (state[ishift] >= State_deferred) {
            // Hits         std::cerr << "Exp_Fimgs_shifted state[ishift] >= State_deferred" << std::endl;
            // Hits         EXIT_ABNORMALLY;
            // Hits     }

            setState(ishift,State_deferred);
            d.init(expCfg_current_size, exp_current_Fsize2, Fimg_aux, shiftx, shifty);

            if (debugging) {
                static int count;
                static int limit = 1;
                static int limitStep = 25000;
                #pragma omp critical
                if (count++ == limit) {
                    if (limit < 10) limit++; else if (limit < limitStep) limit *= 2; else limit += limitStep;
                    std::cerr << "Exp_Fimgs_shifted defering #" << count << " ishift:" << ishift << std::endl;
                }
            }

            omp_unset_lock(&d.lock);
        }

        void capture_Mctf_invsigma2(int nr_images, int current_Fsize2, std::vector<FimgsData> const & Mctf_invsigma2) {
            assert(nr_images      <= max_nr_images);
            assert(current_Fsize2 <= max_current_Fsize2);
			for (int iimage = 0; iimage < nr_images; iimage++) {
				auto p = Mctf_invsigma2[iimage].rptrAll();
				assert(Mctf_invsigma2[iimage].size() >= current_Fsize2);
                ::copy(captured_Mctf_invsigma2+iimage*current_Fsize2, current_Fsize2, p, current_Fsize2);
			}
        }

        void doOrDeferMultipleByMctf_invsigma2(int ishift, int exp_current_Fsize2, int iimage) {
            // This undefers it to modify it, whether or not the result is used.
            // This should be included in the deferal process
            //
            if (state[ishift] == State_deferred) {
                auto & d = deferred[ishift];
                omp_set_lock(&d.lock);

                bool isDeferred = false;
                if (state[ishift] == State_deferred) {  // make sure not undeferred before got the lock
                    d.captured_iimage_within_Mctf_invsigma2 = iimage;
                    isDeferred = true;
                }

                omp_unset_lock(&d.lock);

                if (isDeferred) return;
            }

            doMultipleByMctf_invsigma2(ishift, exp_current_Fsize2, iimage);
        }

        void free() {
            if (debugging) std::cerr << "Exp_Fimgs_shifted::free()" << std::endl;
            if (debugging) countAndResetConsumed();
            vDelete(deferred);
            Heap::freeDoubles(captured_Mctf_invsigma2);
            Heap::freeDoubles(realPtr);
            Heap::freeDoubles(imagPtr);
            vDelete(state);
            max_nr_images = max_ishift = max_current_Fsize2 = 0;
        }

        double* &      refRealPtr_DONT_USE_UNLESS_ABSOLUTELY_NECESSARY() { return realPtr; }
        double* &      refImagPtr_DONT_USE_UNLESS_ABSOLUTELY_NECESSARY() { return imagPtr; }

    private:
        int max_nr_images;
        int max_ishift;
        int max_current_Fsize2;

        int unusedDeferred, unusedWritten, usedWritten;

        void __declspec(noinline) doShiftImageInFourierTransform(int ishift,
			int expCfg_current_size,
            int exp_current_Fsize2, SOAComplexReadonly& Fimg_aux, 
			double shiftx, double shifty, int captured_iimage_within_Mctf_invsigma2) {
            SOAComplexDouble exp_Fimgs_shifted_aux;
            exp_Fimgs_shifted_aux.real = realPtr + offset(ishift,exp_current_Fsize2);
            exp_Fimgs_shifted_aux.imag = imagPtr + offset(ishift,exp_current_Fsize2);

            shiftImageInFourierTransform(Fimg_aux, exp_Fimgs_shifted_aux, expCfg_current_size, shiftx, shifty, ori_size);

            assert(exp_Fimgs_shifted_aux.real[0] * exp_Fimgs_shifted_aux.real[0] >= 0.0);   // not true for negative indefinites
            assert(exp_Fimgs_shifted_aux.real[0] * exp_Fimgs_shifted_aux.real[0] >= 0.0);   // not true for negative indefinites

            if (captured_iimage_within_Mctf_invsigma2 >= 0) doMultipleByMctf_invsigma2(ishift, exp_current_Fsize2, captured_iimage_within_Mctf_invsigma2);

            setState(ishift, State_written);    // Change from State_deferred after written, so other readers will get the correct values
        }

        void __declspec(noinline) doShiftImageInFourierTransform4(int ishift0, int ishift1, int ishift2, int ishift3,
            SOAComplexReadonly& Fimg_aux0,
            SOAComplexReadonly& Fimg_aux1,
            SOAComplexReadonly& Fimg_aux2,
            SOAComplexReadonly& Fimg_aux3,
            int expCfg_current_size, int exp_current_Fsize2, double shiftx, double shifty, int captured_iimage_within_Mctf_invsigma2) {

            ShiftImageInFourierTransform<double,double> shiftImageInFourierTransform(
                expCfg_current_size, shiftx, shifty, ori_size);

            shiftImageInFourierTransform.transform4(
                Fimg_aux0.real, Fimg_aux0.imag, realPtr + offset(ishift0,exp_current_Fsize2), imagPtr + offset(ishift0,exp_current_Fsize2),
                Fimg_aux1.real, Fimg_aux1.imag, realPtr + offset(ishift1,exp_current_Fsize2), imagPtr + offset(ishift1,exp_current_Fsize2),
                Fimg_aux2.real, Fimg_aux2.imag, realPtr + offset(ishift2,exp_current_Fsize2), imagPtr + offset(ishift2,exp_current_Fsize2),
                Fimg_aux3.real, Fimg_aux3.imag, realPtr + offset(ishift3,exp_current_Fsize2), imagPtr + offset(ishift3,exp_current_Fsize2));

            if (captured_iimage_within_Mctf_invsigma2 >= 0) doMultipleByMctf_invsigma2_4(ishift0, ishift1, ishift2, ishift3, exp_current_Fsize2, captured_iimage_within_Mctf_invsigma2);

            setState(ishift0, State_written);   // Change from State_deferred after written, so other readers will get the correct values
            setState(ishift1, State_written);
            setState(ishift2, State_written);
            setState(ishift3, State_written);
        }

        void doMultipleByMctf_invsigma2(int ishift, int exp_current_Fsize2, int iimage) {

            auto real = realPtr + offset(ishift,exp_current_Fsize2);
            auto imag = imagPtr + offset(ishift,exp_current_Fsize2);

            assert(0 <= iimage && iimage < max_nr_images);
            auto ctf_invsigma2 = captured_Mctf_invsigma2 + iimage*exp_current_Fsize2;

            #pragma ivdep
            for (int n = 0; n < exp_current_Fsize2; n++) {
                real[n] *= ctf_invsigma2[n];
                imag[n] *= ctf_invsigma2[n];
            }
        }

        void doMultipleByMctf_invsigma2_4(int ishift0, int ishift1, int ishift2, int ishift3, int exp_current_Fsize2, int iimage) {

            auto real0 = realPtr + offset(ishift0,exp_current_Fsize2);
            auto real1 = realPtr + offset(ishift1,exp_current_Fsize2);
            auto real2 = realPtr + offset(ishift2,exp_current_Fsize2);
            auto real3 = realPtr + offset(ishift3,exp_current_Fsize2);
            auto imag0 = imagPtr + offset(ishift0,exp_current_Fsize2);
            auto imag1 = imagPtr + offset(ishift1,exp_current_Fsize2);
            auto imag2 = imagPtr + offset(ishift2,exp_current_Fsize2);
            auto imag3 = imagPtr + offset(ishift3,exp_current_Fsize2);

            assert(0 <= iimage && iimage < max_nr_images);
            auto ctf_invsigma2 = captured_Mctf_invsigma2 + iimage*exp_current_Fsize2;

            #pragma ivdep
            for (int n = 0; n < exp_current_Fsize2; n++) {
                real0[n] *= ctf_invsigma2[n];
                real1[n] *= ctf_invsigma2[n];
                real2[n] *= ctf_invsigma2[n];
                real3[n] *= ctf_invsigma2[n];
                imag0[n] *= ctf_invsigma2[n];
                imag1[n] *= ctf_invsigma2[n];
                imag2[n] *= ctf_invsigma2[n];
                imag3[n] *= ctf_invsigma2[n];
            }
        }

        bool efsDeferred(int ishift) {
            assert(ishift < max_ishift);
            return state[ishift] == State_deferred;
        }

        int offsetAndChangeState(int lineNumber, int ishift, int exp_current_Fsize2, bool isRead, bool isWrite) {
            assert(ishift             < max_ishift);
            assert(exp_current_Fsize2 <= max_current_Fsize2);

            efsUndeferWkr("offsetAndChangeState", lineNumber, ishift, exp_current_Fsize2, isRead, isWrite);

            if (debugging && !isRead)
            #pragma omp critical
            {
                // preparing to write one - keep stats on whether the previous use was done and used
                //
                if (state[ishift] == State_read) {
                     usedWritten++;
                } else if (state[ishift] == State_deferred) {
                     unusedDeferred++;
                } else if (state[ishift] == State_written) {
                    unusedWritten++;
                    static int limit = 1;
                    if (unusedWritten >= limit) {
                        std::cerr << "unusedDeferred:" << unusedDeferred << " unusedWritten:" << unusedWritten << " usedWritten:" << usedWritten << std::endl;
                        limit = limit < 10000 ? limit*2 : limit+10000;
                    }
                }
            }

            setState(ishift, isWrite ? State_written : State_read); // This includes the read then write as State_written because the result hasn't been read
            return offset(ishift,exp_current_Fsize2);
        }

        void noteUndeferred(const char* why, int lineNumber, int more, int ishift) {
            static DoSomePerIter doSomePerIter;
            doSomePerIter.note(debuggingDeferred, iter, [&](int count){
                std::cerr << "Exp_Fimgs_shifted undefering #" << count << " ishift:" << ishift;
                if (more > 1) std::cerr << " in a group of " << more;
                std::cerr << " why:" << why << "@" << lineNumber << std::endl;
            });
        }

        void __declspec(noinline) offsetAndChangeStateUndefer(const char* why, int lineNumber, int ishift, int exp_current_Fsize2, bool isRead, bool isWrite) {
            // rare and expensive so move out of the offsetAndChangeState to clean that one's path
            auto & d = deferred[ishift];
            omp_set_lock(&d.lock);

            if (state[ishift] == State_deferred) {  // make sure not undeferred before got the lock

                if (!isRead) {
                    setState(ishift,State_uninit);
                } else {
                    doShiftImageInFourierTransform(ishift, d.expCfg_current_size, d.exp_current_Fsize2, d.Fimg_aux, d.shiftx, d.shifty, d.captured_iimage_within_Mctf_invsigma2);
                    assert(state[ishift] == State_written || state[ishift] == State_read);		// Once it is changed to _written, it can be changed to _read without locking
                    noteUndeferred(why, lineNumber, 1, ishift);
                }

                d.forget();
            }

            omp_unset_lock(&d.lock);
        }

        void efsUndeferWkr(const char* why, int lineNumber, int ishift, int exp_current_Fsize2, bool isRead, bool isWrite) {
            if (state[ishift] == State_deferred) {
                offsetAndChangeStateUndefer(why, lineNumber, ishift, exp_current_Fsize2, isRead, isWrite);
            }
        }

        void efsUndeferWkr4(const char* why, int lineNumber, int ishift0, int ishift1, int ishift2, int ishift3, int exp_current_Fsize2, bool isRead, bool isWrite) {
            // caller saw deferred, but might no longer be
            bool done(false);
            const char* whyNot4 = NULL;
            auto checkIfCanDo = [&](bool exp, const char* expstr) {
                if (!exp) whyNot4 = expstr;
                return whyNot4 == NULL;
            };
#define WHYNOT(EXP) checkIfCanDo(EXP,#EXP)
            if (WHYNOT(isRead)) {
                auto & d0 = deferred[ishift0];
                auto & d1 = deferred[ishift1];
                auto & d2 = deferred[ishift2];
                auto & d3 = deferred[ishift3];
                if (true                                     && WHYNOT(omp_test_lock(&d0.lock))) {
                if (WHYNOT(state[ishift0] == State_deferred) && WHYNOT(omp_test_lock(&d1.lock))) {
                if (WHYNOT(state[ishift1] == State_deferred) && WHYNOT(omp_test_lock(&d2.lock))) {
                if (WHYNOT(state[ishift2] == State_deferred) && WHYNOT(omp_test_lock(&d3.lock))) {
                if (WHYNOT(state[ishift3] == State_deferred)) {
#define SAME(X) WHYNOT(d0.X == d1.X) && WHYNOT(d0.X == d2.X) && WHYNOT(d0.X == d3.X)
                    if (SAME(expCfg_current_size) && SAME(exp_current_Fsize2) && SAME(shiftx) && SAME(shifty) && SAME(captured_iimage_within_Mctf_invsigma2)) {
#undef SAME
                        done = true;
                        // Do all four at once
                        static DoSomePerIter doSomePerIter;
                        doSomePerIter.note(debuggingDeferred, iter, [&](int count){
                            std::cerr << "Exp_Fimgs_shifted efsUndeferWkr4 #" << count << " ishift0:" << ishift0 << " why:" << why;
                            std::cerr << "  *********************************** GOOD" << std::endl;
                        });
                        doShiftImageInFourierTransform4(
                            ishift0,
                            ishift1,
                            ishift2,
                            ishift3,
                            d0.Fimg_aux,
                            d1.Fimg_aux,
                            d2.Fimg_aux,
                            d3.Fimg_aux,
							d0.expCfg_current_size, 
                            d0.exp_current_Fsize2, d0.shiftx, d0.shifty, d0.captured_iimage_within_Mctf_invsigma2);
                        assert(state[ishift0] == State_written);
                        assert(state[ishift1] == State_written);
                        assert(state[ishift2] == State_written);
                        assert(state[ishift3] == State_written);
                        noteUndeferred(why, lineNumber, 4, ishift0);
                    }
                }
                omp_unset_lock(&d3.lock);}
                omp_unset_lock(&d2.lock);}
                omp_unset_lock(&d1.lock);}
                omp_unset_lock(&d0.lock);}
            }
            if (done) return;
#undef WHYNOT
            // Do one at a time
            static DoSomePerIter doSomePerIter;
            doSomePerIter.note(debuggingDeferred, iter, [&](int count){
                std::cerr << "Exp_Fimgs_shifted efsUndeferWkr4 doing singly #" << count << " ishift0:" << ishift0 << " because !(" << whyNot4 << ")" << " why:" << why;
                std::cerr << "  *********************************** BAD" << std::endl;
            });
            offsetAndChangeStateUndefer(why, lineNumber, ishift0, exp_current_Fsize2, isRead, isWrite);
            offsetAndChangeStateUndefer(why, lineNumber, ishift1, exp_current_Fsize2, isRead, isWrite);
            offsetAndChangeStateUndefer(why, lineNumber, ishift2, exp_current_Fsize2, isRead, isWrite);
            offsetAndChangeStateUndefer(why, lineNumber, ishift3, exp_current_Fsize2, isRead, isWrite);
        }

        int offset(int ishift, int exp_current_Fsize2) {
            assert(ishift             < max_ishift);
            assert(exp_current_Fsize2 <= max_current_Fsize2);
            return ishift*exp_current_Fsize2;
        }

        void countAndResetConsumed() {
            if (max_ishift == 0) return;
            if (!debugging) {
                for (int i = 0; i < max_ishift; i++) {
                    setState(i, State_uninit);
                }
            } else {
                auto s0 = state[0];
                int counts[State__length];
                int bad(0);
                for (int i = 0; i < State__length; i++) counts[i] = 0;
                for (int i = 0; i < max_ishift; i++) {
                    auto s = state[i];
                    if (s < 0 || s >= State__length) bad++; else counts[s]++;
                    setState(i, State_uninit);
                }
                std::cerr << "Exp_Fimgs_shifted countAndResetConsumed state=" << (void*)state << " state[0]=" << int(s0) << " counts:";
                for (int i = 0; i < State__length; i++) std::cerr << counts[i] << " ";
                if (bad>0) std::cerr << " bad:" << bad;
                std::cerr << std::endl;
            }
        }

        struct Deferred {
            omp_lock_t			lock;
			int					expCfg_current_size;
            int					exp_current_Fsize2;
            SOAComplexReadonly  Fimg_aux;
            double				shiftx;
            double				shifty;
            int					captured_iimage_within_Mctf_invsigma2;
            Deferred() {
                omp_init_lock(&lock);
                forget();
            }
            ~Deferred() {
                omp_destroy_lock(&lock);
            }
            void forget() {
                exp_current_Fsize2 = 0; Fimg_aux.real = Fimg_aux.imag = NULL; shiftx = 0; shifty = 0; captured_iimage_within_Mctf_invsigma2 = -1;
            }
            void init(int expCfg_current_size, int exp_current_Fsize2, SOAComplexReadonly& Fimg_aux, double shiftx, double shifty) {
                this->expCfg_current_size   = expCfg_current_size;
                this->exp_current_Fsize2    = exp_current_Fsize2;
                this->Fimg_aux              = Fimg_aux;
                this->shiftx                = shiftx;
                this->shifty                = shifty;
                this->captured_iimage_within_Mctf_invsigma2 = -1;
            }
        };

        enum State {State_uninit, State_zero, State_deferred, State_written, State_read, State__length};

        char* state;
        void setState(int i, State to, bool jammed = false) {
            if (debugging) {
                if (i < 0 || max_ishift <= i) {
                    std::cerr << "state[" << i << "]" << " out of bounds " << std::endl;
                    EXIT_ABNORMALLY;
                }
                auto was = jammed ? char(0) : state[i];
                // Note: there is a locking 'problem' here - multiple threads could be moving from written to read at the same time
                // We ignore this problem!
                bool wasBad = (was < 0 || State__length <= was);
                bool toBad  = (to  < 0 || State__length <= to );
                if ((jammed && (i == 0)) || (was != to) && (i == 0 || wasBad || toBad)) {
                    std::cerr << "state[" << i << "]@" << (void*)(state+i) << " was:" << int(was) << " to:" << int(to) << (jammed?" jammed":"") << std::endl;
                    if (wasBad) std::cerr << "  Was bad value!" << std::endl;
                    if (toBad) {
                        std::cerr << "  To bad value!" << std::endl;
                        EXIT_ABNORMALLY;
                    }
                }
            }
            state[i] = to;
        }

        double*     realPtr;
        double*     imagPtr;
        double*     captured_Mctf_invsigma2;
        Deferred*   deferred;
    };

    void Exp_Fimgs_shifted::Buffer::flush(const char* why, int lineNumber) {
        if (len==0) return;
        static DoSomePerIter doSomePerIter;
        doSomePerIter.note(debuggingDeferred, iter, [&](int count){
            std::cerr << "Exp_Fimgs_shifted::Buffer::flush count:" << count << " len:" << len << std::endl;
        });
        int i = 0;
        for (; i+4 <= len; i+=4) parent->efsUndefer4(why, lineNumber, ishifts[i+0], ishifts[i+1], ishifts[i+2], ishifts[i+3], exp_current_Fsize2, true, false);
        for (; i+1 <= len; i+=1) parent->efsUndefer (why, lineNumber, ishifts[i],                                             exp_current_Fsize2, true, false);
        len = 0;
    }

    class ThreadData {
    public:
        static const int maxTempImages = 12;    // Should be a multiple of 4, because of loops that are unrolled
    private:
        class FrefctfPortions : public std::vector<double*> {
        public:
            FrefctfPortions() : alloc_current_Fsize2(0) {}
            void alloc(ExpectationConfig& expCfg, int current_Fsize2) {

				if (alloc_current_Fsize2 < current_Fsize2) {
					fini();
					alloc_current_Fsize2 = current_Fsize2;
				}

				auto currentSize = size();
				auto newSize     = square(expCfg.nr_rot());

				if (currentSize == newSize) return;

				for (size_t i = newSize; i < currentSize; i++) {
					Heap::freeDoubles((*this)[i]);
				}

                resize(newSize);
				if (0) std::cerr << "FrefctfPortions resized to :" << newSize << " at " << __FILE__ << ":" << __LINE__ << std::endl;

                for (int i = currentSize; i < size(); i++) {
                    (*this)[i] = Heap::allocDoubles(current_Fsize2, __FILE__, __LINE__);
                }
            }
            ~FrefctfPortions() { 
				fini(); 
			}
            void fini() {
                for (int i = 0; i < size(); i++) {
					Heap::freeDoubles((*this)[i]);
				}
                resize(0);
            }
		private:
			size_t alloc_current_Fsize2;
        };
        FrefctfPortions _Frefctf_real;
        FrefctfPortions _Frefctf_imag;

		size_t     _Fimg_size;
        SOAComplexDouble _Fimg_nomask;
        SOAComplexDouble _Fimg[maxTempImages];

        void setSOANull() {
			_Fimg_size = 0;
            _Fimg_nomask.real = _Fimg_nomask.imag = NULL;
            for (int i = 0; i < maxTempImages; i++) _Fimg[i].real = _Fimg[i].imag = NULL;
        }
		void freeSOA() {
            for (int i = 0; i < maxTempImages; i++)  {
                Heap::freeDoubles(_Fimg[i].real);
                Heap::freeDoubles(_Fimg[i].imag);
            }
            Heap::freeDoubles(_Fimg_nomask.real);
            Heap::freeDoubles(_Fimg_nomask.imag);
			_Fimg_size = 0;
		}

        int _current_Fsize2, _nr_images;
        VectorOfDouble _max_weight, _wsum_norm_correction, _wsum_sigma2_noise, _wsum_pdf_direction, _wsum_pdf_class, _wsum_prior_offsetx_class, _wsum_prior_offsety_class;

    public:
        ThreadData() : _current_Fsize2(0), _nr_images(0) {
            setSOANull();
        }
        ~ThreadData() {
            fini();
        }
        void fini()  {
            freeSOA();
            _Frefctf_real.fini();
            _Frefctf_imag.fini();
            _max_weight.fini();
            _wsum_norm_correction.fini(); _wsum_sigma2_noise.fini();
            _wsum_pdf_direction.fini(); _wsum_pdf_class.fini();
            _wsum_prior_offsetx_class.fini(); _wsum_prior_offsety_class.fini();
            _nr_images = 0; _current_Fsize2 = 0;
        }

        void alloc(ExpectationConfig& expCfg, int current_Fsize2, int nr_images) {
            _current_Fsize2 = current_Fsize2;
            _nr_images      = nr_images;

            _Frefctf_real.alloc(expCfg, current_Fsize2);
            _Frefctf_imag.alloc(expCfg, current_Fsize2);

            //  thread variable for GetImagesFourierTransformsAndCtfs() and PrecalculateShiftedImagesCtfsAndInvSigma2s()
            //  alloc ori_size x ori_size because in GetImagesFourierTransformsAndCtfs() need most ori_size x ori_size for image
			if (_Fimg_size < square(ori_size)) {
				freeSOA();
				_Fimg_size = square(ori_size);
				_Fimg_nomask.real = Heap::allocDoubles(_Fimg_size, __FILE__, __LINE__);
				_Fimg_nomask.imag = Heap::allocDoubles(_Fimg_size, __FILE__, __LINE__);
				for (int i = 0; i < maxTempImages; i++) {
				    _Fimg[i].real = Heap::allocDoubles(_Fimg_size, __FILE__, __LINE__);
				    _Fimg[i].imag = Heap::allocDoubles(_Fimg_size, __FILE__, __LINE__);
				}
			}

            _wsum_norm_correction     .init(nr_images);
            _wsum_sigma2_noise        .init(current_Fsize2);
            _wsum_pdf_direction       .init(nr_classes);
            _wsum_pdf_class           .init(nr_classes);
            _wsum_prior_offsetx_class .init(nr_classes);
            _wsum_prior_offsety_class .init(nr_classes);
            _max_weight               .init(nr_images);
        }

        double* make_zeroed_wsum_norm_correction() { return Heap::allocZeroedDoubles(_nr_images     , __FILE__, __LINE__); }
        double* make_zeroed_wsum_sigma2_noise()    { return Heap::allocZeroedDoubles(_current_Fsize2, __FILE__, __LINE__); }

        double* Frefctf_real            (int rot_over) { return _Frefctf_real[rot_over]; }  // exp_current_Fsize2;
        double* Frefctf_imag            (int rot_over) { return _Frefctf_imag[rot_over]; }  // exp_current_Fsize2;

        double* Fimg_nomask_real        () { return _Fimg_nomask.real;               }  // TBD
        double* Fimg_nomask_imag        () { return _Fimg_nomask.imag;               }  // TBD
        double* Fimg_real               (int index = 0) { assert(0 <= index && index < maxTempImages); return _Fimg[index].real; }  // ori_size*ori_size
        double* Fimg_imag               (int index = 0) { assert(0 <= index && index < maxTempImages); return _Fimg[index].imag; }  // ori_size*ori_size

        double* wsum_sigma2_noise       () { return _wsum_sigma2_noise       .wptrAll(); }  // exp_current_Fsize2
        double* wsum_norm_correction    () { return _wsum_norm_correction    .wptrAll(); }  // expImages.nr_images()
        double* wsum_pdf_class          () { return _wsum_pdf_class          .wptrAll(); }  // nr_classes
        double* wsum_pdf_direction      () { return _wsum_pdf_direction      .wptrAll(); }  // nr_classes
        double* wsum_prior_offsetx_class() { return _wsum_prior_offsetx_class.wptrAll(); }  // nr_classes
        double* wsum_prior_offsety_class() { return _wsum_prior_offsety_class.wptrAll(); }  // nr_classes

        bool written_wsum_norm_correction_sigma_noise() {
            return _wsum_norm_correction.written();
        }

        void fini_wsum_norm_correction_sigma_noise() {
            _wsum_norm_correction.fini();
            _wsum_sigma2_noise.fini();
        }

        bool written_wsum_pdf_prior_offset_max_weight() {
            return _wsum_pdf_direction.written();
        }

        void fini_wsum_pdf_prior_offset_max_weight() {
            _wsum_pdf_direction      .fini();
            _wsum_pdf_class          .fini();
            _wsum_prior_offsetx_class.fini();
            _wsum_prior_offsety_class.fini();
            _max_weight              .fini();
        }

        void zero_wsum_norm_correction_sigma_noise() {
           _wsum_norm_correction.zero();
           _wsum_sigma2_noise   .zero();
        }

        void zero_wsum_pdf_prior_offset_max_weight() {
            _wsum_pdf_direction      .zero();
            _wsum_pdf_class          .zero();
            _wsum_prior_offsetx_class.zero();
            _wsum_prior_offsety_class.zero();
            _max_weight              .zero();
        }
    };
    const int ThreadData::maxTempImages;

	class ThreadDataForEachThread {
    private:
        size_t      _size;
        ThreadData* _ptr;
    public:
        ThreadDataForEachThread() : _size(0), _ptr(NULL) {
			_ptr = NULL;
		}
        ~ThreadDataForEachThread() { fini(); }
        void init() {
            assert(!_ptr);
            assert(maxthreads>0);
            _size = maxthreads;
            _ptr = vNew(ThreadData,_size);
        }
        void alloc(ExpectationConfig& expCfg, int current_Fsize2, int nr_images) {
            assert(_size > 0);
            for (size_t i = 0; i < _size; i++) _ptr[i].alloc(expCfg, current_Fsize2, nr_images);
        }
        void fini() {
            vDelete(_ptr); _size = 0;
        }
        ThreadData* operator[](int i) { assert(0 <= i && i < _size); return &_ptr[i]; }

        void fini_wsum_norm_correction_sigma_noise() {
            for (size_t i = 0; i < _size; i++) {
                _ptr[i].fini_wsum_norm_correction_sigma_noise();
            }
        }

        void fini_wsum_pdf_prior_offset_max_weight() {
            for (size_t i = 0; i < _size; i++) {
                _ptr[i].fini_wsum_pdf_prior_offset_max_weight();
            }
        }
    };

	struct ExpectationIterationData {
		ExpectationIterationData() {
			static bool called;
			called = true;
		}
#define SEP
#define ELT(T,N) T N;
		MAP2DOPTIMIZER_EXPITER_SCALARS
#undef ELT
#undef SEP
#define ELTVoV(ELEMENT_TYPE, NAME)  VectorOfVector<ELEMENT_TYPE> NAME; bool NAME##_comparable;
#define ELTVEC(ELEMENT_TYPE, NAME)  VectorOfScalar<ELEMENT_TYPE> NAME; bool NAME##_comparable;
#define ELTONE(ELEMENT_TYPE, NAME)  ELEMENT_TYPE                 NAME; bool NAME##_comparable;
#define ELTPTR(ELEMENT_TYPE, NAME)  ELEMENT_TYPE*                NAME; bool NAME##_comparable;
#define SEP													     
		MAP2DOPTIMIZER_EXPITER_ARRAYS SEP
		MAP2DOPTIMIZER_EXPITER_ARRAYS_BOOL
#undef SEP
#undef ELTPTR
#undef ELTONE
#undef ELTVEC
#undef ELTVoV
		Exp_FrefsAndWeight FrefsAndWeight;
		Exp_Fimgs_shifted  Fimgs_shifted;

		void init(
			ExpectationConfig& expCfg,
			ExpectationImages& expImages,
			ThreadDataForEachThread&  threadDataForEachThread,
			int current_size) {

#ifdef DATA_STREAM
            // ------------ initialize model and wsum ---------------- //
            global_data_stream.foutDouble(mlModel.tau2_class[0].wptr(ori_size/2+1), ori_size/2+1, "expectationSetup()_tau2_class", __FILE__, __LINE__);
            global_data_stream.foutDouble(mlModel.tau2_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "expectationSetup()_tau2_class", __FILE__, __LINE__);
            global_data_stream.check();global_data_stream.flush();
#endif
            
            // ------------ initialize model and wsum ---------------- //
            mapModel.setFourierTransformMaps(mlModel.tau2_class, current_size);
            
#ifdef DATA_STREAM
            global_data_stream.foutDouble(mlModel.tau2_class[0].wptr(ori_size/2+1), ori_size/2+1, "expectationSetup()_tau2_class", __FILE__, __LINE__);
            global_data_stream.foutDouble(mlModel.tau2_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "expectationSetup()_tau2_class", __FILE__, __LINE__);
            global_data_stream.check();global_data_stream.flush();
#endif
            
            // Initialise all weighted sums to zero
            mlModel.resetZero();

			MapOptimizer_base_new::resetRandomlyPerturbedSampling(iter);

			//the first time we use small nr_pool to initialize model_Iref,
			//in order to peak much performance,set it larger after iter 1.
			//    if(iter == 2) nr_pool = std::min(300,nr_local_images);
			NODE0ONLY calspace(current_size, nr_pool);

			const int current_Fsize2 = current_size*(current_size / 2 + 1);

			// Number of rotational and translational sampling points
			expCfg.init(sampling2d);

			over_rot_psi.init(); //over_rot_psi_comparable = false;
			over_trans_x.init(); //over_trans_x_comparable = false;
			over_trans_y.init(); //over_trans_y_comparable = false;

			// shifted images and ctfs
			Frefs_Rot_len_per_class = expCfg.nr_psi_adaptive_oversampling()* current_Fsize2;
			FrefsAndWeight.init(nr_classes*expCfg.nr_psi_adaptive_oversampling(), current_Fsize2);

			Rot_significant.init(nr_classes);

			// NOTE :
			//  1)We divide nr_global_images across nodes
			//  if we set nr_pool in first iteration(initialize reference) too large
			//  it will cause a big difference for each reference, so during first iteration just use 8
			//  2)make sure nr_pool is larger than 8,in threadData will malloc nr_pool images data
			const int nr_images = nr_pool;

			Fimgs_shifted.malloc(nr_images, expCfg.nr_trans_adaptive_oversampling(), current_Fsize2);
			Mcoarse_significant.init(nr_images, nr_classes, sampling2d.nrPsi, sampling2d.nrTrans);
			threadDataForEachThread.alloc(expCfg, current_Fsize2, nr_images);


			// NOTE : reset all variable to compare
			zerovec(Rot_significant);
			over_rot_psi.zero(); //over_rot_psi_comparable = true;
			over_trans_x.zero(); //over_trans_x_comparable = true;
			over_trans_y.zero(); //over_trans_y_comparable = true;
			fillvec(highres_Xi2_imgs, 0.0);
			fillvec(min_diff2, 0.0);
			fillvec(old_offsetx, 0.0);
			fillvec(old_offsety, 0.0);
			fillvec(sum_weight, 0.0);

			zerovec(power_imgs);					//power_imgs_comparable = true;
			fillvec(significant_weight, 0.0);
			fillvec(wsum_norm_correction, 0.0);
			zerovec(Fimgs_real);					//Fimgs_real_comparable = true;
			zerovec(Fimgs_imag);					//Fimgs_imag_comparable = true;
			zerovec(Fimgs_nomask_real);				//Fimgs_nomask_real_comparable = true;
			zerovec(Fimgs_nomask_imag);				//Fimgs_nomask_imag_comparable = true;

			Fimgs_shifted.zero(nr_images, expCfg.nr_trans_adaptive_oversampling(), current_Fsize2);
            
            particleModel.setup(nr_images, current_size, coarse_size);
		}

		void __declspec(noinline) resize_image_variables(const int size) {
			imgs.init(size);
			power_imgs.init(size);			power_imgs_comparable = false;
			highres_Xi2_imgs.init(size);
			min_diff2.init(size);
			old_offsetx.init(size);
			old_offsety.init(size);
			wsum_norm_correction.init(size);
			Fimgs_real.init(size);			Fimgs_real_comparable = false;
			Fimgs_imag.init(size);			Fimgs_imag_comparable = false;
			Fimgs_nomask_real.init(size);	Fimgs_nomask_real_comparable = false;
			Fimgs_nomask_imag.init(size);	Fimgs_nomask_imag_comparable = false;
			significant_weight.init(size);
			sum_weight.init(size);
		}
		
		void fini_image_variables() {

			resize_image_variables(0);

			// these are in the init rather than the resize
			//
			over_rot_psi.clear();			over_rot_psi_comparable = false;
			over_trans_x.clear();			over_trans_x_comparable = false;
			over_trans_y.clear();			over_trans_y_comparable = false;

			Fimgs_shifted.free();
			FrefsAndWeight.fini();
            
            particleModel.destroy();
		}
	};

}


// Phase support
//
namespace Map2dOptimizer_new {

	class ChooseImageAndClassStepBase {
#define ChooseImageAndClassStepBase_MEMBERS 		\
		ELT(int , numberOfThreads)				SEP \
        ELT(int , numberOfClassesToDo)			SEP \
        ELT(int , numberOfClassSteps)           SEP \
		ELT(int , iclassStep)                   SEP \
		ELT(int , numberOfImageSteps)           SEP \
        ELT(int , iimageStep)                   SEP \
		ELT(int , numberOfRotSteps)             SEP \
		ELT(int , irotStep)						SEP \
		ELT(int , numberOfIterations)           SEP \
		ELT(bool, goodParallelism)					\
		// end of macro

#define SEP
#define ELT(T,N) T _##N;
		ChooseImageAndClassStepBase_MEMBERS
#undef ELT
#undef SEP

	public:

#define SEP
#define ELT(T,N) T N() const { return _##N; }
		ChooseImageAndClassStepBase_MEMBERS
#undef ELT
#undef SEP

			ChooseImageAndClassStepBase(
			ExpectationConfig const & expCfg,
			ExpectationImages const & expImages,
			int iter,
			int initNumberOfThreads)
			: _numberOfThreads(initNumberOfThreads), _numberOfClassesToDo(expImages.iclass_max() + 1 - expImages.iclass_min())
		{
			// There is a real challenge in getting these numbers right, especially in the present of captureUpto
			// The fewer threads, the fewer captures, so the less work to do and the more chances to combine work within an iteration
			// The more threads, the better the load balancing.
			//
			// A chunk is a _irotStep x _iclassStep x _iimageStep brick of the collapsed loops body
			// 5 chunks per thread is good load balancing
			// On KNL with 60 cores, ring11 is struggling to get 2 chunks per thread
			//
			const int numberOfBodys = _numberOfClassesToDo*expImages.nr_images()*expCfg.nr_rot();

			// Need a _irotStep x _iclassStep x _iimageStep that is close to numberOfBodysPerChunk

			for (int estRotSteps = 1; estRotSteps <= expCfg.nr_rot(); estRotSteps++) {
				// Prefer _irotStep as large as possible per thread, to minimise captures
				_irotStep = divRoundUp(expCfg.nr_rot(), estRotSteps);
				_numberOfRotSteps = divRoundUp(expCfg.nr_rot(), _irotStep);

				for (int estClassSteps = 1; estClassSteps <= _numberOfClassesToDo; estClassSteps++) {
					// Prefer to maximize the number of classes, to minimise captures
					// because this is the shorter dimension
					// x+y is the effort to do x*y amount of work, so best solution has x as large as possible when it is known to be less than y
					_iclassStep = divRoundUp(_numberOfClassesToDo, estClassSteps);
					_numberOfClassSteps = divRoundUp(_numberOfClassesToDo, _iclassStep);

					for (int chunksPerThread = 5; chunksPerThread >= 1; chunksPerThread--) {
						// More chunksPerThread is better because it allows load balancing
						// but it is more important to save on captures than to have better load balancing
						const int numberOfChunks = _numberOfThreads * chunksPerThread;
						const int numberOfBodysPerChunk = divRoundUp(numberOfBodys, numberOfChunks);

						_iimageStep = divRoundUp(numberOfBodysPerChunk, _irotStep*_iclassStep);
						_numberOfImageSteps = divRoundUp(expImages.nr_images(), _iimageStep);

						_numberOfIterations = _numberOfRotSteps * _numberOfClassSteps * _numberOfImageSteps;
						_goodParallelism = (_numberOfIterations > 5 * _numberOfThreads);

						if (_goodParallelism) goto Done;
					}
				}
			}
		Done:
			if (omp_get_thread_num() != 0) return;
			static DoSomePerIter doSomePerIter;
			doSomePerIter.note(false, iter, [&](int count) {
#define SEP <<
#define ELT(T,N) " " << #N << ":" << N()
				std::cerr << "ChooseImageAndClassStepBase " << ChooseImageAndClassStepBase_MEMBERS << std::endl;
#undef ELT
#undef SEP

			});
		}
	};


	class ChooseImageAndClassStep : public ChooseImageAndClassStepBase {
	public:
		ChooseImageAndClassStep(ExpectationConfig const & expCfg, ExpectationImages const & expImages, int iter) : ChooseImageAndClassStepBase(expCfg, expImages, iter, omp_get_max_threads()) {}
	};


	class ChooseImageAndClassStepForKNL : public ChooseImageAndClassStepBase {
	public:
		ChooseImageAndClassStepForKNL(ExpectationConfig const & expCfg, ExpectationImages const & expImages, int iter, int lineNumber) : ChooseImageAndClassStepBase(expCfg, expImages, iter, 60) {

			return;	// BEVIN WHY IS THIS HERE?

			if (goodParallelism() || omp_get_thread_num() != 0) return;

			struct S {
#define SEP
#define	ELT(T,N) T N;
				ChooseImageAndClassStepBase_MEMBERS
#undef ELT
#undef SEP
			};
			static const int doneCapacity = 10;
			static S done[doneCapacity];
			static int head = 0;
			static int filtered = 0;

			for (int i = (head + 1) % doneCapacity; i != head; i = (i + 1) % doneCapacity) {
#define SEP &&
#define	ELT(T,N) (done[i].N == N())
				if (ChooseImageAndClassStepBase_MEMBERS) { filtered++; return; }
#undef ELT
#undef SEP
			}

#define SEP
#define	ELT(T,N) done[head].N = N();
			ChooseImageAndClassStepBase_MEMBERS
#undef ELT
#undef SEP
				head = (head + 1) % doneCapacity;

#define SEP <<
#define ELT(T,N) " " << #N << ":" << N()
			std::cerr << "!chooseImageAndClassStepForKNL.goodParallelism@" << lineNumber << ChooseImageAndClassStepBase_MEMBERS << " filtered:" << filtered << " head:" << head << std::endl;
#undef ELT
#undef SEP
			filtered = 0;

		}
	};

}


// Context
//
namespace Map2dOptimizer_new {

	typedef VectorOfVector<FimgsData> VectorOfFimgsData;

	class VectorOfConstFimgsData {
		VectorOfFimgsData * _ptr;
	public:
		VectorOfConstFimgsData() : _ptr(NULL) {}
		void setPtr(VectorOfFimgsData * to) { _ptr = to; }
		FimgsData const & operator[](int i) { assert(_ptr); return (*_ptr)[i]; }
		bool allOnes() const { return !_ptr; }
	};



}

// Globals
//
namespace Map2dOptimizer_new {

	ExpectationConfig expCfg;
	ExpectationImages expImages;

	ThreadDataForEachThread threadDataForEachThread;
}


#include "./map2d_optimizer_kernel.h"

#include "./initialize.h"

#include <algorithm>
#include <cmath>
#include <map>


namespace Map2dOptimizer_new {

    static void TBD(const char* why) {
        std::cout << "Map2dOptimizer_new::TBD " << why << std::endl;
    }


#define SEP
#define ELT(CLASSNAME, BASECLASS, BASE_COMPONENT, SIZE_COMPUTABLE, LENGTH)  \
    inline bool CLASSNAME::isInitable() { return (SIZE_COMPUTABLE); }   \
    inline int  CLASSNAME::length()     { return (LENGTH); }            // end of macro
    RUNTIME_SIZED_ARRAY_CLASSES
        // Generate these sizing expressions after all the static variables are available
#undef ELT
#undef SEP
	//
	//--


	/***      EM-Iteration     ***/
    
	void expectation();
    
    // expectation nr_pool image each time
    void expectationSomeParticles();
    
    // (de)allocate memory space for each expectation step
    void prepareExpData();
    void endExpData();
    
    /***    some function in expectationsomeparticles functon    ***/

    // get Image(with mask and nomask) and CTF in Fourier Space
    // it will downsize the image from ori_size to current_size
	void getImagesFourierTransformsAndCtfs();

	// get all shifted Fimg,Fimg_nomask and CTF,also get inversed sigma^2
	void getShiftedImagesCtfsAndInvSigma2s();
    void getShiftedImagesNomask();
    
    // get all rotated reference  and the significant rotation
    void getReferenceAllOrientations();

	// get all reference and images 's squared differences
	void getAllSquaredDifferences();
    
    // calculates exp_sum_weight and, for adaptive approach, also exp_significant_weight
    void findAllSignificantPoints();
    
	// convert all squared difference to weight(P = exp(-x))
	void convertSquaredDifferencesToWeights();

    // calculate norm_correction,dLL and Pmax
    void storeWeightedSums();
    
	// update other parameters for refine the model:norm_correction,sigma2_noise,pdf_class,pdf_direction,prior_offsetx(y)_class
	void updateOtherParams();
    
    // add all shifted and rotated images back
    void backProjection();
    
    /***     Maximization step     ***/
    
    void maximization();
    
    // Perform the actual reconstructions
    void maximizationReconstructClass(int iclass);
    
    // Updates all other model parameters (besides the reconstructions)
    void maximizationOtherParameters();
    
    /***    Write files to disk     ***/

    void writeClassesAndMetadata();


#if defined(USEMPI)
    int nodeNameLen;
    char nodeName[MPI_MAX_PROCESSOR_NAME];
#endif

    Images* imagesPtr;

	// Array with pointers to the resolution of each point in a Fourier-space FFTW-like array
    IntPerShell			Npix_per_shell;
    VectorOfScalar<int> Mresol_coarse;
    VectorOfScalar<int> Mresol_fine;

	ExpectationIterationData expIterData;

#define SEP
#define ELT(T,N)    auto & exp_##N = expIterData.N;
#define ELTONE(T,N) ELT(T,N)
#define ELTVEC(T,N) ELT(T,N)
#define ELTVoV(T,N) ELT(T,N)
	MAP2DOPTIMIZER_EXPITER_SCALARS SEP
	MAP2DOPTIMIZER_EXPITER_ARRAYS  SEP
	MAP2DOPTIMIZER_EXPITER_ARRAYS_BOOL
#undef ELTVoV
#undef ELTVEC
#undef ELTONE
#undef ELT
#undef SEP

    VectorOfStruct<MetaDataElem>	exp_metadata;

    VectorOfFimgsData exp_Fctfs_writable;
    VectorOfFimgsData exp_local_Fctfs_writable;

    VectorOfConstFimgsData exp_Fctfs_readonly;			// Always acts like exp_Fctfs_writable even when the values are all 1.0
    VectorOfConstFimgsData exp_local_Fctfs_readonly;	// Sometimes equals exp_Fctfs_readonly
														// Sometimes equals exp_local_Fctfs_writable
														// Sometimes allOnes is true, in which case the others crash

    VectorOfFimgsData			exp_local_Minvsigma2s;

    inline bool isSignificantAnyParticleAnyTranslation(int iclass, int irot);

    static const double sigma2_fudge = 1;


void destroyMap2dOptimizer(){
    //
}

void readImages() {
    readImages(imagesPtr, nr_global_images, first_local_image, last_local_image);
}
    
void prepare() {
    
    if (MapOptimizer_base_new::comparing_versions) {
        prepareAlgorithm(NULL)->waypoint("prepare");
    }
    
    // --------------- initialize Particle Model ----------------- //
    if (particle_diameter < 0.)
        particle_diameter = (ori_size - width_mask_edge) * pixel_size;
    
    particleModel.initialize(ori_size, pixel_size, particle_diameter, width_mask_edge,
                             sigma2_fudge, random_seed, do_norm_correction, do_zero_mask,
                             false,maxthreads,&global_data_stream);
    
    // --------------- initialize MAP Model ---------------- //
    // For do_average_unaligned, always use initial low_pass filter
    // By default, use 0.07 dig.freq. low-pass filter
    // See S.H.W. Scheres (2010) Meth Enzym.
    ini_high = 1.0 / MapOptimizer_base_new::getResolution(round(0.07 * ori_size));
    
    // ini_high ,user set
    mapModel.initialize(nr_classes,ori_size,particle_diameter,pixel_size,ini_high,maxthreads,
                        width_mask_edge,width_fmask_edge,5,true,do_map,2);
    
    // ------------ resolution pointer  --------------------- //
    // Npix_per_shell = (int*)aMalloc(sizeof(int)*(ori_size/2+1),64);
    // Mresol_fine = (int*)aMalloc(sizeof(int)*ori_size*(ori_size/2+1),64);
    // Mresol_coarse = (int*)aMalloc(sizeof(int)*ori_size*(ori_size/2+1),64);
    // ------------ initialize model --------------- //
    // do_scale_correction = false;
    // NOTE : not do scale correction for 2D
    assert(sampling2d.NrDir()==1);
    mlModel.initialize(ori_size, nr_classes, metadata.numberOfGroups(), sampling2d.NrDir());
    
    // Mask and sum all the images, and count how many summed
    // Calculate initial sigma noise model from power_class spectra of the masked images
    //
    Image Mavg;
    
#ifdef DATA_STREAM
    //TODO : check why different with relion.
    global_data_stream.foutInt(metadata.numberOfGroups(), "numberOfGroups", __FILE__, __LINE__);
    global_data_stream.foutInt(metadata.numberOfMicrographs(), "numberOfMicrographs", __FILE__, __LINE__);
    global_data_stream.foutInt(metadata.numberOfParticles(), "numberOfOriginalParticles", __FILE__, __LINE__);
    global_data_stream.foutDouble(Mavg.wptrAll(), ori_size*ori_size, "ref_1", __FILE__, __LINE__);// zero
    global_data_stream.foutDouble(Mavg.wptrAll(), ori_size*ori_size, "ref_1", __FILE__, __LINE__);// zero
    global_data_stream.foutInt(0, "mymodel.orientational_prior_mode", __FILE__, __LINE__);
    global_data_stream.foutInt(mapModel.ref_dim, "mymodel.ref_dim", __FILE__, __LINE__);
    global_data_stream.foutInt(sampling2d.NrDir(), "sampling.NrDirections()", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_direction[0].wptr(sampling2d.NrDir()), sampling2d.NrDir(), "mymodel.pdf_direction[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_direction[nr_classes-1].wptr(sampling2d.NrDir()), sampling2d.NrDir(), "mymodel.pdf_direction_nr_class", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
    {
#ifdef TODO
        mlModel.calculateSumOfPowerSpectraAndAverageImage(Mavg, *imagesPtr, do_zero_mask, metadata, mapModel);
#endif
        if (MapOptimizer_base_new::comparing_versions) {
            prepareAlgorithm(NULL)->waypoint("prepare calculateSumOfPowerSpectraAndAverageImage done");
        }
    }
    
    // check the data whether normalized
    if (true/*check_norm*/) {
        checkNormalize(*imagesPtr, particle_diameter, pixel_size);
    }
    
#if defined(USEMPI)
    //TODO use better types
    //reduce data
    double* local_temp = Heap::allocDoubles(std::max(mlModel.nr_groups,ori_size*ori_size), __FILE__, __LINE__);
    memcpy(local_temp, Mavg.rptrAll(), sizeof(double)*ori_size*ori_size);
    MPI::COMM_WORLD.Reduce(local_temp,Mavg.wptrAll(),ori_size*ori_size,MPI::DOUBLE,MPI::SUM,0);
    
    for (int igroup = 0; igroup < mlModel.nr_groups; igroup++){
        memcpy(local_temp, mlModel.sigma2_noise[igroup].rptrAll(), sizeof(double)*(ori_size/2+1));
        MPI::COMM_WORLD.Reduce(local_temp,mlModel.sigma2_noise[igroup].mptrAll(),ori_size/2+1,MPI::DOUBLE,MPI::SUM,0);
    }
    
    memcpy(local_temp, mlModel.wsum_sumw_group.rptrAll(), sizeof(double)*mlModel.nr_groups);
    MPI::COMM_WORLD.Reduce(local_temp,mlModel.wsum_sumw_group.wptrAll(),mlModel.nr_groups,MPI::DOUBLE,MPI::SUM,0);
    
    Heap::freeDoubles(local_temp);
#endif
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(Mavg.wptrAll(), ori_size*ori_size, "Mavg", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sigma2_noise[0].wptr(ori_size/2+1), ori_size/2+1, "wsum_model_sigma2_noise", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sumw_group.wptrAll()[0], "wsum_model_sumw_group1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sumw_group.wptrAll()[mlModel.nr_groups-1], "wsum_model_sumw_groupN", __FILE__, __LINE__);
    global_data_stream.foutInt((int)mlModel.nr_particles_group.wptrAll()[0], "mymodel_sumw_nr_particles_group1", __FILE__, __LINE__);
    global_data_stream.foutInt((int)mlModel.nr_particles_group.wptrAll()[mlModel.nr_groups-1], "mymodel_sumw_nr_particles_groupN", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
    NODE0ONLY{
        // First calculate average image
        Mavg /= mlModel.wsum_sumw_group.rptrAll()[0];
        
        // Use it, or the provided images, as the mlModel IRefs
        //
#ifdef TODO
        mapModel.initializeRef(Mavg);
#endif
    }
    
#ifdef TOOD
    // Set model_sigma2_noise and model_Iref from averaged poser spectra and Mavg
    NODE0ONLY mlModel.setSigmaNoiseEstimates(Mavg);
#endif
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mlModel.sigma2_noise[0].wptr(ori_size/2+1), ori_size/2+1, "sigma2_noise1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.sigma2_noise[mlModel.nr_groups-1].wptr(ori_size/2+1), ori_size/2+1, "sigma2_noiseN", __FILE__, __LINE__);
    global_data_stream.foutInt(ori_size/2+1, "sigma2_size", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
    // First low-pass filter the initial references
    mapModel.applyLowPassFilter();
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mapModel.Irefs[0].wptr(), ori_size*ori_size, "ref_1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.Irefs[nr_classes-1].wptr(), ori_size*ori_size, "ref_N", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
    // Initialise the model_data_versus_prior ratio to get the initial current_size right
    NODE0ONLY mlModel.initialiseDataVersusPrior(mapModel,tau2_fudge_factor);
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mlModel.data_vs_prior_class[0].wptr(ori_size/2+1), ori_size/2+1, "data_vs_prior_class_1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.data_vs_prior_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "data_vs_prior_class_N", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
#if defined(USEMPI)
    for (int igroup = 0; igroup < mlModel.nr_groups; igroup++)
        MPI::COMM_WORLD.Bcast(mlModel.sigma2_noise[igroup].mptrAll(),ori_size/2+1,MPI::DOUBLE,0);
    
    for (int iclass = 0; iclass < nr_classes; iclass++){
        MPI::COMM_WORLD.Bcast(mapModel.Irefs             [iclass].wptr()   ,ori_size*ori_size,MPI::DOUBLE,0);
        MPI::COMM_WORLD.Bcast(mlModel.tau2_class         [iclass].mptrAll(),(ori_size/2+1),   MPI::DOUBLE,0);
        MPI::COMM_WORLD.Bcast(mlModel.data_vs_prior_class[iclass].mptrAll(),(ori_size/2+1),   MPI::DOUBLE,0);
    }
    
#endif
    
    if (MapOptimizer_base_new::comparing_versions) {
        prepareAlgorithm(NULL)->waypoint("prepare done");
    }
    
    threadDataForEachThread.init();
}
    

/** ========================== EM-Iteration  ================================= */

void iterate() {

    // Update the current resolution and image sizes, and precalculate resolution pointers
    // In later iterations this is done after maximization but before writing output files,
    // so that current resolution is in the output files of the current iteration.
	//
    bool set_by_ini_high = ini_high > 0. && (iter == 0);
    static int nr_iter_wo_resol_gain = 0;
    mapModel.updateCurrentResolution(mlModel.data_vs_prior_class,set_by_ini_high,nr_iter_wo_resol_gain);
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mapModel.current_resolution, "updateCurrentResolution()_current_resolution", __FILE__, __LINE__);
#endif
    
    NODE0ONLY std::cout<<"current_resolution = "<<current_resolution<<std::endl;

	// Init recovery support
	//
	setupStatusTracer();

    // if restarting from a later iteration...
	//
    if (iter != 0) {
        NODE0ONLY std::cout<<"------------------ from "<<iter<<"  --------------------"<<std::endl;
        std::string status_fn = pathRemoveSuffix(star_fn);
        statusTracer.recoveryStatus(status_fn);
	}

    // After the first iteration the references are always CTF-corrected
    refs_are_ctf_corrected = (iter != 0) && do_ctf_correction;

	// do the iterations
	//
    int remainingToEmitTestOutput = 10;                     // BEVIN
    bool has_already_reached_convergence = false;
    for (set_iter(iter + 1); iter <= nr_iter; set_iter(iter+1)) {
        TUNING_FLUSH

#ifdef DATA_STREAM
        global_data_stream.foutInt(iter, "iterate()_iter", __FILE__, __LINE__);
#endif
        
		if (emit_test_output()) {
		    testos() << "~~~~TEST OUTPUT: iterate iter:" << iter << " current_resolution " << current_resolution << std::endl;
		}

		ml_original_waypoint("loop start");
        if (--remainingToEmitTestOutput <= 0) all_versions_agree_emit_test_output_continues = false;

        const double starttime = dtime();

        Npix_per_shell.init(ori_size/2+1); Npix_per_shell.zero();
        Mresol_fine.init(ori_size*(ori_size/2+1)/*current_Fsize2*/);Mresol_coarse.init(ori_size*(ori_size/2+1)/*coarse_Fsize2*/);
        double angularSampler = sampling2d.getAngularSampling();
        mapModel.updateImageSizeAndResolutionPointers(Npix_per_shell, Mresol_coarse, Mresol_fine,
                                                      coarse_size, current_size,
                                                      adaptive_oversampling, angularSampler,
                                                      mlModel.ave_Pmax);

#ifdef DATA_STREAM
        global_data_stream.foutInt(current_size, "updateImageSizeAndResolutionPointers()_current_size", __FILE__, __LINE__);
        global_data_stream.foutInt(coarse_size, "updateImageSizeAndResolutionPointers()_coarse_size", __FILE__, __LINE__);
        global_data_stream.foutInt(Npix_per_shell.wptrAll(), (ori_size/2+1), "updateImageSizeAndResolutionPointers()_Npix_per_shell", __FILE__, __LINE__);
        global_data_stream.foutInt(Mresol_fine.wptrAll(), current_size*(current_size/2+1), "updateImageSizeAndResolutionPointers()_Mresol_fine", __FILE__, __LINE__);
        global_data_stream.foutInt(Mresol_coarse.wptrAll(), coarse_size*(coarse_size/2+1), "updateImageSizeAndResolutionPointers()_Mresol_coarse", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
        
        NODE0ONLY std::cout<<"current_size = "<<current_size<<" coarse_size = "<<coarse_size<<std::endl;

        //set the exp_metadata and exp_image_data and some other data
        //allocate the maximum data
        prepareExpData();
        TUNING_FLUSH

        expectation();
        TUNING_FLUSH

#if defined(USEMPI)
        // TODO

        MPI::COMM_WORLD.Barrier();

#endif

        maximization();
        TUNING_FLUSH

        // Apply masks to the reference images
        // At the last iteration, do not mask the map for validation purposes
        if (do_solvent && !has_converged)
            mapModel.applySolventFlatten("NULL");

        // Re-calculate the current resolution, do this before writing to get the correct values in the output files
        bool set_by_ini_high = false;
        mapModel.updateCurrentResolution(mlModel.data_vs_prior_class,set_by_ini_high,nr_iter_wo_resol_gain);
        TUNING_FLUSH
        
        NODE0ONLY {
			if (guidance_fn != "NULL") {
                statusTracer.checkStatus(guidance_fn+"_iter"+num2str(iter));
			}
        }
        
        // Write output files
        NODE0ONLY {
			if (write_path != "NULL") {
                std::string iterText = num2str(iter);
                std::string fn_class    = write_path+write_fn+"_iter"+iterText;
                std::string fn_metadata = write_path+write_fn+"_iter"+iterText;
				assert(nr_classes == mapModel.Irefs.size());
				// ListOfImages listOfImages(model_Irefs);
				// MapOptimizer_base_new::writeClassesAndMetadata(fn_class, fn_metadata, metadata, listOfImages, testos());
				std::cerr << __FILE__ << ":" << __LINE__ << " **** writeClassesAndMetadata TRACING_DATA NOT WRITTEN ****" << std::endl;
                std::string status_fn = write_path+write_fn+"_iter"+iterText;
                statusTracer.backupStatus(status_fn);
            }
        }

		//update the metadata,free exp_metadata and exp_image_data
		//
        endExpData();

        const double endtime = dtime();
        NODE0ONLY {
			std::cout<<"**** iteration "<<iter<<" completed in "<<endtime-starttime<<" seconds **** census:" << census() <<std::endl<<std::flush;
			showPopulation();
		}
        TUNING_FLUSH
    }
    ml_original_waypoint("loop exit");
}


void __declspec(noinline) resize_exp_image_variables(const int size) {
	exp_metadata            .init(size);
	expIterData.resize_image_variables(size);
	exp_Fctfs_writable.init(size);
	exp_Fctfs_readonly.setPtr(nullptr);

	exp_local_Fctfs_writable.init(size);
	exp_local_Minvsigma2s   .init(size);
	exp_local_Fctfs_readonly.setPtr(nullptr);
}


void prepareExpData() {
	ml_original_waypoint("prepareExpData before.");
	expIterData.init(expCfg, expImages, threadDataForEachThread, current_size);
    ml_original_waypoint("prepareExpData after");
}


void endExpData()
{
	// Reset things in the opposite order to how they are set
	// Not essential, but more consistent in case some depend on others
	// This frees some data
	//
    exp_local_Fctfs_writable.init(0);
    exp_local_Fctfs_readonly.setPtr(nullptr);
    exp_local_Minvsigma2s.init(0);

    exp_metadata.init(0);

	expIterData.fini_image_variables();

    exp_Fctfs_writable.init(0);
    exp_Fctfs_readonly.setPtr(nullptr);

    exp_Rot_significant.init(0);

    exp_Mcoarse_significant.fini();

	// TODO - check these are done at the right time
	Npix_per_shell.init(0);
    Mresol_coarse .init(0);
    Mresol_fine   .init(0);
}


void expectation()
{
    ml_original_waypoint("expectation before");

    NODE0ONLY std::cout << " Expectation iteration " << iter<< " of " << nr_iter<<std::endl;

	for (int nr_images_done = 0; nr_images_done < nr_local_images;) {

		if (0) NODE0ONLY std::cout << " nr_images_done : " <<nr_images_done<<std::endl;
        //
		const int my_first_image  = first_local_image + nr_images_done;
        int my_last_image;
		
        //first time divided the images by nr_classes piece, make the initialized reference different
		{
	        if (iter == 1) {
                // divied the images that can be equally split by nr_classes
                // only use in first iteration to generate reference seed
                int iclass = divide_equally_which_group(nr_global_images, nr_classes, my_first_image);
                int first,last;
                divide_equally(nr_global_images, nr_classes, iclass, first, last);
                int suitable_pool = std::min(nr_pool, last-my_first_image+1);
                my_last_image = std::min(first_local_image+nr_local_images - 1, my_first_image + suitable_pool - 1);
			}
            else{
                my_last_image = std::min(first_local_image+nr_local_images - 1, my_first_image + nr_pool - 1);
            }
            expImages.init(my_first_image, my_last_image);
            
            int iclass_min = 0;
            int iclass_max = nr_classes - 1;
            if (iter == 1){
                iclass_min = iclass_max = divide_equally_which_group(nr_global_images, nr_classes, my_first_image);
            }
			expImages.setClass(iclass_min, iclass_max);
        }
        
		exp_nr_images = my_last_image - my_first_image + 1;
        assert(expImages.nr_images() == exp_nr_images);
		resize_exp_image_variables(expImages.nr_images());

#ifdef DATA_STREAM
        global_data_stream.foutInt(expImages.first_image(), "expectation()_my_first_ori_particle", __FILE__, __LINE__);
        global_data_stream.foutInt(expImages.last_image(), "expectation()_my_last_ori_particle", __FILE__, __LINE__);
        global_data_stream.foutInt(expImages.nr_images(), "expectation()_nr_pool", __FILE__, __LINE__);
        global_data_stream.foutInt(0, "expectationOneParticle()_do_firstiter_cc", __FILE__, __LINE__);
        global_data_stream.foutInt(expImages.iclass_min(), "expectationOneParticle()_exp_iclass_min", __FILE__, __LINE__);
        global_data_stream.foutInt(expImages.iclass_max(), "expectationOneParticle()_exp_iclass_max", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
        
        // if(iter ==1 && node == 0) std::cout<<"expImages.first_image() : "<<expImages.first_image()<<",expImages.last_image() : "<<expImages.last_image()<<" expImages.iclass_min() : "<<expImages.iclass_min()<<std::endl;
        //prepare exp_image_data and exp_metadata
        //each node keep nr_local_images's images data and all(nr_global_images) metadata
		auto & images = *imagesPtr;
        
        auto const nr_images = expImages.nr_images();
        for (int iimage = 0; iimage < nr_images; iimage++) {
			copy(exp_imgs[iimage], images[iimage+my_first_image-first_local_image]);
			auto & md = metadata[iimage+my_first_image];
            exp_metadata[iimage] = md;
        }

		//get Image(with mask and nomask) and CTF in Fourier Space
        //it will downsize the image from ori_size to current_size
#ifdef TODO
        particleModel.prepare(exp_imgs, exp_metadata, exp_first_image, expImages.nr_images());
#endif
        //
        ml_original_waypoint("getImagesFourierTransformsAndCtfs before");
        exp_Fctfs_readonly.setPtr(nullptr);
        particleModel.setFourierTransforms(exp_power_imgs, exp_highres_Xi2_imgs, exp_old_offsetx, exp_old_offsety, mlModel, do_ctf_correction);
        //
        exp_Fctfs_readonly.setPtr((VectorOfFimgsData*)&particleModel.Fctfs);
        ml_original_waypoint("getImagesFourierTransformsAndCtfs after");
        
        // perform the actual expectation step on several particles
        expectationSomeParticles();

        // Also monitor the changes in the optimal orientations and classes
        // TOOD.... use HiddenVarMonitor
        // monitorHiddenVariableChanges(my_first_image, my_last_image);
		ml_original_waypoint("monitorHiddenVariableChanges");
        //update the metadata
        update_metadata(metadata, my_first_image, exp_metadata.mptrAll(), expImages.nr_images());

        //
        nr_images_done += my_last_image - my_first_image + 1;

    }
    ml_original_waypoint("expectation");
}

void __declspec(noinline) expectationSomeParticles()
{
    TUNING_SCOPE_STEP(expectationSomeParticles)
    ml_original_waypoint("expectationSomeParticles before");
    // Only perform a second pass when using adaptive oversampling
    int nr_sampling_passes = (adaptive_oversampling > 0) ? 2 : 1;

    // Pass twice through the sampling of the entire space of rot, tilt and psi
    // The first pass uses a coarser angular sampling and possibly smaller FFTs than the second pass.
    // Only those sampling points that contribute to the highest x% of the weights in the first pass are oversampled in the second pass
    // Only those sampling points will contribute to the weighted sums in the third loop below
    for (int ipass = 0; ipass < nr_sampling_passes; ipass++) {

		// Anything that depended on the current size can no longer be compared
		//

        // Use smaller images in the first pass, larger ones in the second pass
        // Use coarse sampling in the first pass, oversampled one the second pass
		expCfg.setIpass(sampling2d, ipass, current_size, coarse_size);

        //initialize sampling data
#ifdef TODO
        sampling2d.getAllOrientations(expCfg.current_oversampling(),exp_over_rot_psi.wptrAll());
        sampling2d.getAllTranslations(expCfg.current_oversampling(),exp_over_trans_x.wptrAll(),exp_over_trans_y.wptrAll());
#endif
        {   TUNING_SCOPE_STEP(serial_exp_Mweight_init)
            auto capacityPerClass = sampling2d.NrPoints(adaptive_oversampling);
            exp_Mweight.init(
                exp_Mcoarse_significant.nr_images(),
                exp_Mcoarse_significant.nr_classes(),
                capacityPerClass);
        }

        // NOTE : reset all variable to compare
        // NOTE : make the exp_rot same as new version
        expIterData.FrefsAndWeight.zero(__LINE__);

        //get all shifted Fimg,Fimg_nomask and CTF,also get inversed sigma^2
        getShiftedImagesCtfsAndInvSigma2s();

        //get all rotated reference  and the significant rotation
        getReferenceAllOrientations();

        //get all reference and images 's squared differences
        getAllSquaredDifferences();

        //convert all squared differences to weight,and find significant(maximum) weight for each image
        convertSquaredDifferencesToWeights();

        findAllSignificantPoints();

        // last ipass iteration, update some parameter for maximization
        if (expCfg.ipass() == nr_sampling_passes-1) {

            expCfg.setCurrentSize(current_size);
	            // For the reconstruction step use model current_size, but the other expCfg members are unchanged

            updateOtherParams();

            getShiftedImagesNomask();

            backProjection();

            storeWeightedSums();
        }

        exp_Mweight.fini();

    }// end loop over 2 expCfg.ipass() iterations
	expCfg.setIpass(nr_sampling_passes);

    ml_original_waypoint("expectationSomeParticles after");
}


struct GetShiftedImagesCtfsAndInvSigma2s_new_Debugging{
    int done;
    int saved;
    int next;
    GetShiftedImagesCtfsAndInvSigma2s_new_Debugging() : done(0), saved(0), next(1) {}
} getShiftedImagesCtfsAndInvSigma2s_new_Debugging;


void __declspec(noinline) getShiftedImagesCtfsAndInvSigma2s()
{   TUNING_SCOPE_STEP(getShiftedImagesCtfsAndInvSigma2s)

    ml_original_waypoint("getShiftedImagesCtfsAndInvSigma2s before.");

#ifdef DATA_STREAM
    global_data_stream.foutDouble(0, "##################start_precalculateShiftedImagesCtfsAndInvSigma2s#####################", __FILE__, __LINE__);
#endif
    
    const int current_Fsize2     = current_size*(current_size/2+1);
    const int exp_current_Fsize2 = expCfg.current_size()*(expCfg.current_size()/2+1);

    if (current_size == expCfg.current_size()) {
        exp_local_Fctfs_readonly = exp_Fctfs_readonly;
    } else {
        exp_local_Fctfs_readonly.setPtr(&exp_local_Fctfs_writable);
    }

    // defer when in fine searching step,do not windows Fctfs and FImags
    const bool canDefer = (current_size == expCfg.current_size());

	auto const nr_images = expImages.nr_images();
    TUNING_SCOPE_PARA_BEGIN(getShiftedImagesCtfsAndInvSigma2s,nr_images)
	#pragma omp parallel for  // TODO MORE PARALLELISM
    for (int iimage = 0; iimage < nr_images; iimage++)
    {   TUNING_SCOPE_ITER_BEGIN(getShiftedImagesCtfsAndInvSigma2s)

#ifdef DATA_STREAM
        global_data_stream.foutInt(iimage, "precalculateShiftedImagesCtfsAndInvSigma2s()_start------", __FILE__, __LINE__);
#endif
        
        int tid = omp_get_thread_num();
        auto td = threadDataForEachThread[tid];

        SOAComplexReadonly exp_Fimgs_aux_ro;														// input of the windowFourierTransform
#ifdef TODO
        exp_Fimgs_aux_ro.real = particleModel.Fimages_mask_fine_real[iimage].wptrAll();//exp_Fimgs_real[iimage].rptrAll();
        exp_Fimgs_aux_ro.imag = particleModel.Fimages_mask_fine_imag[iimage].wptrAll();//exp_Fimgs_imag[iimage].rptrAll();
#endif
        auto exp_Fctfs_aux = exp_Fctfs_readonly[iimage].rptrAll();								// input

        SOAComplexReadonly Fimg_aux_ro;															// output of the windowFourierTransform

        // Downsize Fimg and Fctf to exp_current_image_size, also initialise Fref and Fimg_shift to the right size
        // In the second pass of the adaptive approach this will have no effect,
        // since then exp_current_image_size will be the same as the size of exp_Fctfs
        if (false) {
			if (current_size != expCfg.current_size()) {	// this happens
				static int count;
				if (count++ < 10) std::cerr << "current_size != expCfg.current_size() at line " << __LINE__ << std::endl;
			}
        }

        if (false && omp_get_thread_num() == 0)
        #pragma omp critical
        {
            // This measurement shows that we are saving huge numbers of Transforms
            auto & saved = getShiftedImagesCtfsAndInvSigma2s_new_Debugging.saved;
            auto & done  = getShiftedImagesCtfsAndInvSigma2s_new_Debugging.done;
            auto & next  = getShiftedImagesCtfsAndInvSigma2s_new_Debugging.next;
            if (current_size == expCfg.current_size()) saved++; else done++;
            if ((done + saved) == next) { next *= 2; std::cerr << "getShiftedImagesCtfsAndInvSigma2s_new wt done:" << done << " saved:" << saved << " on:" << omp_get_thread_num() << std::endl; }
        }

        if (canDefer) {
            Fimg_aux_ro.real = exp_Fimgs_aux_ro.real;
            Fimg_aux_ro.imag = exp_Fimgs_aux_ro.imag;
        } else {
            static DoSomePerIter doSomePerIter;
            doSomePerIter.note(expIterData.Fimgs_shifted.debuggingDeferred, iter, [&](int count){
                std::cerr << "getShiftedImagesCtfsAndInvSigma2s !canDefer count:" << count
                    << "current_size:" << current_size << " != " << "expCfg.current_size():" << expCfg.current_size()
                    << std::endl;
            });

            // TODO This maybe could also be deferred if it is expensive
			SOAComplexDouble Fimg_aux;
            Fimg_aux_ro.real = Fimg_aux.real = td->Fimg_real();
            Fimg_aux_ro.imag = Fimg_aux.imag = td->Fimg_imag();

            windowFourierTransform(exp_Fimgs_aux_ro, current_size, Fimg_aux, expCfg.current_size());

			double* exp_local_Fctfs_aux = exp_local_Fctfs_writable[iimage].wptrAll();
            windowTransform(exp_Fctfs_aux, current_size, exp_local_Fctfs_aux, expCfg.current_size());
        }

#ifdef DATA_STREAM
        global_data_stream.foutDouble(const_cast<double*>(Fimg_aux_ro.real), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_Fimg_real", __FILE__, __LINE__);
        global_data_stream.foutDouble(const_cast<double*>(Fimg_aux_ro.imag), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_Fimg_imag", __FILE__, __LINE__);
        global_data_stream.foutDouble(const_cast<double*>(exp_local_Fctfs_readonly[iimage].rptrAll()), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_local_Fctfs", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
        
        // Store all translated variants of Fimg	// TODO MORE PARALLELISM MAYBE POSSIBLE HERE IF !canDefer
        int my_trans_image = iimage * expCfg.nr_trans() * expCfg.nr_over_trans();
        for (int itrans = 0; itrans < expCfg.nr_trans(); itrans++)//100(offset_range=10) or 400(offset_range=20),offset_step = 2
        {
            for (int iover_trans = 0; iover_trans < expCfg.nr_over_trans(); iover_trans++)//1(ipass=0) and 4(ipass=1)
            {
                // Shift through phase-shifts in the Fourier transform
                // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
                double shiftx = exp_over_trans_x[itrans*expCfg.nr_over_trans()+iover_trans];
                double shifty = exp_over_trans_y[itrans*expCfg.nr_over_trans()+iover_trans];

                expIterData.Fimgs_shifted.doOrDeferShiftImageInFourierTransform(
					my_trans_image,
					expCfg.current_size(), exp_current_Fsize2, Fimg_aux_ro, shiftx, shifty, canDefer);
                    // Can't really defer until the Fimg_aux is saved

                my_trans_image++;
            }
        }

        auto const         Minvsigma2 = exp_local_Minvsigma2s[iimage].wptrAll();
        int  const * const myMresol   = (expCfg.current_size() == coarse_size) ? Mresol_coarse.rptrAll() : Mresol_fine.rptrAll();

        // With group_id and relevant size of Fimg, calculate inverse of sigma^2 for relevant parts of Mresol
        for (int n = 0; n < exp_current_Fsize2; n++) {
            auto const ires = myMresol[n];
            // Exclude origin (ires==0) from the Probability-calculation
            // This way we are invariant to additive factors
            if (ires > 0) {
                Minvsigma2[n] = CHECK_NOT_NAN_OR_INF(1. / (sigma2_fudge * mlModel.sigma2_noise[iimage][ires]));
			} else {
                Minvsigma2[n] = 0;
			}
        }

		if (emit_test_output_prolific() && iimage == 0) {
			testos() << "getShiftedImagesCtfsAndInvSigma2s iimage:0" 
				<< " sigma2_fudge:" << sigma2_fudge
				<< std::endl;
			for (int n = 0; n < 5; n++) {
				int ires = myMresol[n];
				testos() << "n:" << n 
					<< " Minvsigma2[n:"<<n<<"]:" << Minvsigma2[n]
					<< " ires:" << ires;
				if (ires > 0) testos() << " model_sigma2_noise[ires:"<<ires<<"]:" << mlModel.sigma2_noise[iimage][ires];
				testos() << std::endl;
			}
		}
                    
#ifdef DATA_STREAM
        global_data_stream.foutDouble(Minvsigma2, exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_local_Minvsigma2s", __FILE__, __LINE__);
//        global_data_stream.foutDouble(exp_Fimgs_shifted_real.wptr(iimage, 0, exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_Fimgs_shifted_real1", __FILE__, __LINE__);
//        global_data_stream.foutDouble(exp_Fimgs_shifted_imag.wptr(iimage, 0, exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_Fimgs_shifted_imagn", __FILE__, __LINE__);
//        global_data_stream.foutDouble(exp_Fimgs_shifted_real.wptr(iimage, exp_nr_trans*exp_nr_over_trans-1, exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_Fimgs_shifted_realN", __FILE__, __LINE__);
//        global_data_stream.foutDouble(exp_Fimgs_shifted_imag.wptr(iimage, exp_nr_trans*exp_nr_over_trans-1, exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_Fimgs_shifted_imagN", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif

        TUNING_SCOPE_ITER_END
    }
    TUNING_SCOPE_PARA_END

    ml_original_waypoint("getShiftedImagesCtfsAndInvSigma2s after.");

}

void __declspec(noinline) getShiftedImagesNomask()
{   TUNING_SCOPE_STEP(getShiftedImagesNomask)

    int current_Fsize2 = current_size*(current_size/2+1);
    int exp_current_Fsize2 = expCfg.current_size()*(expCfg.current_size()/2+1);

	auto const nr_images     = expImages.nr_images();
	auto const nr_trans      = expCfg.nr_trans();
	auto const nr_over_trans = expCfg.nr_over_trans();
    TUNING_SCOPE_PARA_BEGIN(getShiftedImagesNomask, nr_images*nr_trans*nr_over_trans)
	#pragma	omp for collapse(3)
    for (int iimage = 0; iimage < nr_images; iimage++)
    {
        for (int itrans = 0; itrans < nr_trans; itrans++)//100(offset_range=10) or 400(offset_range=20),offset_step = 2
        {
            for (int iover_trans = 0; iover_trans < nr_over_trans; iover_trans++)//4(ipass=1)
            {   TUNING_SCOPE_ITER_BEGIN(getShiftedImagesNomask)

                int my_trans_image = (iimage*nr_trans + itrans)*nr_over_trans + iover_trans;

                SOAComplexReadonly exp_Fimgs_nomask_aux;
                exp_Fimgs_nomask_aux.real = exp_Fimgs_nomask_real[iimage].rptrAll();  // The input
                exp_Fimgs_nomask_aux.imag = exp_Fimgs_nomask_imag[iimage].rptrAll();

                SOAComplexReadonly Fimg_nomask_aux_ro;

                // expCfg.current_size() is same as current_size
                bool canDefer = false;
                if (current_size != expCfg.current_size())
                #pragma omp critical
                {
                    static int count;
                    if (count++ < 10) std::cerr << "current_size != expCfg.current_size() at line " << __LINE__ << std::endl;
                    EXIT_ABNORMALLY;
                }

                if (current_size != expCfg.current_size()) {
                    int tid = omp_get_thread_num();                                                                     // use a temp
                    auto td = threadDataForEachThread[tid];
                    auto td_Fimg_real = td->Fimg_real();
                    auto td_Fimg_imag = td->Fimg_imag();

					SOAComplexDouble Fimg_nomask_aux;
                    Fimg_nomask_aux_ro.real = Fimg_nomask_aux.real = td_Fimg_real;
                    Fimg_nomask_aux_ro.imag = Fimg_nomask_aux.imag = td_Fimg_imag;

                    windowFourierTransform(exp_Fimgs_nomask_aux, current_size, Fimg_nomask_aux, expCfg.current_size());      // THIS IS ALWAYS A REDUNDANT COPY, so never used
                } else {
                    canDefer = true;
                    Fimg_nomask_aux_ro.real = exp_Fimgs_nomask_aux.real;                                                   // Just use from the original location
                    Fimg_nomask_aux_ro.imag = exp_Fimgs_nomask_aux.imag;
                }

                // Shift through phase-shifts in the Fourier transform
                // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
                double shiftx = exp_over_trans_x[itrans*nr_over_trans+iover_trans];
                double shifty = exp_over_trans_y[itrans*nr_over_trans+iover_trans];

                expIterData.Fimgs_shifted.doOrDeferShiftImageInFourierTransform(
					my_trans_image, 
					expCfg.current_size(), exp_current_Fsize2, Fimg_nomask_aux_ro, shiftx, shifty, canDefer);
                    // Can't really defer until the Fimg_nomask is saved
                    // which I suspect it already is...

                TUNING_SCOPE_ITER_END
            }
        }
    }
    TUNING_SCOPE_PARA_END

    ml_original_waypoint("getShiftedImagesNomask after.");

}


void getReferenceAllOrientations()
{   TUNING_SCOPE_STEP(getReferenceAllOrientations)

	auto const iclass_min  = expImages.iclass_min(); 
	auto const iclass_max  = expImages.iclass_max();
    auto const nr_rot      = expCfg.nr_rot();
	auto const nr_over_rot = expCfg.nr_over_rot();

    int  const exp_current_Fsize2 = expCfg.current_size()*(expCfg.current_size()/2+1);

    ml_original_waypoint("getReferenceAllOrientations before.");

	TUNING_SCOPE_PARA_BEGIN(getReferenceAllOrientations,(iclass_max+1-iclass_min)*nr_rot)
	#pragma omp parallel for collapse(2)
    for (int iclass = iclass_min; iclass <= iclass_max; iclass++)//50 or 100
    {
        for (int ipsi = 0; ipsi < nr_rot; ipsi++)//36(psi_step=10)
        {

            if (mlModel.pdf_class.rptrAll()[iclass] <= 0.)
                continue;

            int iorientclass_offset = iclass * expCfg.nr_rot();

            int iorientclass = iorientclass_offset + ipsi;
            //In the first pass, always proceed
            //In the second pass, check whether one of the translations for this orientation of any of
            //the particles had a significant weight in the first pass
            //if so, proceed with projecting the reference in that direction
            if (expCfg.ipass() == 0) {
                exp_Rot_significant[iclass][ipsi] = true;
            }
            else
                exp_Rot_significant[iclass][ipsi] = isSignificantAnyParticleAnyTranslation(iclass,ipsi);
        }//end loop of ipsi
    }// end loop of iclass
	TUNING_SCOPE_PARA_END

	TUNING_SCOPE_PARA_BEGIN(getReferenceAllOrientations2,(iclass_max+1-iclass_min)*nr_rot*nr_over_rot)
	#pragma omp parallel for collapse(3)
    for (int iclass = iclass_min; iclass <= iclass_max; iclass++)//50 or 100
    {
        for (int ipsi = 0; ipsi < nr_rot; ipsi++)//36
        {
            for (int iover_rot = 0; iover_rot < nr_over_rot; iover_rot++)//1,2
            {
                double pdf_orientation = mlModel.pdf_direction[iclass][0];

                int iorientclass_offset = iclass * nr_rot;

                int iorientclass = iorientclass_offset + ipsi;

                if (!exp_Rot_significant[iclass].rptrAll()[ipsi] || pdf_orientation <= 0. || mlModel.pdf_class.rptrAll()[iclass] <= 0.)
                    continue;

                //transposition matrix
                double A[3][3];
                // Take tilt-series into account
#ifdef TODO
                Euler_angles2matrix(0, 0, exp_over_rot_psi[ipsi*nr_over_rot+iover_rot], A);
#endif
                SOAComplexDouble exp_Frefs_Rot_aux;
                exp_Frefs_Rot_aux.real = expIterData.FrefsAndWeight.frefs_Rot_real_write1st(__LINE__, iorientclass*nr_over_rot+iover_rot, exp_current_Fsize2);
                exp_Frefs_Rot_aux.imag = expIterData.FrefsAndWeight.frefs_Rot_imag_write1st(__LINE__, iorientclass*nr_over_rot+iover_rot, exp_current_Fsize2);
#ifdef TODO
                mapModel.get2DFourierTransform(iclass, exp_Frefs_Rot_aux.real, exp_Frefs_Rot_aux.imag, exp_current_size, A, false);
#endif
                for (int i = 0; i < expCfg.current_size(); i++) {
                    checkFrefRotValue(exp_Frefs_Rot_aux.real[i]);
                    checkFrefRotValue(exp_Frefs_Rot_aux.imag[i]);
                }
            }// end loop of iover_rot
        }//end loop of ipsi
    }// end loop of iclass
	TUNING_SCOPE_PARA_END

    ml_original_waypoint("getReferenceAllOrientations after.");

}

static void doMweightTestOutput() {
    if (!emit_test_output_prolific()) return;
	for (int iimage = 0; iimage < expImages.nr_images(); iimage++) {
		for (int iclass = expImages.iclass_min(); iclass <= expImages.iclass_max(); iclass++) {
			if (!TestOutputMweight::shouldDo(iclass,iimage)) continue;
			TestOutputMweight testOutputMweight;
	        for (int irot = 0; irot < expCfg.nr_rot(); irot++) {
				for (int iover_rot = 0; iover_rot < expCfg.nr_over_rot(); iover_rot++) {
					int rot_over = irot * expCfg.nr_over_rot() + iover_rot;
	                for (int itrans = 0; itrans < expCfg.nr_trans(); itrans++) {
						for (int iover_trans = 0; iover_trans < expCfg.nr_over_trans(); iover_trans++) {
	                        int rot_trans_over = (rot_over*expCfg.nr_trans() + itrans)*expCfg.nr_over_trans() + iover_trans;
							auto diff2 = exp_Mweight.get(iimage, iclass, rot_trans_over);
							testOutputMweight.note(rot_trans_over, diff2);
						}
					}
				}
			}
			testos() << "~~~~ TEST OUTPUT: ";
			testOutputMweight.print(testos());
		}
	}
}

static void getAllSquaredDifferences()
{
	auto const iclass_min	 = expImages.iclass_min(); 
	auto const iclass_max	 = expImages.iclass_max();
    auto const nr_rot		 = expCfg.nr_rot();
	auto const nr_over_rot	 = expCfg.nr_over_rot();
    auto const nr_images	 = expImages.nr_images(); 
    auto const nr_trans		 = expCfg.nr_trans();
	auto const nr_over_trans = expCfg.nr_over_trans();

    TUNING_SCOPE_STEP(getAllSquaredDifferences)
    ml_original_waypoint("getAllSquaredDifferences before");

    const int exp_current_Fsize2 = expCfg.current_size()*(expCfg.current_size()/2+1);

    // Trick!
    //      The used coarse maybe smaller than exp_Mweight_xsize depend on two thing:
    //      1) in first EM-interation,  not all class used in big while loop
    //      2) some not significant rot,direction,class,trans will be escapsed
    // so we need some minimum integer number to flag this
    //
    exp_Mweight.invalidateAll();
    exp_Mweight.setSizePerClass(nr_rot*nr_over_rot*nr_trans*nr_over_trans);

    static const double undefWeight = -777.777;
    static const double unsetWeight = Exp_Mweight_unsetWeight;

    static const int interesting_iter        = 1 ;
    static const int interesting_ipass       = 0 ;
    static const int interesting_iimage      = 0 ;
    static const int interesting_iclass      = 0 ;
    static const int interesting_irot        = 36;
    static const int interesting_iover_rot   = 0 ;
    static const int interesting_itrans      = 0 ;
    static const int interesting_iover_trans = 0 ;

    static const bool checkForUnsetUsed = false;

#if defined(LOADBALANCEANALYZER)
    LoadBalanceAnalyzer lba("getAllSquaredDifferences", __FILE__, __LINE__, nr_images*(iclass_max+1-iclass_min));
#endif

    if (false) {
        static int count;
        if (count < 2 // This test really cuts down on the number of stalls the #pragma omp critical causes
            && nr_images*(iclass_max+1-iclass_min) < 4*omp_get_max_threads())
        #pragma omp critical
        {
            if (count++ < 2) {
                std::cerr << "getAllSquaredDifferences - not enough iterations for good parallelism" << std::endl
                    #define P(N) << " " << #N << ":" << N
                    P(nr_images) P(iclass_max) P(iclass_min) P(omp_get_max_threads())
                    #undef P
                    << std::endl;
            }
        }
    }


    // The kernel can only process one image at a time, so its loop must have images as the outermost case
    // The transforms want to process multiple images but one transform, so its loop must have the images as the innermost loop

    // FIRST LOOP PREDICTS WHETHER A KERNEL APPEND WILL HAPPEN IN THE SECOND LOOP
    //
    {
        // This is surprisingly easy and fast
        // However it is outside the second loop, because the second loop divides the images into too small a set to get the parallelism we want here
        // Put a range of images with the same itrans and iover_trans into one thread to minimize the number of sharing clashes
        //
        const int iimageStepForKNL = std::max(1,nr_images/60/5);
        const int iimageStep = std::max(1,nr_images/omp_get_num_threads()/5);

		TUNING_SCOPE_PARA_BEGIN(gsd_newer_FIRST,(nr_images+iimageStepForKNL-1)/iimageStepForKNL*nr_trans*nr_over_trans)
        #pragma omp parallel for collapse(3) schedule(dynamic)
        for (int iimageBegin = 0; iimageBegin < nr_images; iimageBegin += iimageStep) {		
            for (int itrans = 0; itrans < nr_trans; itrans++) {
				for (int iover_trans = 0; iover_trans < nr_over_trans; iover_trans++) {

					for (int iclass = iclass_min; iclass < iclass_max + 1; iclass++) {
						if (!(mlModel.pdf_direction[iclass].rptrAll()[0] <= 0. || mlModel.pdf_class.rptrAll()[iclass] <= 0.)) {
							for (int irot = 0; irot < nr_rot; irot++) {

								//int iclass_rot = iclass * nr_rot + irot;
								if (!exp_Rot_significant[iclass][irot]) continue;

                                Exp_Fimgs_shifted::Buffer buffer(&expIterData.Fimgs_shifted, "getAllSquaredDifferences", __LINE__, exp_current_Fsize2);

								const int iimageEnd = std::min(iimageBegin + iimageStep, nr_images);
                                for (int iimage = iimageBegin; iimage < iimageEnd; iimage++) {

									if (expCfg.ipass() != 0 && !exp_Mcoarse_significant.get(iimage, iclass, irot, itrans)) {
                                        continue;
                                    }

                                    int ishift = iimage * nr_trans * nr_over_trans + itrans * nr_over_trans + iover_trans;

                                    if (expIterData.Fimgs_shifted.efsDeferred(ishift, exp_current_Fsize2)) 
										buffer.append("getAllSquaredDifferences", __LINE__, ishift);
                                }
                            }
                        }
                    }
                }
            }
        }
		TUNING_SCOPE_PARA_END
    }

    // SECOND LOOP
	ChooseImageAndClassStep chooseImageAndClassStep(expCfg,expImages,iter);
    const int iclassStep = chooseImageAndClassStep.iclassStep();
    const int iimageStep = chooseImageAndClassStep.iimageStep();
    const int irotStep   = chooseImageAndClassStep.irotStep();

    ChooseImageAndClassStepForKNL chooseImageAndClassStepForKNL(expCfg,expImages,iter,__LINE__);

    TUNING_SCOPE_PARA_BEGIN(getAllSquaredDifferences,chooseImageAndClassStepForKNL.numberOfIterations())

	// Didn't get enough parallelism for KNL without the irotBegin loop on ring11 - only getting 150 steps on a 60 core system, 
	// would prefer 600 for balancing
	//
	#pragma omp parallel for collapse(3) schedule(dynamic)
    for (int iclassBegin = iclass_min; iclassBegin <= iclass_max; iclassBegin += iclassStep) {
      for (int iimageBegin = 0; iimageBegin < nr_images; iimageBegin += iimageStep) {
		  for (int irotBegin = 0; irotBegin < nr_rot; irotBegin += irotStep) {
				TUNING_SCOPE_ITER_BEGIN(getAllSquaredDifferences)

#if defined(LOADBALANCEANALYZER)
				lba.iterationBegin();
#endif

				const int iclassEnd = std::min(iclass_max + 1, iclassBegin + iclassStep);
				for (int iclass = iclassBegin; iclass < iclassEnd; iclass++) {
				    const int iimageEnd = std::min(nr_images, iimageBegin + iimageStep);

				    for (int iimage = iimageBegin; iimage < iimageEnd; iimage++) {

				        const bool interesting_class = false    // using this for debugging
				            && interesting_iter == iter
				            && interesting_ipass == expCfg.ipass()
				            && interesting_iimage == iimage
				            && interesting_iclass == iclass;
				        if (interesting_class) {
				            std::cerr << "interesting_class" << std::endl;
				        }


				        const auto mweight_for_class = exp_Mweight.ptrForWrite(iimage, iclass, 0);
				        // This makes this [iimage, iclass] one valid but does not set the values

				        if (checkForUnsetUsed) exp_Mweight.unsetAll(iimage, iclass, undefWeight);
				        // NOTE : This was the hottest function on KNL, but it is not needed

				        const int  tid = omp_get_thread_num();
				        const auto td  = threadDataForEachThread[tid];

				        const auto Minvsigma2_aux = exp_local_Minvsigma2s[iimage].rptrAll();

				        struct DeleteOnScopeExit {
				            Map2dOptimizer_Kernel::GetAllSquaredDifferences_Kernel* p;
				            ~DeleteOnScopeExit() {
				                p->release();
				                p = NULL;
				            }
				        } kernel;
				        kernel.p =
				            Map2dOptimizer_Kernel::GetAllSquaredDifferences_Kernel::acquire(
								iter,
								nr_rot   * nr_over_rot,
								nr_trans * nr_over_trans,
								exp_current_Fsize2,
								Minvsigma2_aux);

				        if (interesting_class) {
				            if (interesting_irot >= nr_rot || interesting_iover_rot >= nr_over_rot) {
				                std::cerr << "interesting_irot >= " << nr_rot << " or " << " interesting_iover_rot >= " << nr_over_rot << std::endl;
				                std::cerr << "exp_Mweight.sizePerClass() = " << exp_Mweight.sizePerClass() << std::endl;
				                std::cerr << "nr_rot*nr_over_rot*nr_trans*nr_over_trans = " << nr_rot*nr_over_rot*nr_trans*nr_over_trans << std::endl;
				            }
				        }

						int irotEnd = std::min(nr_rot, irotBegin + irotStep);
				        for (int irot = irotBegin; irot < irotEnd; irot++) // 36 but spread over 4, so 8
				        {
				            for (int iover_rot = 0; iover_rot < nr_over_rot; iover_rot++) //1,2
				            {
				                const bool interesting_rot = interesting_class
				                    && interesting_irot == irot
				                    && interesting_iover_rot == iover_rot;
				                if (interesting_rot) {
				                    std::cerr << "interesting_rot" << std::endl;
				                }

				                // exp_iclass loop does not always go from 0 to nr_classes!
				                int iclass_rot = iclass * nr_rot + irot;
				                int iclass_rot_over = iclass_rot * nr_over_rot + iover_rot;
				                int        rot_over = irot * nr_over_rot + iover_rot;

				                if (!exp_Rot_significant[iclass].rptrAll()[irot] || mlModel.pdf_direction[iclass].rptrAll()[0] <= 0. || mlModel.pdf_class.rptrAll()[iclass] <= 0.) {
				                    if (interesting_rot) {
				                        std::cerr << "setting all mweight_for_class to " << unsetWeight << std::endl;
				                    }
				                    int rot_trans_over_begin = (rot_over*nr_trans + 0)*nr_over_trans + 0;
				                    int rot_trans_over_end = rot_trans_over_begin + nr_trans*nr_over_trans;
				                    for (int rot_trans_over = rot_trans_over_begin; rot_trans_over < rot_trans_over_end; rot_trans_over++) {
				                        mweight_for_class[rot_trans_over] = unsetWeight;
				                        // This sets all the entries for [iimage, iclass] [rotation]
				                    }
				                    continue;
				                }

				                SOAComplexReadonly exp_Frefs_Rot_aux;
				                exp_Frefs_Rot_aux.real = expIterData.FrefsAndWeight.frefs_Rot_real_readonly(__LINE__, iclass_rot_over, exp_current_Fsize2);
				                exp_Frefs_Rot_aux.imag = expIterData.FrefsAndWeight.frefs_Rot_imag_readonly(__LINE__, iclass_rot_over, exp_current_Fsize2);

				                SOAComplexReadonly Frefctf_aux;

				                // Apply CTF to reference projection
				                //after first iteration refs_are_ctf_corrected = true
				                if (do_ctf_correction && refs_are_ctf_corrected && (!exp_local_Fctfs_readonly.allOnes())) {
				                    const double* exp_local_Fctfs_aux = exp_local_Fctfs_readonly[iimage].rptrAll();

				                    double* outReal;
				                    double* outImag;
				                    Frefctf_aux.real = outReal = td->Frefctf_real(rot_over);
				                    Frefctf_aux.imag = outImag = td->Frefctf_imag(rot_over);
				                    for (int n = 0; n < exp_current_Fsize2; n++) {
				                        outReal[n] = checkFrefCtfValue(exp_Frefs_Rot_aux.real[n] * exp_local_Fctfs_aux[n]);
				                        outImag[n] = checkFrefCtfValue(exp_Frefs_Rot_aux.imag[n] * exp_local_Fctfs_aux[n]);
				                    }
				                } else {
				                    // There was no reason to copy these elements when the pointers can be copied
				                    Frefctf_aux.real = exp_Frefs_Rot_aux.real;
				                    Frefctf_aux.imag = exp_Frefs_Rot_aux.imag;
				                    for (int n = 0; n < exp_current_Fsize2; n++) {
				                        /* Frefctf_aux.real[n] = */ checkFrefCtfValue(exp_Frefs_Rot_aux.real[n]);
				                        /* Frefctf_aux.imag[n] = */ checkFrefCtfValue(exp_Frefs_Rot_aux.imag[n]);
				                    }
				                }


				                for (int itrans = 0; itrans < nr_trans; itrans++)//100 or 400
				                {
				                    for (int iover_trans = 0; iover_trans < nr_over_trans; iover_trans++)//1,4
				                    {
				                        const bool interesting = interesting_rot
				                            && interesting_itrans == itrans
				                            && interesting_iover_trans == iover_trans
				                            ;
				                        if (interesting) {
				                            std::cerr << "interesting" << std::endl;
				                        }
				                        typedef void interesting_rot;

				                        int rot_trans_over = (rot_over*nr_trans + itrans)*nr_over_trans + iover_trans;

				                        // In the first pass, always proceed
				                        // In the second pass, check whether this translations (&orientation) had a significant weight in the first pass
				                        if (expCfg.ipass() != 0 && !exp_Mcoarse_significant.get(iimage, iclass, irot, itrans)) {
				                            if (interesting) {
				                                std::cerr << "setting  mweight_for_class for this iover_trans to " << unsetWeight << std::endl;
				                            }
				                            mweight_for_class[rot_trans_over] = unsetWeight;
				                            continue;
				                        }

				                        // Get the shifted image
				                        int ishift = iimage * nr_trans * nr_over_trans + itrans * nr_over_trans + iover_trans;

				                        if (expIterData.Fimgs_shifted.efsDeferred(ishift, exp_current_Fsize2)) {
				                            static int count;
				                            if (count++ < 10) {
				                                std::cerr << "getAllSquaredDifferences_newer failed to undefer an exp_Fimgs_shifted" << std::endl;
				                            }
				                        }

				                        SOAComplexReadonly Fimg_shift_aux;
				                        Fimg_shift_aux.real = expIterData.Fimgs_shifted.efsRealConst(__LINE__, ishift, exp_current_Fsize2);
				                        Fimg_shift_aux.imag = expIterData.Fimgs_shifted.efsImagConst(__LINE__, ishift, exp_current_Fsize2);

				                        if (interesting) {
				                            std::cerr << "kernel appending this iover_trans" << std::endl;
				                        }

				                        kernel.p->append(
				                            irot   * nr_over_rot + iover_rot,
				                            itrans * nr_over_trans + iover_trans,
				                            rot_trans_over,
				                            Frefctf_aux.real, Frefctf_aux.imag,
				                            Fimg_shift_aux.real, Fimg_shift_aux.imag,
				                            interesting);

				                    } // end loop iover_trans
				                } // end loop itrans
				            }// end loop iover_rot
				        } // end loop irot

				        kernel.p->compute();

				        double kernel_diff2;
				        int    rot_trans_over;
				        for (int i = 0; kernel.p->getI(kernel_diff2, rot_trans_over, i); i++) {
				            auto& shouldBe = mweight_for_class[rot_trans_over];
				            double diff2 = (kernel_diff2 + exp_highres_Xi2_imgs[iimage]) / 2;
				            shouldBe = diff2;
				        }

				    }//end loop iimage
				} // end loop iclass

#if defined(LOADBALANCEANALYZER)
				lba.iterationEnd();
#endif
				TUNING_SCOPE_ITER_END
			} // end loop irotBegin
        } // end loop iimageBegin
    } // end loop iclassBegin
    TUNING_SCOPE_PARA_END

	doMweightTestOutput();

    TUNING_SCOPE_PARA_BEGIN(getAllSquaredDifferences3,nr_images)
	#pragma omp parallel for
    for (int iimage = 0; iimage < nr_images; iimage++) {
        TUNING_SCOPE_ITER_BEGIN(getAllSquaredDifferences)
        exp_min_diff2[iimage] = (std::numeric_limits<double>::max)();

        for (int iclass = iclass_min; iclass < iclass_max+1; iclass++) {

            for (int withinClass = 0; withinClass < exp_Mweight.sizePerClass(); withinClass++) {
                double diff2 = exp_Mweight.get(iimage, iclass, withinClass, true);

                if (checkForUnsetUsed && diff2 == undefWeight) {

                    int rot_trans_over          = withinClass;                          // (rot_over*nr_trans + itrans)*nr_over_trans + iover_trans
                    int iover_trans             = rot_trans_over % nr_over_trans;
                    int rot_over_trans_itrans   = rot_trans_over / nr_over_trans;   // (rot_over*nr_trans + itrans)
                    int itrans                  = rot_over_trans_itrans % nr_trans;
                    int rot_over                = rot_over_trans_itrans / nr_trans;
                    int iover_rot               = rot_over % nr_over_rot;
                    int irot                    = rot_over / nr_over_rot;

                    // std::cerr << "iover_trans            :" << iover_trans           << " = " << rot_trans_over        << "%" << nr_over_trans   << std::endl;
                    // std::cerr << "rot_over_trans_itrans :" << rot_over_trans_itrans << " = " << rot_trans_over         << "/" << nr_over_trans   << std::endl;
                    // std::cerr << "itrans             :" << itrans                << " = " << rot_over_trans_itrans << "%" << nr_trans        << std::endl;
                    // std::cerr << "rot_over               :" << rot_over              << " = " << rot_over_trans_itrans << "/" << nr_trans        << std::endl;
                    // std::cerr << "iover_rot             :" << iover_rot             << " = " << rot_over           << "%" << nr_over_rot << std::endl;
                    // std::cerr << "irot                  :" << irot                  << " = " << rot_over           << "/" << nr_over_rot << std::endl;
                    //
                    std::cerr << "undefWeight fetched from exp_Mweight"
                        << " iter:"        << iter
                        << " ipass:"       << expCfg.ipass()
                        << " iimage:"      << iimage
                        << " iclass:"      << iclass
                        << " irot:"        << irot
                        << " iover_rot:"   << iover_rot
                        << " itrans:"      << itrans
                        << " iover_trans:" << iover_trans
                        << " withinClass:" << withinClass
                        << " diff2:"       << diff2
                        << " next value:"  << exp_Mweight.get(iimage, iclass, withinClass+1, true)
                        << std::endl;
                    EXIT_ABNORMALLY;
                }

                if (diff2 > 0. && exp_min_diff2[iimage] > diff2) {
                    exp_min_diff2[iimage] = diff2;
                }
            }
        }

        if (emit_test_output_prolific() && iimage == 0)
        #pragma omp critical
        {
            testos() << "~~~~TEST OUTPUT: exp_min_diff2[iimage:" << iimage << "]:"
                << exp_min_diff2[iimage]
                << std::endl;
        }
        TUNING_SCOPE_ITER_END
    }
    TUNING_SCOPE_PARA_END
    ml_original_waypoint("getAllSquaredDifferences after");

}


void findAllSignificantPoints()
{   TUNING_SCOPE_STEP(findAllSignificantPoints)

    ml_original_waypoint("findAllSignificantPoints before");

    const int sorted_weight_len = exp_Mweight.nr_classes()*exp_Mweight.sizePerClass();
    std::vector<Exp_Mweight::Elt*> sorted_weight_per_thread(omp_get_max_threads());
    for (int i = 0; i < sorted_weight_per_thread.size(); i++) sorted_weight_per_thread[i] = NULL;
        // Shift the malloc outside the loop and share it to reduce contention on 60+ core systems

    // For each image,find the exp_significant_weight that encompasses adaptive_fraction of exp_sum_weight
    TUNING_SCOPE_PARA_BEGIN(findAllSignificantPoints, expImages.nr_images())
	auto const nr_images = expImages.nr_images();
	#pragma omp parallel for if (!emit_test_output_prolific())
    for (int iimage = 0; iimage < nr_images; iimage++)
    {
        TUNING_SCOPE_ITER_BEGIN(findAllSignificantPoints)
        auto& sorted_weight_ptr = sorted_weight_per_thread[omp_get_thread_num()];
        if (!sorted_weight_ptr) {
            sorted_weight_ptr = mallocCacheAligned(Exp_Mweight::Elt,sorted_weight_len);
        }

        auto sorted_weight       = sorted_weight_ptr;
        int  sorted_weight_count = 0;

        double frac_weight = 0.;
        int my_nr_significant_coarse_samples = 0;

        for (int iclass = 0; iclass < exp_Mweight.nr_classes(); iclass++) {
            if (!exp_Mweight.isSet(iimage,iclass)) continue;
            for (int withinClass = 0; withinClass < exp_Mweight.sizePerClass(); withinClass++) {
                auto w = exp_Mweight.get(iimage,iclass,withinClass,true);
                if (w > 0.0) {
                    assert(sorted_weight_count < sorted_weight_len);
                    sorted_weight[sorted_weight_count] = w;
                    sorted_weight_count++;
                }
            }
        }
        assert(sorted_weight_count > 0);

        //sort from small to large
        std::sort(sorted_weight, sorted_weight + sorted_weight_count);

        if (emit_test_output_prolific() && iimage == 0)
        #pragma omp critical
        {
            testos() << "~~~~TEST OUTPUT: sorted_weights[iimage] "
                << " iimage:" << iimage << " sorted_weight_count:" << sorted_weight_count
                << std::endl;
            for (int i = sorted_weight_count - 1; i >= 0; i--) {
                if (i > 10 && i < sorted_weight_count - 10) i = 10; // skip the middle ones
                testos() << "~~~~    i:" << i << " => " << sorted_weight[i] << std::endl;
            }
        }

        double my_significant_weight;
        for (int i = sorted_weight_count - 1; i >= 0; i--)
        {
            if (expCfg.ipass()==0) my_nr_significant_coarse_samples++;
            my_significant_weight = double(sorted_weight[i]);
            frac_weight += my_significant_weight;
            if (frac_weight > adaptive_fraction * exp_sum_weight[iimage])
                break;
        }
        if (my_significant_weight == 0.0) {
            static int count; static const int limit = 5;
            if (count++ < limit) {
                std::cerr << "my_significant_weight = 0.0 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
                std::cout << "my_significant_weight = 0.0 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
            }
            if (count == limit) {
                std::cout << "Suppressing future reports" << std::endl;
                std::cerr << "Suppressing future reports" << std::endl;
            }
        }

        if (expCfg.ipass()==0 && my_nr_significant_coarse_samples == 0)
        #pragma omp critical
        {
            std::cout << " iimage= " << iimage << " adaptive_fraction= " << adaptive_fraction << std::endl;
            std::cout << " frac-weight= " << frac_weight << std::endl;
            std::cout << " exp_sum_weight[iimage] = " << exp_sum_weight[iimage] << std::endl;
            if (exp_sum_weight[iimage] < 0) std::cout << "     *** NOTE: Negative!" << std::endl;
            std::cout << " sorted_weight_count = " << sorted_weight_count << std::endl;
            if (sorted_weight_count > 0)
            {
                std::cout<<"sum of abs(sorted_weight) = "<<sumvec(sorted_weight, sorted_weight_count)<<std::endl;
                std::cout << "written sorted_weight.spi" << std::endl;
            }
            ERROR_REPORT("my_nr_significant_coarse_samples == 0");
        }

        if (expCfg.ipass()==0)
        {
            exp_metadata[iimage].NR_SIGN = (double)my_nr_significant_coarse_samples;
            // Keep track of which coarse samplings were significant for this particle
            //in ipass = 0,exp_Mcoarse_significant_xsize equal to exp_Mweight_xsize
            exp_Mcoarse_significant.set(exp_Mweight, iimage, my_significant_weight);
        }

        exp_significant_weight[iimage] = my_significant_weight;
        if (emit_test_output_prolific())
        #pragma omp critical
        {
            testos() << "~~~~TEST OUTPUT: exp_significant_weight[iimage] "
                << " iimage:" << iimage << " my_significant_weight:" << my_significant_weight
                << std::endl;
        }

        TUNING_SCOPE_ITER_END
    } // end loop iimage
    TUNING_SCOPE_PARA_END

    for (int i = 0; i < sorted_weight_per_thread.size(); i++) {
        auto & ptr = sorted_weight_per_thread[i];
        if (!ptr) continue;
        Heap::freeDoubles(ptr);
    }

    ml_original_waypoint("findAllSignificantPoints after");

}


static double EncodeIntoCLASS(int iclass, int withinClass) {
    return iclass * expCfg.nr_rot()*expCfg.nr_over_rot()*expCfg.nr_trans()*expCfg.nr_over_trans() + withinClass;
}

static void DecodeFromCLASS(double CLASS, int & iclass, int & ipsi, int& iover_rot, int & itrans, int & iover_trans) {
    int ihidden_over = CLASS;
    int denominator  = expCfg.nr_rot()*expCfg.nr_over_rot()*expCfg.nr_trans()*expCfg.nr_over_trans();
    iclass = ihidden_over / denominator;

    ihidden_over = ihidden_over % denominator;
    denominator = expCfg.nr_over_rot()*expCfg.nr_trans()*expCfg.nr_over_trans();
    ipsi = ihidden_over / denominator;

    ihidden_over = ihidden_over % denominator;
    denominator = expCfg.nr_trans()*expCfg.nr_over_trans();
    iover_rot = ihidden_over / denominator;

    ihidden_over = ihidden_over % denominator;
    denominator = expCfg.nr_over_trans();
    itrans = ihidden_over / denominator;
    iover_trans = ihidden_over % denominator;
}

void convertSquaredDifferencesToWeights()
{
	auto const iclass_min	 = expImages.iclass_min(); 
	auto const iclass_max	 = expImages.iclass_max();
    auto const nr_rot		 = expCfg.nr_rot();
	auto const nr_over_rot	 = expCfg.nr_over_rot();
    auto const nr_images	 = expImages.nr_images(); 
    auto const nr_trans		 = expCfg.nr_trans();
	auto const nr_over_trans = expCfg.nr_over_trans();

    TUNING_SCOPE_STEP(convertSquaredDifferencesToWeights)
    ml_original_waypoint("convertSquaredDifferencesToWeights before");

    {
#if defined(LOADBALANCEANALYZER)
        LoadBalanceAnalyzer lba("convertSquaredDifferencesToWeights - new way", __FILE__, __LINE__, nr_images*(iclass_max+1-iclass_min));
#endif

        auto showSituation = [&](std::ostream & os) {
            os << "convertSquaredDifferencesToWeights - new way"
#define P(N)    << " " << #N << ":" << N
                P(iclass_min) P(iclass_max) P(nr_images) P(omp_get_max_threads())
#undef P
                << std::endl;
        };
        if (false) showSituation(std::cout);
        if (false) showSituation(std::cerr);

        TUNING_SCOPE_PARA_BEGIN(convertSquaredDifferencesToWeights,nr_images*(iclass_max+1-iclass_min))
		#pragma omp parallel for collapse(2) schedule(dynamic) if (!emit_test_output_prolific())
        for (int iimage = 0; iimage < nr_images; iimage++)//nr_pool
        {
            for (int iclass = iclass_min; iclass <= iclass_max; iclass++)//50 or 100
            {
                TUNING_SCOPE_ITER_BEGIN(convertSquaredDifferencesToWeights)
#if defined(LOADBALANCEANALYZER)
                lba.iterationBegin();
#endif
                float pdf_orientation = mlModel.pdf_direction[iclass][0];
                auto exp_min_diff2_iimage = exp_min_diff2[iimage];

                size_t testOutputsDone = 0;
				size_t count_of_nonzero_pdf_offsets = 0;

                for (int itrans = 0; itrans < nr_trans; itrans++)//100 or 400
                {
					const bool doTestOutput = 
						(emit_test_output_prolific() 
						 && iimage == 0 
						 && testOutputsDone++ < 20);

					FDOUBLE offset_x,offset_y;
                    sampling2d.getTranslation(itrans, offset_x, offset_y);

                    float pdf_offset =
                        mlModel.calculatePdfOffset(
                            	exp_old_offsetx[iimage] + offset_x,
                            	exp_old_offsety[iimage] + offset_y,
                            	mlModel.prior_offsetx_class.rptrAll()[iclass],
                            	mlModel.prior_offsety_class.rptrAll()[iclass]);
					
					if (pdf_offset != 0.0) count_of_nonzero_pdf_offsets++;

                    const float unadjustedWeight = pdf_orientation * pdf_offset;

                    if (unadjustedWeight < 0 || isnan(unadjustedWeight)) {
                        std::cerr << "unadjustedWeight bad:" << unadjustedWeight << std::endl;
                        std::cout << "unadjustedWeight bad:" << unadjustedWeight << std::endl;
                        EXIT_ABNORMALLY;
                    }

                    auto inner = [&](int known_nr_over_rot, int known_nr_over_trans) {
                        for (int irot = 0; irot < expCfg.nr_rot(); irot++)//36
                        {
                            for (int iover_rot = 0; iover_rot < known_nr_over_rot; iover_rot++)//1,2
                            {
                                int rot_over         = irot*expCfg.nr_over_rot() + iover_rot;
                                int rot_over_trans_0 = (rot_over*nr_trans + itrans)*nr_over_trans + 0;
                                auto mweight         = exp_Mweight.ptrForModify(iimage,iclass, rot_over_trans_0);

                                for (int iover_trans = 0; iover_trans < known_nr_over_trans; iover_trans++)//1,4
                                {
                                    auto & mw = mweight[iover_trans];
                                    auto mwBefore = mw;

                                    // Only exponentiate for determined values of exp_Mweight
                                    // (this is always true in the first pass, but not so in the second pass)
                                    // Only deal with this sampling point if its weight was significant
                                    if (mw < 0.0) {
                                        mw = 0.0;
                                    } else {
                                        double diff2 = mw - exp_min_diff2_iimage;
                                        // next line because of numerical precision of exp-function

                                        if (diff2 < 0.0 || isnan(diff2)) {
                                            std::cerr << "convertSquaredDifferencesToWeights diff2 bad, is " 
												<< diff2 << " = mw:" << mw << " - exp_min_diff2_iimage:" << exp_min_diff2_iimage << std::endl;
                                            std::cout << "convertSquaredDifferencesToWeights diff2 bad, is " 
												<< diff2 << " = mw:" << mw << " - exp_min_diff2_iimage:" << exp_min_diff2_iimage << std::endl;
                                            EXIT_ABNORMALLY;
                                        }
                                        else if (diff2 > 700.0) mw = 0.0;
                                        // TODO: use tabulated exp function?
                                        else mw = unadjustedWeight * exp(-diff2);
                                    }

									if (doTestOutput && irot == 0 && iover_rot == 0 && iover_trans == 0)
									#pragma omp critical
									{
									    testos() << "~~~~TEST OUTPUT: "
									        << " iimage:" << iimage
									        << " iclass:" << iclass
									        << " itrans:" << itrans
											<< " getWeight:" << mwBefore
									        << std::endl;
										if (mw != 0.0 && mw != Exp_Mweight_unsetWeight) {
											testos() << "~~~~TEST OUTPUT: pdf_offset:" << pdf_offset
												<< " iimage:" << iimage
												<< " iclass:" << iclass
												<< " itrans:" << itrans
												<< " weight:" << mw
												<< std::endl;
										}
									}
                                }
                            }
                        }
                    };

                    // really help the compiler unroll those inner loops
                    if (nr_over_rot == 1 && nr_over_trans == 1) inner(1, 1); else
                    if (nr_over_rot == 2 && nr_over_trans == 1) inner(2, 1); else
                    if (nr_over_rot == 1 && nr_over_trans == 4) inner(1, 4); else
                    if (nr_over_rot == 2 && nr_over_trans == 4) inner(2, 4); else {
                        std::cerr << "convertSquaredDifferencesToWeights did not optimize expCfg.nr_over_rot():" << expCfg.nr_over_rot() << " nr_over_trans:" << nr_over_trans << std::endl;
                        inner(expCfg.nr_over_rot(), nr_over_trans);
                    }
                }

                if (count_of_nonzero_pdf_offsets == 0) {
                    std::cout << "all pdf offset are zero!" << std::endl;
                    ERROR_REPORT("ERROR: all pdf offset are zero!");
                }

#if defined(LOADBALANCEANALYZER)
                lba.iterationEnd();
#endif
                TUNING_SCOPE_ITER_END
            } //end loop iclass
        }//end loop iimage
		TUNING_SCOPE_PARA_END
    }

	doMweightTestOutput();


    //get sum of weight and maximum weight for each image,when ipass=0,the max_weight may not need
    // NOTE : change the order of weight add-operation may get different of sum_weight,because the
    // weight is much small and the add-operation take in some error.
    // this is why my baseline testing case:data6 get different .NR_SIGN and .DLL..................
    TUNING_SCOPE_PARA_BEGIN(convertSquaredDifferencesToWeights,nr_images)
#pragma omp parallel for if (!emit_test_output_prolific())
    for (int iimage = 0; iimage < nr_images; iimage++)
    {
        TUNING_SCOPE_ITER_BEGIN(convertSquaredDifferencesToWeights)
        double max_weight             = 0;
        int    max_weight_iclass      = 0;
        int    max_weight_withinClass = 0;

        const int usedPortionOfSizePerClass = expCfg.nr_rot()*expCfg.nr_over_rot()*nr_trans*nr_over_trans;
        assert(exp_Mweight.sizePerClass() >= usedPortionOfSizePerClass);

        double sum_weight = 0.0;
        for (int iclass = iclass_min; iclass < iclass_max+1; iclass++) {
            for (int withinClass = 0; withinClass < usedPortionOfSizePerClass; withinClass++) {

                double weight = exp_Mweight.get(iimage,iclass,withinClass);
                if (weight < 0.0) {
                    std::cerr << "weight < 0! "
#define P(X) << " " << #X << ":" << X
                        P(weight)
                        P(iimage)
                        P(iclass) P(iclass_min) P(iclass_max)
                        P(withinClass) P(expCfg.nr_rot()) P(expCfg.nr_over_rot()) P(nr_trans) P(nr_over_trans)
#undef P
                        << std::endl;
					EXIT_ABNORMALLY;
                }
                sum_weight += weight;
                if (isnan(sum_weight)) {
                    std::cerr << "weight isnan, added " << weight << std::endl;
					EXIT_ABNORMALLY;
                }
                //notice! if some maximum weight is equal,than choice the first one
                if (max_weight < weight) {
                    max_weight = weight;
                    max_weight_iclass = iclass;
                    max_weight_withinClass = withinClass;
                }
            }
        }

        exp_sum_weight[iimage] = sum_weight;

        if (emit_test_output_prolific()) {
            testos() << "convertSquaredDifferencesToWeights exp_sum_weight[iimage:"<<iimage<<"]=" << exp_sum_weight[iimage]
                << " max_weight:" << max_weight << std::endl;
        }
        if (sum_weight <= 0.0) {
            std::cout << "convertSquaredDifferencesToWeights sum_weight<=0"
#define P(X)    << " " << #X << ":" << X
                P(sum_weight) P(iclass_min) P(iclass_max) P(usedPortionOfSizePerClass) P(max_weight)
#undef P
                << std::endl;
            int count = 0;
            for (int iclass = iclass_min; iclass < iclass_max+1; iclass++) {
                for (int withinClass = 0; count < 10 && withinClass < usedPortionOfSizePerClass; withinClass++) {
                    auto w = exp_Mweight.get(iimage,iclass,withinClass);
                    if (w >= 0.0) continue;
                    count++;
                    std::cout << "exp_Mweight[" << iclass <<"][" << withinClass << "]:" << w << std::endl;
                }
            }
        }

        // TODO - what should 0/0 be here?  This is a guess...
        // YongBei : need to debug to find out why this is zero...
        if (max_weight == 0) max_weight = sum_weight = 1.0;

        exp_metadata[iimage].PMAX  = max_weight/sum_weight;
        //trick!code the position index in metadata CLASS,the position information will
        //decode in StoreWeightedSumsAllOrientations()
        exp_metadata[iimage].CLASS = EncodeIntoCLASS(max_weight_iclass, max_weight_withinClass);
        TUNING_SCOPE_ITER_END
    }
    TUNING_SCOPE_PARA_END

    ml_original_waypoint("convertSquaredDifferencesToWeights after");

}


void __declspec(noinline) storeWeightedSums()
{   TUNING_SCOPE_STEP(storeWeightedSums)

    int current_Fsize = current_size/2+1;
    int current_Fsize2 = current_size*(current_size/2+1);
    int ori_Fsize = ori_size/2+1;

    // Extend norm_correction and sigma2_noise estimation to higher resolutions
	auto const nr_images = expImages.nr_images();
    TUNING_SCOPE_PARA_BEGIN(sws,nr_images)
	#pragma omp parallel for
    for (int iimage = 0; iimage < nr_images; iimage++) {

        int igroup = exp_metadata[iimage].GROUP_NO-1;
        auto exp_power_imgs_iimage = exp_power_imgs[iimage].rptrAll();
        auto wsum_sigma2_noise_igroup = mlModel.wsum_sigma2_noise[igroup].mptrAll();
        // If the current images were smaller than the original size, fill the rest of wsum_model.sigma2_noise with the power_class spectrum of the images
        for (int ires = current_Fsize; ires < ori_Fsize; ires++) {
            wsum_sigma2_noise_igroup[ires] += exp_power_imgs_iimage[ires];
            // Also extend the weighted sum of the norm_correction
            exp_wsum_norm_correction[iimage] += exp_power_imgs_iimage[ires];
        }

        // Store norm_correction
        // Multiply by old value because the old norm_correction term was already applied to the image
        if (do_norm_correction)
        {
            double old_norm_correction = exp_metadata[iimage].NORM / mlModel.avg_norm_correction;
            // Now set the new norm_correction in the relevant position of exp_metadata22
            // The factor two below is because exp_wsum_norm_correctiom is similar to sigma2_noise, which is the variance for the real/imag components
            // The variance of the total image (on which one normalizes) is twice this value!
            exp_metadata[iimage].NORM = old_norm_correction * sqrt(exp_wsum_norm_correction[iimage] * 2.);
            if (exp_metadata[iimage].NORM > 10.)
            {
#pragma omp critical
                std::cerr<<"Warning in storeWeightedSums(),please debug this function."<<std::endl;
            }
#pragma omp atomic
            mlModel.wsum_avg_norm_correction += old_norm_correction * sqrt(exp_wsum_norm_correction[iimage] * 2.);
        }

        // Some analytics...
        // Calculate normalization constant for dLL
        // loop over all particles inside this ori_particle
        double logsigma2 = 0.;
        auto sigma2_noise_igroup = mlModel.sigma2_noise[igroup].mptrAll();
        for (int n = 0; n < current_Fsize2; n++)
        {
            int ires = Mresol_fine[n];
            // Note there is no sqrt in the normalisation term because of the 2-dimensionality of the complex-plane
            // Also exclude origin from logsigma2, as this will not be considered in the P-calculations
            if (ires > 0)
                logsigma2 += CHECK_NOT_NAN_OR_INF(log( 2. * PI * sigma2_noise_igroup[ires]));
        }


        if (exp_sum_weight[iimage]<0)
        {
            ERROR_REPORT("ERROR: exp_sum_weight[iimage]<=0");
        }

        // TODO what to do when exp_sum_weight[iimage] == 0
        // YONGBEI : it will be bug if exp_sum_weight[iimage] == 0
        double dLL = 0.0;
        if (exp_sum_weight[iimage] > 0) {
            dLL = log(exp_sum_weight[iimage]) - exp_min_diff2[iimage] - logsigma2;
        }
        if (emit_test_output_prolific()) std::cout << "dLL:" << dLL << std::endl;

        // Also store dLL of each image in the output array
        exp_metadata[iimage].DLL = dLL;

#pragma omp atomic
        mlModel.wsum_LL += dLL;
#pragma omp atomic
        mlModel.wsum_ave_Pmax += exp_metadata[iimage].PMAX;

    } // end loop iimage
    TUNING_SCOPE_PARA_END

    ml_original_waypoint("storeWeightedSums after");

}


void __declspec(noinline) updateOtherParams()
{   TUNING_SCOPE_STEP(updateOtherParams)

	auto const iclass_min	 = expImages.iclass_min(); 
	auto const iclass_max	 = expImages.iclass_max();
    auto const nr_rot		 = expCfg.nr_rot();
	auto const nr_over_rot	 = expCfg.nr_over_rot();
    auto const nr_images	 = expImages.nr_images(); 
    auto const nr_trans		 = expCfg.nr_trans();
	auto const nr_over_trans = expCfg.nr_over_trans();

    ml_original_waypoint("updateOtherParams before");

    const int exp_current_Fsize2 = expCfg.current_size()*(expCfg.current_size()/2+1);

    fillvec(exp_wsum_norm_correction, 0);

    threadDataForEachThread.fini_wsum_norm_correction_sigma_noise();
        // These should be fini'ed from any previous use, and only init'ed if the tid appears below

{
#if defined(LOADBALANCEANALYZER)
    #define LOADBALANCEANALYZER_updateOtherParams
#else
    // #define LOADBALANCEANALYZER_updateOtherParams
#endif

#if defined(LOADBALANCEANALYZER_updateOtherParams)
    LoadBalanceAnalyzer lba("updateOtherParams", __FILE__, __LINE__, nr_images*(iclass_max+1-iclass_min));
#endif

#pragma omp parallel
{   TUNING_SCOPE_PARA_BEGIN(updateOtherParams, (iclass_max+1-iclass_min)*nr_rot*nr_over_rot())

#if defined(LOADBALANCEANALYZER_updateOtherParams)
    lba.note();
#endif
    const int tid = omp_get_thread_num();
    auto td = threadDataForEachThread[tid];

    // The following loop computes this data
    //
    td->zero_wsum_norm_correction_sigma_noise();

    auto newWay_td_wsum_norm_correction = td->wsum_norm_correction();
    auto newWay_td_wsum_sigma2_noise    = td->wsum_sigma2_noise();

    if (nr_images*(iclass_max+1-iclass_min) < 4*omp_get_max_threads())
    #pragma omp critical
    if (false) {
        static int count;
        if (count++ < 2) {
            std::cerr << "UpdateOtherParams - not enough iterations for good parallelism" << std::endl
                #define P(N) << " " << #N << ":" << N
                P(nr_images) P(iclass_max) P(iclass_min) P(omp_get_max_threads())
                #undef P
                << std::endl;
        }
    }

    // The kernel can only process one image at a time, so its loop must have images as the outermost case
    // The transforms want to process multiple images but one transform, so its loop must have the images as the innermost loop

    // FIRST LOOP PREDICTS WHETHER A KERNEL APPEND WILL HAPPEN IN THE SECOND LOOP, and if so does the deferred translations
    // This is surprisingly easy and fast
    // But it may be a waste of time, because apparently getAllSquaredDifferences undefers them first
    //
	{
		// BEVIN ANNOTATE THIS
		#pragma omp for collapse(3) schedule(dynamic)
		for (int iclass = iclass_min; iclass <= iclass_max; iclass++) {
		    for (int irot = 0; irot < nr_rot; irot++) {
		        for (int iover_rot = 0; iover_rot < nr_over_rot; iover_rot++) {

		            int iclass_rot = iclass * nr_rot + irot;
		            int   rot_over = irot * nr_over_rot + iover_rot;

		            if (!exp_Rot_significant[iclass][irot])
		                continue;

		            for (int itrans = 0; itrans < nr_trans; itrans++) {
		                for (int iover_trans = 0; iover_trans < nr_over_trans; iover_trans++) {
		                    int rot_trans_over = (rot_over*nr_trans + itrans)*nr_over_trans+iover_trans;

		                    Exp_Fimgs_shifted::Buffer buffer(&expIterData.Fimgs_shifted, "updateOtherParams", __LINE__, exp_current_Fsize2);

		                    for (int iimage = 0; iimage < nr_images; iimage++) {
		                        // Only deal with this sampling point if its weight was significant
		                        double weight = exp_Mweight.get(iimage,iclass,rot_trans_over);
		                        if (weight >= exp_significant_weight[iimage]) {
		                            int ishift = iimage * nr_over_trans * nr_trans + itrans * nr_over_trans + iover_trans;
		                            if (expIterData.Fimgs_shifted.efsDeferred(ishift, exp_current_Fsize2)) 
										buffer.append("UpdateOtherParams", __LINE__, ishift);
		                        }
		                    }
		                }
		            }
		        }
		    }
		}
	}

	// BEVIN ANNOTATE THIS
	#pragma omp for collapse(2) schedule(dynamic)//----- remove this to produce same result as new one -----
    for (int iimage = 0; iimage < nr_images; iimage++)//nr_pool
    {
        for (int iclass = iclass_min; iclass <= iclass_max; iclass++)//50 or 100
        {   TUNING_SCOPE_ITER_BEGIN(updateOtherParams)
#if defined(LOADBALANCEANALYZER_updateOtherParams)
            lba.iterationBegin();
#endif

			struct DeleteOnScopeExit {
				Map2dOptimizer_Kernel::UpdateOtherParams_Kernel* p;
				~DeleteOnScopeExit() {
					if (!p) return;
					p->release();
					p = NULL;
				}
			}* kernelPtr = nullptr;	// TODO - remove the indirection

			{
				kernelPtr = sNew(DeleteOnScopeExit);

				// double* const Minvsigma2_aux = exp_local_Minvsigma2s[iimage].data(); // not use this

				kernelPtr->p =
					Map2dOptimizer_Kernel::UpdateOtherParams_Kernel::acquire(
						iter,
						nr_rot   * nr_over_rot,
						nr_trans * nr_over_trans,
						exp_current_Fsize2,
						Mresol_fine.mptrAll(), newWay_td_wsum_sigma2_noise,
						nr_rot*nr_over_rot, nr_trans*nr_over_trans);
			}

            for (int irot = 0; irot < nr_rot; irot++)//36
            {
                for (int iover_rot = 0; iover_rot < nr_over_rot; iover_rot++)//2
                {
                    int iclass_rot = iclass * nr_rot + irot;

                    if (!exp_Rot_significant[iclass][irot])
                        continue;

                    int iclass_rot_over = iclass_rot * nr_over_rot + iover_rot;
                    int        rot_over =       irot * nr_over_rot + iover_rot;

                    SOAComplexReadonly exp_Frefs_Rot_aux;

                    exp_Frefs_Rot_aux.real = expIterData.FrefsAndWeight.frefs_Rot_real_readonly(__LINE__, iclass_rot_over, exp_current_Fsize2);
                    exp_Frefs_Rot_aux.imag = expIterData.FrefsAndWeight.frefs_Rot_imag_readonly(__LINE__, iclass_rot_over, exp_current_Fsize2);

                    // Apply CTF to reference
                    SOAComplexReadonly Frefctf_aux;
                    if (do_ctf_correction && refs_are_ctf_corrected && (!exp_local_Fctfs_readonly.allOnes())) {
                        double* outReal;
                        double* outImag;
                        Frefctf_aux.real = outReal = td->Frefctf_real(rot_over);
                        Frefctf_aux.imag = outImag = td->Frefctf_imag(rot_over);

                        auto exp_local_Fctfs_aux = exp_local_Fctfs_readonly[iimage].rptrAll();
                        #pragma ivdep
                        for (int n = 0; n < exp_current_Fsize2; n++) {
                            outReal[n] = exp_Frefs_Rot_aux.real[n] * exp_local_Fctfs_aux[n];
                            outImag[n] = exp_Frefs_Rot_aux.imag[n] * exp_local_Fctfs_aux[n];
                        }
                    } else {
                        Frefctf_aux.real = exp_Frefs_Rot_aux.real;
                        Frefctf_aux.imag = exp_Frefs_Rot_aux.imag;
                    }


                    for (int itrans = 0; itrans < nr_trans; itrans++)//100 or 400
                    {
                        for (int iover_trans = 0; iover_trans < nr_over_trans; iover_trans++)//4
                        {
                            int rot_trans_over = (rot_over*nr_trans + itrans)*nr_over_trans+iover_trans;

                            // Only deal with this sampling point if its weight was significant
                            double weight = exp_Mweight.get(iimage,iclass,rot_trans_over);

                            // Only sum weights for non-zero weights
                            // Must be > to avoid 0 >= 0 which would lead to a divide by zero below
                            // NOTE : Must be >=
                            if (weight >= exp_significant_weight[iimage]) {

                                // Normalise the weight (do this after the comparison with exp_significant_weight!)
                                CHECK_NOT_NAN_OR_INF(weight /= exp_sum_weight[iimage]);

                                // Get the shifted image
                                int ishift = iimage * nr_over_trans * nr_trans + itrans * nr_over_trans + iover_trans;

                                SOAComplexReadonly Fimg_shift_aux;
                                Fimg_shift_aux.real = expIterData.Fimgs_shifted.efsRealConst(__LINE__, ishift,exp_current_Fsize2);
                                Fimg_shift_aux.imag = expIterData.Fimgs_shifted.efsImagConst(__LINE__, ishift,exp_current_Fsize2);

								kernelPtr->p->append(
								    irot   * nr_over_rot   + iover_rot,
								    itrans * nr_over_trans + iover_trans,
								    Frefctf_aux   .real,
								    Frefctf_aux   .imag,
								    Fimg_shift_aux.real,
								    Fimg_shift_aux.imag,
								    weight,
								    false);
                            } // end if weight >= exp_significant_weight
                        } // end loop iover_trans
                    } // end loop itrans

                }// end if iover_rot
            } // end loop iorient


			kernelPtr->p->compute();
			{
			    double kernel_diff2;
			    for (int i = 0; kernelPtr->p->getI(kernel_diff2,i); i++) {
			        // Partially reduce, into a thread private
			        newWay_td_wsum_norm_correction[iimage] += kernel_diff2;
			    }
			}

			sDelete(kernelPtr);

#if defined(LOADBALANCEANALYZER_updateOtherParams)
            lba.iterationEnd();
#endif
            TUNING_SCOPE_ITER_END
        }// end of iclass
    } // end loop iimage

#if defined(LOADBALANCEANALYZER_updateOtherParams)
    lba.note();
#endif

	TUNING_SCOPE_PARA_END
}  // end parallel
}  // end LOADBALANCEANALYZER

    for (int thread = 0; thread < maxthreads; thread++) {
        auto td = threadDataForEachThread[thread];
        if (td->written_wsum_norm_correction_sigma_noise()) {
            auto td_wsum_norm_correction = td->wsum_norm_correction();
            for (int iimage = 0; iimage < nr_images; iimage++) {
                exp_wsum_norm_correction[iimage] += td_wsum_norm_correction[iimage];
            }
        }
        td = threadDataForEachThread[thread];
        if (td->written_wsum_norm_correction_sigma_noise()) {
            // TODO....TODO
            int igroup = exp_metadata[0/*iimage*/].GROUP_NO-1;
            auto wsum_sigma2_noise_igroup = mlModel.wsum_sigma2_noise[igroup].mptrAll();
            auto td_wsum_sigma2_noise = td->wsum_sigma2_noise();
            for (int n = 0; n < exp_current_Fsize2; n++) {
                int ires = Mresol_fine[n];
                if (ires > -1) {
                    CHECK_NOT_NAN_OR_INF(wsum_sigma2_noise_igroup[ires] += td_wsum_sigma2_noise[n]);
                }
            }
        }
    }

    threadDataForEachThread.fini_wsum_norm_correction_sigma_noise();
        // these should not be used after here

    threadDataForEachThread.fini_wsum_pdf_prior_offset_max_weight();
        // These should be fini'ed from any previous use, and only init'ed if the tid appears below

#pragma omp parallel if(!emit_test_output_prolific())
{   TUNING_SCOPE_PARA_BEGIN(updateOtherParams, nr_images*(iclass_max+1-iclass_min)*nr_rot*nr_over_rot())

    int tid = omp_get_thread_num();

    double thread_wsum_sigma2_offset = 0.;
    double thread_sumw_group = 0.;

    threadDataForEachThread[tid]->zero_wsum_pdf_prior_offset_max_weight();

    //before,this openmp region is together with prevous one
    //after,break that one openmp region to two.
    //this could give more parallelism threads

	// BEVIN ANNOTATE THIS
#pragma omp for collapse(4) schedule(dynamic)//no load imbalance,overhead high
    for (int iimage = 0; iimage < nr_images; iimage++)//nr_pool
    {
        for (int iclass = iclass_min; iclass <= iclass_max; iclass++)//50 or 100
        {
            for (int irot = 0; irot < nr_rot; irot++)//36
            {
                for (int iover_rot = 0; iover_rot < nr_over_rot; iover_rot++)//2
                {   TUNING_SCOPE_ITER_BEGIN(updateOtherParams)
                    size_t skipped(0);
                    size_t significant(0);

                    auto td = threadDataForEachThread[tid];
                    auto td_wsum_pdf_class           = td->wsum_pdf_class          ();
                    auto td_wsum_pdf_direction       = td->wsum_pdf_direction      ();
                    auto td_wsum_prior_offsetx_class = td->wsum_prior_offsetx_class();
                    auto td_wsum_prior_offsety_class = td->wsum_prior_offsety_class();

                    auto thread_wsum_sigma2_offset_before = thread_wsum_sigma2_offset;

                    for (int itrans = 0; itrans < nr_trans; itrans++)//100 or 400
                    {
                        for (int iover_trans = 0; iover_trans < nr_over_trans; iover_trans++)//4
                        {
                            int iclass_rot = iclass * nr_rot + irot;

                            if (!exp_Rot_significant[iclass][irot])
                                continue;

                            int rot_over       = irot * nr_over_rot + iover_rot;
                            int rot_trans_over = (rot_over*nr_trans + itrans)*nr_over_trans + iover_trans;

                            // Only deal with this sampling point if its weight was significant
                            double weight = exp_Mweight.get(iimage,iclass,rot_trans_over);

                            // Only sum weights for non-zero weights
                            // This must be <= to catch the 0 <= 0 case
                            // NOTE : this must be <
                            if (weight < exp_significant_weight[iimage]) {
                                skipped++;
                            } else {
                                significant++;

                                if (exp_sum_weight[iimage] == 0.0) {
                                    std::cout << "exp_sum_weight[iimage] == 0 yet signicant because exp_significant_weight[iimage]:"
                                        << exp_significant_weight[iimage]
                                        << std::endl;
                                    EXIT_ABNORMALLY;
                                }

                                // Normalise the weight (do this after the comparison with exp_significant_weight!)
                                weight /= exp_sum_weight[iimage];
                                assert(weight >= 0.0);

                                // Store sum of weights for this group
                                thread_sumw_group += weight;

                                // Store weights for this class and orientation
                                td_wsum_pdf_class          [iclass] += weight;
                                td_wsum_pdf_direction      [iclass] += weight;
                                td_wsum_prior_offsetx_class[iclass] += weight * (exp_old_offsetx[iimage] + exp_over_trans_x[itrans*nr_over_trans+iover_trans]);
                                td_wsum_prior_offsety_class[iclass] += weight * (exp_old_offsety[iimage] + exp_over_trans_y[itrans*nr_over_trans+iover_trans]);

                                double offsetx = mlModel.prior_offsetx_class.rptrAll()[iclass] - exp_old_offsetx.rptrAll()[iimage] - exp_over_trans_x[itrans*nr_over_trans+iover_trans];
                                double offsety = mlModel.prior_offsety_class.rptrAll()[iclass] - exp_old_offsety.rptrAll()[iimage] - exp_over_trans_y[itrans*nr_over_trans+iover_trans];

                                //this may cause some false share!
                                CHECK_NOT_NAN_OR_INF(thread_wsum_sigma2_offset += weight * (offsetx*offsetx+offsety*offsety));

                            } // end if weight >= exp_significant_weight
                        } // end loop iover_trans
                    } // end loop itrans

                    if (emit_test_output_prolific() && significant == 0)
                    #pragma omp critical
                    {
						static size_t count;
						if (count++ < 10) {
							std::cout << "updateOtherParams no weights exceed " << exp_significant_weight[iimage] << " in"
#define P(X) << " " << #X << ":" << X
								P(iimage) P(iclass) P(irot) P(iover_rot)
#undef P
								<< std::endl;
						}
                    }

                    if (emit_test_output_prolific())
                    #pragma omp critical
                    {
                        testos() << "~~~~TEST OUTPUT: thread_wsum_sigma2_offset "
                            << " iimage:" << iimage << " exp_significant_weight[iimage]:" << exp_significant_weight[iimage]
                            << " iclass:" << iclass
                            << " irot:" << irot
                            << " iover_rot:" << iover_rot
                            << " thread_wsum_sigma2_offset:" << thread_wsum_sigma2_offset
                            << ", was " << thread_wsum_sigma2_offset_before
                            << ".  Skipped " << skipped
                            << std::endl;
                    }
                    TUNING_SCOPE_ITER_END
                }// end loop iover_rot
            } // end loop irot
        }// end of iclass
    } // end loop iimage

#pragma omp atomic // TODO TODO TODO TODO
    mlModel.wsum_sumw_group.wptrAll()[0/*iimage*/] += thread_sumw_group;

#pragma omp critical
    {
        CHECK_NOT_NAN_OR_INF(mlModel.wsum_sigma2_offset += thread_wsum_sigma2_offset);

        if (emit_test_output_prolific()) {
            testos() << "~~~~TEST OUTPUT: wsum_sigma2_offset " << mlModel.wsum_sigma2_offset  << " from adding " << thread_wsum_sigma2_offset << std::endl;
        }
    }
    TUNING_SCOPE_PARA_END
}   // omp parallel if(...)

    for (int thread = 0; thread < maxthreads; thread++) {
        auto td = threadDataForEachThread[thread];
        if (td->written_wsum_pdf_prior_offset_max_weight()) {
            auto td_wsum_pdf_class           = td->wsum_pdf_class          ();
            auto td_wsum_pdf_direction       = td->wsum_pdf_direction      ();
            auto td_wsum_prior_offsetx_class = td->wsum_prior_offsetx_class();
            auto td_wsum_prior_offsety_class = td->wsum_prior_offsety_class();
            for (int n = 0; n < nr_classes; n++) {
                mlModel.wsum_pdf_class.wptrAll()[n]           += checkNonNegative(td_wsum_pdf_class          [n]);
                mlModel.wsum_pdf_direction[n]       +=                  td_wsum_pdf_direction      [n];
                mlModel.wsum_prior_offsetx_class.wptrAll()[n] +=                  td_wsum_prior_offsetx_class[n];
                mlModel.wsum_prior_offsety_class.wptrAll()[n] +=                  td_wsum_prior_offsety_class[n];
            }
        }
    }

    threadDataForEachThread.fini_wsum_pdf_prior_offset_max_weight();
        // these should not be used after here

    TUNING_SCOPE_PARA_BEGIN(updateOtherParams,nr_images)
	#pragma omp parallel for
    for (int iimage = 0; iimage < nr_images; iimage++)
    {   TUNING_SCOPE_ITER_BEGIN(updateOtherParams)
        //decode the real position,class->rot->trans->over_rot->over_trans
        //TODO : change the order to class->rot->over_rot->trans->over_trans maybe better
        int iclass, ipsi, iover_rot, itrans, iover_trans;
        DecodeFromCLASS(exp_metadata[iimage].CLASS, iclass, ipsi, iover_rot, itrans, iover_trans);
        int meta_class = iclass + 1;

        double A[3][3];
        double rot,tilt;
        double meta_psi = exp_over_rot_psi[ipsi*nr_over_rot+iover_rot];
#ifdef TODO
        Euler_angles2matrix(0, 0, meta_psi, A);
        Euler_matrix2angles(A, rot, tilt, meta_psi);
#endif
        double meta_xoff = exp_old_offsetx[iimage] + exp_over_trans_x[itrans*nr_over_trans+iover_trans];
        double meta_yoff = exp_old_offsety[iimage] + exp_over_trans_y[itrans*nr_over_trans+iover_trans];

        exp_metadata[iimage].ROT = 0;
        exp_metadata[iimage].TILT = 0;
        exp_metadata[iimage].PSI = meta_psi;
        exp_metadata[iimage].XOFF = meta_xoff;
        exp_metadata[iimage].YOFF = meta_yoff;
        exp_metadata[iimage].CLASS = meta_class;
        TUNING_SCOPE_ITER_END
    }
    TUNING_SCOPE_PARA_END

    ml_original_waypoint("updateOtherParams after");
}

// This is a fundamentally a simple algorithm
// It is
//      for a subset of {(class, irot, iover_rot, iimage, itrans, iover_trans)}
//          weight = function of the specific element of the subset
//          for (int n = 0; n < exp_current_Fsize2; n++)
//              exp_FrefsAndWeight.frefs_Rot_real()[class, irot, iover_rot, n] += Fimg_shift_nomask.real[iimage, itrans, iover_trans, n] * weight
//              exp_FrefsAndWeight.frefs_Rot_imag()[class, irot, iover_rot, n] += Fimg_shift_nomask.imag[iimage, itrans, iover_trans, n] * weight
//              exp_Fweight       [class, irot, iover_rot, n] += Mctf2_invsigma2                                    [n] * weight
//
// The following code goes parallel over the classes and images, because that is what a kernel object currently needs - it can currently only do one class and image at a time.
// It avoids having to lock the LHS for every add, by summing into accumulators, and only adding the accumulators into the LHS at the end.
//
// However measurements are showing that the Heap::allocZeroedDoubles is causing a performance problem
// which must be due to the sub-subset that share a class and image must be getting very small!
//
// The solution to this is to initially zero the accumulators, clear them as they are added into the original, and reuse the zeroed accumulators over multiple runs.
//
class BPK_ThreadAddTos {
public:
    BPK_ThreadAddTos() :
        addTo_exp_capacity(0),
        addTo_exp_Frefs_Rot_real(NULL),
        addTo_exp_Frefs_Rot_imag(NULL),
        addTo_exp_Fweight_Rot   (NULL) {}
    ~BPK_ThreadAddTos() { deallocate(); }
    void get(
        int needed,
        double* & addTo_exp_Frefs_Rot_real,
        double* & addTo_exp_Frefs_Rot_imag,
        double* & addTo_exp_Fweight_Rot) {
        allocate(needed);
        addTo_exp_Frefs_Rot_real = this->addTo_exp_Frefs_Rot_real;
        addTo_exp_Frefs_Rot_imag = this->addTo_exp_Frefs_Rot_imag;
        addTo_exp_Fweight_Rot    = this->addTo_exp_Fweight_Rot   ;
    }
private:
    // TODO CHANGE THESE TO FLOATS?
    int     addTo_exp_capacity;
    double* addTo_exp_Frefs_Rot_real;
    double* addTo_exp_Frefs_Rot_imag;
    double* addTo_exp_Fweight_Rot;

    void allocate(int needed) {
        if (needed <= addTo_exp_capacity) return;
        deallocate();
        needed += needed/2;     // grow exponentially
        TUNING_SCOPE_STEP(bpk_ThreadAddTos_allocate)
        addTo_exp_capacity = needed;
        addTo_exp_Frefs_Rot_real = Heap::allocZeroedDoubles(addTo_exp_capacity, __FILE__, __LINE__); // Zero here and after use
        addTo_exp_Frefs_Rot_imag = Heap::allocZeroedDoubles(addTo_exp_capacity, __FILE__, __LINE__); // to avoid having to zero before each use
        addTo_exp_Fweight_Rot    = Heap::allocZeroedDoubles(addTo_exp_capacity, __FILE__, __LINE__);
    }

private:
    void deallocate() {
        if (addTo_exp_capacity == 0) return;
        TUNING_SCOPE_STEP(bpk_ThreadAddTos_deallocate)
        Heap::freeDoubles(addTo_exp_Frefs_Rot_real); addTo_exp_Frefs_Rot_real = NULL;
        Heap::freeDoubles(addTo_exp_Frefs_Rot_imag); addTo_exp_Frefs_Rot_imag = NULL;
        Heap::freeDoubles(addTo_exp_Fweight_Rot   ); addTo_exp_Fweight_Rot    = NULL;
    }
};
std::vector<BPK_ThreadAddTos> bpk_threadAddTos_per_thread(omp_get_max_threads());


class BackProjectionKernel_replacement {
    const int exp_current_Fsize2;

    // TODO OBSOLETE - Undo - To keep the code simple, when updating the original, the copy points into the original
    //
    // The data for expImages.iclass_min() starts at &copy[0] and &original[copy_minus_original]
    //
    Exp_FrefsAndWeight saved_Exp_FrefsAndWeight;

public:
	static BackProjectionKernel_replacement* make(const int exp_current_Fsize2)
	{
#include "./util_heap_undefs.h"
		return sNewA(BackProjectionKernel_replacement, (exp_current_Fsize2));
#include "./util_heap_defs.h"
	}

	BackProjectionKernel_replacement(const int exp_current_Fsize2)
      : exp_current_Fsize2(exp_current_Fsize2) {
    }

    ~BackProjectionKernel_replacement() {
    }

    void __declspec(noinline) compute(const int exp_current_Fsize2) {
        //std::cout << "BackProjectionKernel_replacement::compute" << std::endl;

		auto const iclass_min	 = expImages.iclass_min(); 
		auto const iclass_max	 = expImages.iclass_max();
		auto const nr_rot		 = expCfg.nr_rot();
		auto const nr_over_rot	 = expCfg.nr_over_rot();
		auto const nr_images	 = expImages.nr_images(); 
		auto const nr_trans		 = expCfg.nr_trans();
		auto const nr_over_trans = expCfg.nr_over_trans();

        {
            ChooseImageAndClassStepForKNL chooseImageAndClassStepForKNL(expCfg,expImages,iter,__LINE__);
            ChooseImageAndClassStep		  chooseImageAndClassStep(expCfg,expImages,iter);
            const int iclassStep  = chooseImageAndClassStep.iclassStep();
            const int iimageStep  = chooseImageAndClassStep.iimageStep();
            const int irotStep    = chooseImageAndClassStep.irotStep();

            TUNING_SCOPE_STEP(backProjectionKernel)

#if defined(LOADBALANCEANALYZER)
            LoadBalanceAnalyzer lba("backProjectionKernel::compute", __FILE__, __LINE__, nr_images*(iclass_max+1-iclass_min));
#endif

            // Note: the innermost loop updates exp_Frefs_Rot[class,rotation]
            // but it does so purely as a sum reduction.
            //
            // However all those reductions need to be done in parallel
            // so a lock is used for each one
            //
            std::vector<omp_lock_t> rot_locks((chooseImageAndClassStep.numberOfClassesToDo() * nr_rot)*nr_over_rot);
            for (auto i = rot_locks.begin(); i < rot_locks.end(); i++) omp_init_lock(&*i);
            auto rot_locks_index = [&](int iclass, int irot, int iover_rot) {
                auto result = ((iclass-iclass_min) * nr_rot + irot)*nr_over_rot + iover_rot;
                assert(0 <= result && result < rot_locks.size());
                return result;
            };

            auto showSituation = [&](std::ostream & os) {
                os << "backProjectionKernel::compute "
					<< " nr Class steps:" << chooseImageAndClassStep.numberOfClassSteps()
                    << " nr Image steps:"  << (nr_images + iimageStep - 1)/iimageStep
                    << " nr rotations  :"  << nr_rot*nr_over_rot
                    << std::endl;
            };
            if (false) showSituation(std::cout);
            if (false) showSituation(std::cerr);

            TUNING_SCOPE_PARA_BEGIN(backProjectionKernel, chooseImageAndClassStepForKNL.numberOfIterations())

            #pragma omp parallel for collapse(3) schedule(dynamic)
            for (int iclassBegin = iclass_min; iclassBegin <= iclass_max; iclassBegin += iclassStep)    //50 or 100
            {
                for (int iimageBegin = 0; iimageBegin < nr_images; iimageBegin += iimageStep)
                {
					for (int irotBegin = 0; irotBegin < nr_rot; irotBegin += irotStep)
					{
						TUNING_SCOPE_ITER_BEGIN(backProjectionKernel)
#if defined(LOADBALANCEANALYZER)
						lba.iterationBegin();
#endif
						const int iclassEnd = std::min(iclass_max+1,iclassBegin + iclassStep);

						Map2dOptimizer_Kernel::RotationState rotationState(
							iclass_max+1 - iclass_min,
							nr_rot*nr_over_rot);

                        int const addTo_exp_len = (iclassEnd-iclassBegin)*exp_Frefs_Rot_len_per_class;
                        double* addTo_exp_Frefs_Rot_real;
                        double* addTo_exp_Frefs_Rot_imag;
                        double* addTo_exp_Fweight_Rot   ;
                        bpk_threadAddTos_per_thread[omp_get_thread_num()].get(addTo_exp_len,addTo_exp_Frefs_Rot_real,addTo_exp_Frefs_Rot_imag,addTo_exp_Fweight_Rot);
					    
                        // zeroed on creation and after each use
                        //  for (int i = 0; i < addTo_exp_len; i++) {
                        //      if (0 != (addTo_exp_Frefs_Rot_real[i] + addTo_exp_Frefs_Rot_imag[i] + addTo_exp_Fweight_Rot[i])) {
                        //          std::cerr << "Not zeroed" << std::endl; // BEVIN
                        //          EXIT_ABNORMALLY;
                        //      }
                        //  }
					    
                        auto manageAnAccumulator = [&](bool addingBack, int iclass, int irot, int iover_rot) {
					    
                            int const rotation          = irot * nr_over_rot + iover_rot;
                            int const iclass_rot        = iclass * nr_rot + irot;
                            int const iclass_rot_over   = iclass_rot * nr_over_rot + iover_rot;
					    
                            int const offsetTofirstSampleForRotationInAddTo =
                                (iclass - iclassBegin) * exp_Frefs_Rot_len_per_class +
                                rotation * exp_current_Fsize2;
					    
                            assert(0 <= offsetTofirstSampleForRotationInAddTo);
                            assert(offsetTofirstSampleForRotationInAddTo + exp_current_Fsize2 <= addTo_exp_len);
					    
                            auto arr = &addTo_exp_Frefs_Rot_real[offsetTofirstSampleForRotationInAddTo];
                            auto ari = &addTo_exp_Frefs_Rot_imag[offsetTofirstSampleForRotationInAddTo];
                            auto af  = &addTo_exp_Fweight_Rot   [offsetTofirstSampleForRotationInAddTo];
					    
                            if (!addingBack) {
                                SOAComplexDouble exp_Frefs_Rot_aux;
                                double* exp_Fweight_Rot_aux;
					    
                                // get the address where they will be when addingBack is true but don't change their state
                                exp_Frefs_Rot_aux.real = expIterData.FrefsAndWeight.frefs_Rot_real_noAccess (__LINE__,iclass_rot_over,exp_current_Fsize2);
                                exp_Frefs_Rot_aux.imag = expIterData.FrefsAndWeight.frefs_Rot_imag_noAccess (__LINE__,iclass_rot_over,exp_current_Fsize2);
                                exp_Fweight_Rot_aux    = expIterData.FrefsAndWeight.fweight_Rot_noAccess    (__LINE__,iclass_rot_over,exp_current_Fsize2);
					    
                                //  for (int n = 0; n < exp_current_Fsize2; n++) {
                                //      if (0 != (arr[n] + ari[n] + af[n])) {
                                //          std::cerr << "Not zeroed and not adding back" << std::endl; // BEVIN
                                //          EXIT_ABNORMALLY;
                                //      }
                                //  }
                                rotationState.setRotationState(iclass - iclassBegin, rotation,
                                    rotationState.RS_zeroedBuffer, exp_Frefs_Rot_aux.real, arr);
					    
                                return;
                            }
					    
                            auto lock = &rot_locks[rot_locks_index(iclass, irot, iover_rot)];
                            omp_set_lock(lock);
                                // Must not change from being zero to not being zero while we do this
					    
                            bool relevantFrefsAndWeightsAreZero = true
                                && expIterData.FrefsAndWeight.frefs_Rot_real_isZeroed(__LINE__,iclass_rot_over,exp_current_Fsize2)
                                && expIterData.FrefsAndWeight.frefs_Rot_imag_isZeroed(__LINE__,iclass_rot_over,exp_current_Fsize2)
                                && expIterData.FrefsAndWeight.fweight_Rot_isZeroed   (__LINE__,iclass_rot_over,exp_current_Fsize2);
					    
                            SOAComplexDouble exp_Frefs_Rot_aux;
                            double* exp_Fweight_Rot_aux;
                            if (relevantFrefsAndWeightsAreZero) {
                                // avoid having to zero them and then read the zeros and add to them!
                                exp_Frefs_Rot_aux.real = expIterData.FrefsAndWeight.frefs_Rot_real_write1st(__LINE__,iclass_rot_over,exp_current_Fsize2);
                                exp_Frefs_Rot_aux.imag = expIterData.FrefsAndWeight.frefs_Rot_imag_write1st(__LINE__,iclass_rot_over,exp_current_Fsize2);
                                exp_Fweight_Rot_aux    = expIterData.FrefsAndWeight.fweight_Rot_write1st   (__LINE__,iclass_rot_over,exp_current_Fsize2);
                            } else {
                                exp_Frefs_Rot_aux.real = expIterData.FrefsAndWeight.frefs_Rot_real_readWrite(__LINE__,iclass_rot_over,exp_current_Fsize2);
                                exp_Frefs_Rot_aux.imag = expIterData.FrefsAndWeight.frefs_Rot_imag_readWrite(__LINE__,iclass_rot_over,exp_current_Fsize2);
                                exp_Fweight_Rot_aux    = expIterData.FrefsAndWeight.fweight_Rot_readWrite   (__LINE__,iclass_rot_over,exp_current_Fsize2);
                            }
					    
                            auto rs = rotationState.getRotationState(iclass - iclassBegin, rotation);
                            if (rs == rotationState.RS_zeroedBuffer) {
                                //  for (int n = 0; n < exp_current_Fsize2; n++) {
                                //      if (0 != (arr[n] + ari[n] + af[n])) {
                                //          std::cerr << "Not zeroed and not zeroedBuffer" << std::endl;    // BEVIN
                                //          EXIT_ABNORMALLY;
                                //      }
                                //  }
                                if (relevantFrefsAndWeightsAreZero) {
                                    expIterData.FrefsAndWeight.zero(__LINE__,iclass_rot_over,exp_current_Fsize2);
                                    expIterData.FrefsAndWeight.zero(__LINE__,iclass_rot_over,exp_current_Fsize2);
                                    expIterData.FrefsAndWeight.zero(__LINE__,iclass_rot_over,exp_current_Fsize2);
                                }
                            } else {
                                assert(rs == rotationState.RS_addPending);
					    
                                TUNING_SCOPE_STEP(manageAnAccumulator_addPending)
					    
                                {
                                    // If this is a bottleneck, can change the caller to test the lock and do an available one
                                    if (relevantFrefsAndWeightsAreZero) {
                                        #pragma ivdep
                                        for (int n = 0; n < exp_current_Fsize2; n++) {
                                            checkFrefRotValue(exp_Frefs_Rot_aux.real[n] = arr[n]);
                                            checkFrefRotValue(exp_Frefs_Rot_aux.imag[n] = ari[n]);
                                                              exp_Fweight_Rot_aux   [n] = af [n];
                                            arr[n] = ari[n] = af[n] = 0.0;  // zero in preparation for next use
                                        }
                                    } else {
                                        #pragma ivdep
                                        for (int n = 0; n < exp_current_Fsize2; n++) {
                                            checkFrefRotValue(exp_Frefs_Rot_aux.real[n] += arr[n]);
                                            checkFrefRotValue(exp_Frefs_Rot_aux.imag[n] += ari[n]);
                                                              exp_Fweight_Rot_aux   [n] += af [n];
                                            arr[n] = ari[n] = af[n] = 0.0;  // zero in preparation for next use
                                        }
                                    }
                                }
                            }
                            rotationState.setRotationState(iclass - iclassBegin, rotation,
                                rotationState.RS_added, exp_Frefs_Rot_aux.real, arr);
					    
                            omp_unset_lock(lock);
                        };
					    
                        auto manageTheAccumulators = [&](bool addingBack) {
                            for (int iclass = iclassBegin; iclass < iclassEnd; iclass++) {
                                for (int irot = 0; irot < nr_rot; irot++) {
                                    for (int iover_rot = 0; iover_rot < nr_over_rot; iover_rot++) {
                                        manageAnAccumulator(addingBack, iclass, irot, iover_rot);
                                    }
                                }
                            }
                        };
					    
                        manageTheAccumulators(false);
					    
                        for (int iclass = iclassBegin; iclass < iclassEnd; iclass++)
                        {
                            TUNING_SCOPE_STEP(bpt_iclass)
					    
                            const int iimageEnd = std::min(nr_images, iimageBegin+iimageStep);
					    
                            for (int iimage = iimageBegin; iimage < iimageEnd; iimage++)
                            {
                                struct DeleteOnScopeExit {
                                    Map2dOptimizer_Kernel::BackProjection_Kernel* p;
                                    ~DeleteOnScopeExit() {
                                        fini();
                                    }
                                    void fini() {
                                        TUNING_SCOPE_STEP(bpk_freeKernel)
                                        p->releaseCapturedOutputs();
                                        p->release();
                                        p = NULL;
                                    }
                                } kernel;
					    
                                {
                                    TUNING_SCOPE_STEP(bpk_makeKernel)
                                    kernel.p =
                                        Map2dOptimizer_Kernel::BackProjection_Kernel::acquire(
											iter,
                                            rotationState,
                                            iclass - iclassBegin,
                                            nr_trans * nr_over_trans,
                                            exp_current_Fsize2);
                                }
					    
								int irotEnd = std::min(nr_rot, irotBegin + irotStep);
                                for (int irot = irotBegin; irot < irotEnd; irot++) // 36 / irotStep
                                {
                                    int iclass_rot = iclass * nr_rot + irot;
					    
                                    for (int iover_rot = 0; iover_rot < nr_over_rot; iover_rot++)	//2
                                    {
                                        if (!exp_Rot_significant[iclass][irot])
                                            continue;
					    
//                                      const int iclass_rot_over  = iclass_rot * nr_over_rot + iover_rot;
                                        const int        rot_over  =       irot * nr_over_rot + iover_rot;
                                        const int add_Frefs_offsetIndex = (((iclass-iclassBegin) * nr_rot + irot)*nr_over_rot + iover_rot);
                                        const int add_Frefs_offset      = add_Frefs_offsetIndex * exp_current_Fsize2;
					    
                                        // Inside the loop over all translations and all part_id sum all shift Fimg's and their weights
                                        // Then outside this loop do the actual backprojection
                                        /// Now that reference projection has been made loop over someParticles!
                                        for (int itrans = 0; itrans < nr_trans; itrans++) {
                                            for (int iover_trans = 0; iover_trans < nr_over_trans; iover_trans++) {

                                                auto Mctf2_invsigma2 = exp_local_Minvsigma2s[iimage].rptrAll();
					    
                                                int rot_trans_over = (rot_over*nr_trans + itrans)*nr_over_trans + iover_trans;
                                                double weight = exp_Mweight.get(iimage,iclass,rot_trans_over);
					    
                                                // Only deal with this sampling point if its weight was significant
                                                // Must be > to deal with 0 >= 0
                                                // NOTE : Must be >=
                                                if (weight >= exp_significant_weight[iimage])
                                                {
                                                    weight /= exp_sum_weight[iimage];
                                                    // Get the shifted image
                                                    int ishift = iimage * nr_over_trans * nr_trans + itrans * nr_over_trans + iover_trans;
					    
                                                    auto rotation    = irot*nr_over_rot + iover_rot;
                                                    auto translation = itrans * nr_over_trans + iover_trans;
					    
                                                    // TBD - BATCH THE TRANSFORMS EXPLOITING THE EXISTING BATCHING OF THE IMAGES
					    
                                                    kernel.p->append(
                                                        rotation,
                                                        translation,
                                                        addTo_exp_Frefs_Rot_real + add_Frefs_offset,
                                                        addTo_exp_Frefs_Rot_imag + add_Frefs_offset,
                                                        expIterData.Fimgs_shifted.efsRealConst(__LINE__, ishift,exp_current_Fsize2),
                                                        expIterData.Fimgs_shifted.efsImagConst(__LINE__, ishift,exp_current_Fsize2),
                                                        addTo_exp_Fweight_Rot    + add_Frefs_offset,
                                                        Mctf2_invsigma2,
                                                        weight);
                                                }
                                            }
                                        }
                                    }
                                }   // irot
                                {
                                    TUNING_SCOPE_STEP(bpk_compute)
                                    kernel.p->compute();
                                }
                            }
                        }
					    
                        manageTheAccumulators(true);

#if defined(LOADBALANCEANALYZER)
                        lba.iterationEnd();
#endif
                        TUNING_SCOPE_ITER_END
					}
                }
            }
            TUNING_SCOPE_PARA_END

            for (auto i = rot_locks.begin(); i < rot_locks.end(); i++) omp_destroy_lock(&*i);

        }   // LoadBalanceAnalyzer
    }   // compute()
};


void __declspec(noinline) backProjection()
{
	auto const iclass_min	 = expImages.iclass_min(); 
	auto const iclass_max	 = expImages.iclass_max();
	auto const nr_rot		 = expCfg.nr_rot();
	auto const nr_over_rot	 = expCfg.nr_over_rot();
	auto const nr_images	 = expImages.nr_images(); 
	auto const first_image	 = expImages.first_image(); 
	auto const last_image    = expImages.last_image(); 
	auto const nr_trans		 = expCfg.nr_trans();
	auto const nr_over_trans = expCfg.nr_over_trans();

    TUNING_SCOPE_STEP(backProjection)

    ml_original_waypoint("backProjection before");

    int exp_current_Fsize2 = expCfg.current_size()*(expCfg.current_size() / 2 + 1);

    // In doThreadPrecalculateShiftedImagesCtfsAndInvSigma2s() the origin of the exp_local_Minvsigma2s was omitted.
    // Set those back here
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        int igroup = exp_metadata[iimage].GROUP_NO-1;
        auto exp_local_Minvsigma2s_iimage = exp_local_Minvsigma2s[iimage].wptr(1);
        exp_local_Minvsigma2s_iimage[0] = CHECK_NOT_NAN_OR_INF( 1. / (sigma2_fudge * mlModel.sigma2_noise[igroup][0]) );
    }

	#pragma omp parallel
    {   TUNING_SCOPE_PARA_BEGIN(backProjection,nr_images)

        if (do_map == false || do_ctf_correction == false) {
            static bool laterTime = false;
            if (!laterTime) { std::cerr << "backProjection  do_map == false || do_ctf_correction == false" << std::endl; laterTime = true; }
			// BEVIN ANNOTATE THIS
            #pragma omp for
            for (int iimage = 0; iimage < nr_images; iimage++) {
                TUNING_SCOPE_ITER_BEGIN(backProjection)
                if (do_map == false) {
                    exp_local_Minvsigma2s[iimage].fill(1.0);
                }
                //Apply CTF to reference
                if (do_ctf_correction == false) {
                    exp_local_Fctfs_writable[iimage].fill(1.0);
                }
                TUNING_SCOPE_ITER_END
            }
        }

        // replace exp_local_Minvsigma2s with exp_local_Fctfs*exp_local_Minvsigma2s
		// BEVIN ANNOTATE THIS
        if (!exp_local_Fctfs_readonly.allOnes()) {
	        #pragma omp for
		    for (int iimage = 0; iimage < nr_images; iimage++) {
			    TUNING_SCOPE_ITER_BEGIN(backProjection)
				auto Mctf_invsigma2 = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
				auto ctfs        = exp_local_Fctfs_readonly[iimage].rptr(exp_current_Fsize2);
				auto minvsigma2s = exp_local_Minvsigma2s   [iimage].rptr(exp_current_Fsize2);
				assert(exp_local_Fctfs_readonly[iimage].size() == exp_current_Fsize2);
				assert(exp_local_Minvsigma2s   [iimage].size() == exp_current_Fsize2);
                for (int n = 0; n < exp_current_Fsize2; n++) {
                    Mctf_invsigma2[n] = ctfs[n] * minvsigma2s[n];
                }
				TUNING_SCOPE_ITER_END
            }
        }

        // multiply images by ctf and sigma2
        // This is one of the highest remaining sources of LLC misses
        // However it is very hard to improve locally
        //      Each exp_Fimgs_shifted is processed exactly once so comes from main memory
        //      Each Mctf_invsigma2 is loop invariant over the inner loops so comes from L1 or L2 and the loads and multiplies are completely hidden by the exp_Fimgs_shifted loads
        //
        // The only other possibility is to do this operation
        //      1 - when it is made in getShiftedImagesCtfsAndInvSigma2s
        //      2 - when they are used in BackProjectionKernel, updateOtherParams, other places?
        // This second might be a good idea because it is not clear to me that all the shifted images are used, so some may be being done unnecessarily
        //
        {   expIterData.Fimgs_shifted.capture_Mctf_invsigma2(nr_images,exp_current_Fsize2,exp_local_Minvsigma2s);
			// BEVIN ANNOTATE THIS
            #pragma omp for collapse(3)
            for (int iimage = 0; iimage < nr_images; iimage++)//nr_pool
            {
                for (int itrans = 0; itrans < nr_trans; itrans++)//100 ro 400
                {
                    for (int iover_trans = 0; iover_trans < nr_over_trans; iover_trans++)//4
                    {
                        TUNING_SCOPE_ITER_BEGIN(backProjection)

                        int ishift = iimage * nr_over_trans * nr_trans + itrans * nr_over_trans + iover_trans;

                        expIterData.Fimgs_shifted.doOrDeferMultipleByMctf_invsigma2(ishift,exp_current_Fsize2,iimage);
                            // This assumes that the Mctf_invsigma2 will survive until the undeferal

                        TUNING_SCOPE_ITER_END
                    }
                }
            }
        }

        // replace exp_local_Minvsigma2s with exp_local_Fctfs^2*exp_local_Minvsigma2s
		// BEVIN ANNOTATE THIS
        if (!exp_local_Fctfs_readonly.allOnes()) {
	        #pragma omp for
		    for (int iimage = 0; iimage < nr_images; iimage++) {
			    TUNING_SCOPE_ITER_BEGIN(backProjection)
				auto Mctf2_invsigma2 = exp_local_Minvsigma2s   [iimage].mptr(exp_current_Fsize2);
				auto ctfs            = exp_local_Fctfs_readonly[iimage].rptr(exp_current_Fsize2);
                for (int n = 0; n < exp_current_Fsize2; n++) {
                    Mctf2_invsigma2[n] *= ctfs[n];
                }
				TUNING_SCOPE_ITER_END
            }
        }

		// BEVIN ANNOTATE THIS
        #pragma omp for collapse(3)
        for (int iclass = iclass_min; iclass <= iclass_max; iclass++)//50 or 100
        {
            for (int irot = 0; irot < nr_rot; irot++)//36
            {
                for (int iover_rot = 0; iover_rot < nr_over_rot; iover_rot++)//2
                {
                    TUNING_SCOPE_ITER_BEGIN(backProjection)

                    int iclass_rot      = iclass * nr_rot + irot;
                    int iclass_rot_over = iclass_rot * nr_over_rot + iover_rot;
                    if (false && iter == 2 && iclass_rot_over == 1)
                    #pragma omp critical
                    {   // BEVIN
                        std::cout << "BackProjection zeroing iclass:" << iclass << " subsection:" << iclass_rot_over << std::endl;
                    }
                    expIterData.FrefsAndWeight.zero(__LINE__, iclass_rot_over, exp_current_Fsize2);

                    TUNING_SCOPE_ITER_END
                }
            }
        }

        TUNING_SCOPE_PARA_END
    }   // omp parallel

    ml_original_waypoint("backProjection middle1");

    {
        auto backProjectionKernel_replacement =
            BackProjectionKernel_replacement::make(exp_current_Fsize2);	// TODO eliminate the indirection

        double time_start = dtime();
        backProjectionKernel_replacement->compute(exp_current_Fsize2);
        double time_end = dtime();
        //NODE0ONLY std::cout<<" backProjectionKernel_replacement time : "<<(time_end-time_start)<<std::endl;

        sDelete(backProjectionKernel_replacement);
    }

    ml_original_waypoint("backProjection middle2");

	static volatile size_t count;
	count++;
	if (false && count > 0) {
		std::cerr << "NEW backProjection count:" 
			<< count 
			<< "iclass_min:"		<< iclass_min
			<< "iclass_max:"		<< iclass_max
			<< "wsum_data_real.size()"	<< mapModel.backprojector[0].pad_size*mapModel.backprojector[0].pad_Fsize /* wsum_data_real.size() */
			<< std::endl;
	}

    TUNING_SCOPE_PARA_BEGIN(backProjection,(iclass_max+1-iclass_min)*nr_rot*nr_over_rot)
    #pragma omp parallel for if (!emit_test_output())
    for (int iclass = iclass_min; iclass <= iclass_max; iclass++)//50 or 100
    {   TUNING_SCOPE_ITER_BEGIN(backProjection)

		std::ofstream* debug_of = nullptr;
		const char*    debug_fn = 
			// (iclass == 0) ? "C:/temp/map2d_backprojection_0_new.txt" :
			// (iclass == 7) ? "C:/temp/map2d_backprojection_7_new.txt" :
							nullptr;
		if (count >= 1 && debug_fn) 
			debug_of = ofstreamCheckingCreated::make(debug_fn);

        int tid = omp_get_thread_num();
        
        for (int ipsi = 0; ipsi < nr_rot; ipsi++)//36
        {
            for (int iover_rot = 0; iover_rot < nr_over_rot; iover_rot++)//2
            {
                int iorientclass_offset = iclass * nr_rot;
                int iorientclass = iorientclass_offset + ipsi;

                if (!exp_Rot_significant[iclass][ipsi])
                    continue;

                const int subsectionIndex = iorientclass*nr_over_rot + iover_rot;

                SOAComplexReadonly exp_Frefs_Rot_aux;
                exp_Frefs_Rot_aux.real            = expIterData.FrefsAndWeight.frefs_Rot_real_readonly(__LINE__, subsectionIndex, exp_current_Fsize2);
                exp_Frefs_Rot_aux.imag            = expIterData.FrefsAndWeight.frefs_Rot_imag_readonly(__LINE__, subsectionIndex, exp_current_Fsize2);
                const double* exp_Fweight_Rot_aux = expIterData.FrefsAndWeight.fweight_Rot_readonly   (__LINE__, subsectionIndex, exp_current_Fsize2);

                double A[3][3];
#ifdef TODO
                Euler_angles2matrix(0, 0, exp_over_rot_psi[ipsi*nr_over_rot + iover_rot], A);
                mapModel.set2DFourierTransform(tid, iclass, exp_Frefs_Rot_aux.real, exp_Frefs_Rot_aux.imag, exp_current_size, A, false, exp_Fweight_Rot_aux);
#endif
				if (debug_of) {
                    auto data_size = mapModel.backprojector[iclass].pad_size;
					auto dataLen = data_size*(data_size / 2 + 1);
					for (int i = 0; i < dataLen; i++) {
						*debug_of << "backProject data_real["<<i<<"]:" << mapModel.backprojector[iclass].data(0,0,i).real << std::endl;
						*debug_of << "backProject data_imag["<<i<<"]:" << mapModel.backprojector[iclass].data(0,0,i).imag << std::endl;
						*debug_of << "backProject weight   ["<<i<<"]:" << mapModel.backprojector[iclass].weight(0,0,i) << std::endl;
					}
				}

            }// end loop iover_rot
        } // end loop iorient

		auto print = [&](std::ostream& os, const char* prefix) {

            auto firstNonzero = [](double const* p, int l) { for (int i = 0; i < l; i++) if (p[i]!=0.0) return i; return 0; };
            auto firstNonzeroFloat = [](float const* p, int l) { for (int i = 0; i < l; i++) if (p[i]!=0.0) return i; return 0; };
            int wsum_data_size = mapModel.backprojector[iclass].pad_size;
			int wsum_data_Fsize  = wsum_data_size/2 + 1;
            int wsum_data_Fsize2 = wsum_data_size*wsum_data_Fsize;
            int wsum_weight_size = wsum_data_size*wsum_data_Fsize;

            auto exp_Frefs_Rot_aux_real = expIterData.FrefsAndWeight.frefs_Rot_real_readonly(__LINE__, 0, exp_current_Fsize2);
            //to exp_Frefs_Rot_aux_imag = expIterData.FrefsAndWeight.frefs_Rot_imag_readonly(__LINE__, 0, exp_current_Fsize2);
            auto exp_Fweight_Rot_aux	= expIterData.FrefsAndWeight.fweight_Rot_readonly   (__LINE__, 0, exp_current_Fsize2);

            auto wsum_data_aux = mapModel.backprojector[iclass].data.wptr();
            //to wsum_data_aux_imag = wsum_data_imag[iclass].rptrAll();
            auto wsum_weight_aux    = mapModel.backprojector[iclass].weight.wptr();

            auto ir = firstNonzero(exp_Frefs_Rot_aux_real, exp_current_Fsize2);
            auto id = firstNonzero((double*)wsum_data_aux,     wsum_data_Fsize2*2);
#if defined(FLOAT_PRECISION)
            auto iw = firstNonzeroFloat(wsum_weight_aux,        wsum_weight_size);
#else
            auto iw = firstNonzero(wsum_weight_aux,        wsum_weight_size);
#endif
            os << prefix << "backProjection"
#define P(X)   << " " << #X << ": " << X
               P(ir                        )
               P(exp_Frefs_Rot_aux_real[ir])
               P(id                        )
               P(wsum_data_aux[id].real    )
               P(iw                        )
               P(wsum_weight_aux[iw]       )
#undef P
               << std::endl;
		};

        if (false && iclass == 0) {	// BEVIN
			std::ofstream of("C:/temp/map2d_backrotate2D_new.txt");
			print(of,"");
			std::cerr << "Wrote C:/temp/map2d_backrotate2D_new.txt" << std::endl;
		}

        if (emit_test_output() && iclass == 0) {
			print(testos(),"~~~~TEST OUTPUT: ");
		}

		if (debug_of) {
			std::cerr << "Wrote " << debug_fn << std::endl;
			sDelete(debug_of);
		}

        TUNING_SCOPE_ITER_END
    } //end loop iclass
    TUNING_SCOPE_PARA_END

    //    static int flag = 0;
    //    if (iter == 2) {
    //        flag++;
    //    }
    //    if (flag >= 5) {
    //        EXIT_ABNORMALLY;
    //    }
    ml_original_waypoint("backProjection after");

}


void maximization()
{
    ml_original_waypoint("maximization before");

    NODE0ONLY std::cout << " Maximization ..." << std::endl;

	if (emit_test_output()) {
		for (int iclass = 0; iclass < nr_classes; iclass++) {
			testos() << "~~~~ maximization wsum_pdf_class[" << iclass << "] = " << mlModel.wsum_pdf_class.rptrAll()[iclass] << std::endl;
		}
	}

    // load balance by filtering those not being done
    std::vector<int> iclassToDo(0);
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        if (mlModel.wsum_pdf_class.rptrAll()[iclass] > 0.) iclassToDo.push_back(iclass);
		else mapModel.Irefs[iclass].zero();	// This is a memory-bound operation that should be optimized
    }
    // do only those that should be
	TUNING_SCOPE_PARA_BEGIN(maximization,iclassToDo.size())
    int iclassToDo_size = iclassToDo.size();
    #pragma omp parallel for if (!emit_test_output())
    for (int iclassIndex = 0; iclassIndex < iclassToDo_size; iclassIndex++)
    {
        const int iclass = iclassToDo[iclassIndex];

		if (emit_test_output()) {
            testos() << "~~~~TEST OUTPUT: maximization iclass:" << iclass << " wsum_pdf_class[iclass]:" << mlModel.wsum_pdf_class.rptrAll()[iclass] << std::endl;
        }

        if (mlModel.wsum_pdf_class.rptrAll()[iclass] <= 0.) {
            assert(false);// will not call
			mapModel.Irefs[iclass].zero();
        } else {
            auto sigma2_iclass           	= mlModel.sigma2_class[iclass].mptrAll();//update
            mlModel.sigma2_class[iclass].zero();
            auto tau2_iclass             	= mlModel.tau2_class[iclass].mptrAll();//update or not depend on update_tau2_with_fsc
            auto data_vs_prior_class_iclass = mlModel.data_vs_prior_class[iclass].mptrAll();//update or not depend on update_tau2_with_fsc
            mlModel.data_vs_prior_class[iclass].zero();
			//
            auto fsc_halves_class_iclass = mlModel.fsc_halves_class[iclass].wptr(ori_size/2+1);

            FDOUBLE maxWsumData   = 0;
            FDOUBLE maxWsumWeight = 0;
            int wsum_data_size = mapModel.backprojector[0].pad_size;
            for (int i = 0; i < wsum_data_size*(wsum_data_size/2+1); i++) {
                checkDoubleValue(mapModel.backprojector[iclass].data(0,0,i).real, "wsum_data_aux.real");
                checkDoubleValue(mapModel.backprojector[iclass].data(0,0,i).imag, "wsum_data_aux.imag");
                checkDoubleValue(mapModel.backprojector[iclass].weight(0,0,i), "wsum_weight_aux");
                maxWsumData   = std::max(maxWsumData,   std::abs(mapModel.backprojector[iclass].data(0,0,i).real));
                maxWsumData   = std::max(maxWsumData,   std::abs(mapModel.backprojector[iclass].data(0,0,i).imag));
                maxWsumWeight = std::max(maxWsumWeight, std::abs(mapModel.backprojector[iclass].weight(0,0,i)));
            }
            // reconstruct produces infinities unless it is given something to chew on
            assert(maxWsumData   > 0);
            assert(maxWsumWeight > 0);
            for (int i = 0; i < ori_size/2+1; i++) {
                checkDoubleValue(sigma2_iclass[i], "sigma2_aux");
                checkDoubleValue(data_vs_prior_class_iclass[i], "data_vs_prior_class_aux");
                checkDoubleValue(tau2_iclass[i], "tau2_aux");
            }
            checkDoubleValue(tau2_fudge_factor, "tau2_fudge_factor");

			static size_t count = 0;
			if (false && ++count >= 1) {	// BEVIN
				std::ofstream of("C:/temp/map2d_reconstruct_new.txt");
#define P(X)	<< #X << ":" << X << std::endl
				of P(iclass) P(ori_size) P(tau2_fudge_factor) P(do_map) P(gridding_nr_iter) P(mlModel.wsum_pdf_class.rptrAll()[iclass]) P(mapModel.minres_map);
#undef P
				auto pDouble = [&](const char* name, double* v, size_t s) {
					of << name;
					size_t count = 0;
					for (size_t i = 0; i < s; i++) {
						if (v[i] == 0.0) continue;
						of << " [" << i << "]:" << v[i];
						if (++count > 10) break;
					}
					of << std::endl;
				};
                auto pFloat = [&](const char* name, float* v, size_t s) {
                    of << name;
                    size_t count = 0;
                    for (size_t i = 0; i < s; i++) {
                        if (v[i] == 0.0) continue;
                        of << " [" << i << "]:" << v[i];
                        if (++count > 10) break;
                    }
                    of << std::endl;
                };
#if defined(FLOAT_PRECISION)
				pFloat("wsum_data_aux.real_image", (float*)mapModel.backprojector[iclass].data.wptr(), 2*wsum_data_size*(wsum_data_size/2+1));
				pFloat("wsum_weight_aux"   		, mapModel.backprojector[iclass].weight.wptr()   , wsum_data_size*(wsum_data_size/2+1));
				pFloat("sigma2_aux"             ,sigma2_iclass             ,ori_size/2+1);
				pFloat("data_vs_prior_class_aux",data_vs_prior_class_iclass,ori_size/2+1);
				pFloat("tau2_aux"               ,tau2_iclass               ,ori_size/2+1);
#else
                pDouble("wsum_data_aux.real_image", (double*)mapModel.backprojector[iclass].data.wptr(), 2*wsum_data_size*(wsum_data_size/2+1));
                pDouble("wsum_weight_aux"   		, mapModel.backprojector[iclass].weight.wptr()   , wsum_data_size*(wsum_data_size/2+1));
                pDouble("sigma2_aux"             ,sigma2_iclass             ,ori_size/2+1);
                pDouble("data_vs_prior_class_aux",data_vs_prior_class_iclass,ori_size/2+1);
                pDouble("tau2_aux"               ,tau2_iclass               ,ori_size/2+1);
#endif
				std::cerr << "Wrote C:/temp/map2d_reconstruct_new.txt" << std::endl;
			}

            mapModel.reconstruction(iclass, gridding_nr_iter, do_map, tau2_fudge_factor, tau2_iclass, sigma2_iclass,
                                    data_vs_prior_class_iclass, fsc_halves_class_iclass,
                                    mlModel.wsum_pdf_class.rptrAll()[iclass], false, false);

            if (false && emit_test_output() && iclass == 0) {   // BEVIN
                auto firstNonzeroDouble = [](double const* p, int l) { for (int i = 0; i < l; i++) if (p[i]!=0.0) return i; return 0; };
                auto firstNonzeroFloat = [](float const* p, int l) { for (int i = 0; i < l; i++) if (p[i]!=0.0) return i; return 0; };
                int wsum_data_Fsize = wsum_data_size/2 + 1;
                int wsum_data_Fsize2 = wsum_data_size*wsum_data_Fsize;
                int wsum_weight_size = wsum_data_size*wsum_data_Fsize;
                auto id = firstNonzeroDouble((double*)mapModel.backprojector[iclass].data.wptr(), 2*wsum_data_Fsize2);
#if defined(FLOAT_PRECISION)
                auto iw = firstNonzeroFloat(mapModel.backprojector[iclass].weight.wptr(), wsum_weight_size);
#else
                auto iw = firstNonzeroDouble(mapModel.backprojector[iclass].weight.wptr(), wsum_weight_size);
#endif
                testos() << "~~~~TEST OUTPUT: maximization"
#define P(X)        << " " << #X << ": " << X
                    P(iclass                    )
                    P(id                        )
                    P(mapModel.backprojector[iclass].data(0,0,id).real    )
                    P(iw                        )
                    P(mapModel.backprojector[iclass].weight(0,0,iw)       )
                    P(maxWsumData               )
                    P(maxWsumWeight             )
                    P(data_vs_prior_class_iclass[0])
#undef P
                    << std::endl;
            }

            for (int i = 0; i < ori_size*ori_size; i++) {
                checkIrefsValue(mapModel.Irefs[iclass](0,0,i));
            }
        }

    }
	TUNING_SCOPE_PARA_END

	ml_original_waypoint("maximization during");

    // Then perform the update of all other model parameters
    maximizationOtherParameters();

    // Keep track of changes in hidden variables
    // TOOD TODO use HiddenVarMonitor...
    // updateOverallChangesInHiddenVariables();

    ml_original_waypoint("maximization after");
}


void maximizationOtherParameters()
{
    // Calculate total sum of weights, and average CTF for each class (for SSNR estimation)
    double sum_weight = 0.;
    for (int iclass = 0; iclass < nr_classes; iclass++)
        sum_weight += mlModel.wsum_pdf_class.rptrAll()[iclass];//sum_weight += wsum_model.pdf_class[iclass];

    // TODO what is the right fix for this?
    // YongBei: if this is 0.0 menas all class is empty image
    // it is impossible in real case,so fix this may prevent us to detect bug...
    if (sum_weight == 0.0) sum_weight = 1.0;

    // Update average norm_correction
    if (do_norm_correction)
    {
        assert(sum_weight > 0.0);
        mlModel.avg_norm_correction = mlModel.wsum_avg_norm_correction / sum_weight;
    }


    // Update model.pdf_class vector (for each k)
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
        mlModel.pdf_class.wptrAll()[iclass] = mlModel.wsum_pdf_class.rptrAll()[iclass] / sum_weight;

        // for 2D also update priors of translations for each class!

        if (mlModel.wsum_pdf_class.rptrAll()[iclass] > 0.){
            mlModel.prior_offsetx_class.wptrAll()[iclass] = mlModel.wsum_prior_offsetx_class.rptrAll()[iclass] / mlModel.wsum_pdf_class.rptrAll()[iclass];
            mlModel.prior_offsety_class.wptrAll()[iclass] = mlModel.wsum_prior_offsety_class.rptrAll()[iclass] / mlModel.wsum_pdf_class.rptrAll()[iclass];
        }
        else{
            mlModel.prior_offsetx_class.wptrAll()[iclass] = 0.;
            mlModel.prior_offsety_class.wptrAll()[iclass] = 0.;
        }

		assert(mlModel.nr_directions==1);
        mlModel.pdf_direction[iclass][0] = mlModel.wsum_pdf_direction[iclass][0] / sum_weight;
        if (emit_test_output() && mlModel.pdf_direction[iclass][0] == 0.0) {
            static int count; static const int limit = 5;
            if (count++ < limit) {
                testos() << "~~~~ model_pdf_direction[" << iclass
                    << "] set to zero because  wsum_pdf_direction[iclass]:" <<  mlModel.wsum_pdf_direction[iclass][0]
                    << "/" << "sum_weight:" << sum_weight << std::endl;
            }
            if (count == limit) {
                testos() << "Suppressing future reports" << std::endl;
            }
        }

    }


    set_sigma2_offset(CHECK_NOT_NAN_OR_INF(mlModel.wsum_sigma2_offset / (2. * sum_weight)));

    if (emit_test_output()) {
        testos() << "~~~~TEST OUTPUT: sigma2_offset " << sigma2_offset << " from wsum_sigma2_offset:" << mlModel.wsum_sigma2_offset << "/(2. * sum_weight:" << sum_weight << ")" << std::endl;
    }

    // TODO: update estimates for sigma2_rot, sigma2_tilt and sigma2_psi!
    for (int igroup = 0; igroup < mlModel.nr_groups; igroup++)
    {
        auto sigma2_noise_igroup = mlModel.sigma2_noise[igroup].wptr(mlModel.ori_Fsize);
        auto wsum_sigma2_noise_igroup = mlModel.wsum_sigma2_noise[igroup].wptr(mlModel.ori_Fsize);
        // Factor 2 because of the 2-dimensionality of the complex-plane
        for (int n = 0; n < (ori_size/2+1); n++)
        {
            auto dividend = CHECK_NOT_NAN_OR_INF(wsum_sigma2_noise_igroup[n]);
            auto divisor  = (2. * mlModel.wsum_sumw_group.rptrAll()[igroup] * Npix_per_shell.rptrAll()[n]);

            // TODO what to do if this happens?
            if (divisor == 0.0) {
                std::cerr << "2. * wsum_sumw_group * Npix_per_shell[n] zero!"
                    << " wsum_sigma2_noise[n]:" << wsum_sigma2_noise_igroup[n]
                    << " wsum_sumw_group:" << mlModel.wsum_sumw_group.rptrAll()[igroup]
                    << " Npix_per_shell[n]:" << Npix_per_shell[n]
                    << std::endl;
                std::cout << "2. * wsum_sumw_group * Npix_per_shell[n] zero!"
                    << " wsum_sigma2_noise[n]:" << wsum_sigma2_noise_igroup[n]
                    << " wsum_sumw_group:" << mlModel.wsum_sumw_group.rptrAll()[igroup]
                    << " Npix_per_shell[n]:" << Npix_per_shell[n]
                    << std::endl;
                EXIT_ABNORMALLY;
            }

            sigma2_noise_igroup[n] = CHECK_NOT_NAN_OR_INF(dividend / divisor);

        }
    }


    // After the first iteration the references are always CTF-corrected
    refs_are_ctf_corrected = do_ctf_correction;

    // Some statistics to output

    mlModel.LL = mlModel.wsum_LL;

    mlModel.ave_Pmax = mlModel.wsum_ave_Pmax / sum_weight;
}


void writeClassesAndMetadata() {
    NODE0ONLY
    if (write_fn != "NULL") {
        std::string fn_class = write_path+write_fn;
        std::string fn_metadata = write_path+write_fn;
        // TOOD
		// ListOfImages listOfImages(model_Irefs);
        // MapOptimizer_base_new::writeClassesAndMetadata(fn_class, fn_metadata, metadata, listOfImages, testos());
    }
}


inline bool isSignificantAnyParticleAnyTranslation(int iclass, int irot)
{
    for (int iimage = 0; iimage < expImages.nr_images(); iimage++) {
        if (exp_Mcoarse_significant.isAnyTransSet(iimage, iclass, irot)) return true;
    }
    return false;
}

void calspace(int current_size,int set_nr_pool)
{

}

void debugStoreWeightedSums(){

    std::cout << " WARNING: norm_correction : "<< exp_metadata[0].NORM  << " for particle " << 0 <<std::endl;
    std::cout << " mymodel.current_size : " << current_size << " mymodel.ori_size= " << ori_size <<std::endl;
    std::cout << " coarse_size : " << coarse_size << std::endl;
    std::cout << " DIRECT_A2D_ELEM(exp_metadata2, my_image_no, exp_nr_imagas-1) : " <<exp_metadata[expImages.nr_images()-1].NORM << std::endl;
    std::cout << " mymodel.avg_norm_correction : " << mlModel.avg_norm_correction << std::endl;
    std::cout << " exp_wsum_norm_correction[ipart] : " << exp_wsum_norm_correction[0] << std::endl;
//    std::cout << " old_norm_correction : " << old_norm_correction << std::endl;
    std::cout << " wsum_model.avg_norm_correction : " << mlModel.wsum_avg_norm_correction << std::endl;
//    std::cout << " group_id : " << group_id << " mymodel.scale_correction[group_id] : " << mymodel.scale_correction[group_id] << std::endl;
    std::cout << " mymodel.sigma2_noise[group_id = 0] : " << mlModel.sigma2_noise[0][0] << std::endl;
    std::cout << " wsum_model.sigma2_noise[group_id = 0] : " << mlModel.wsum_sigma2_noise[0][0] << std::endl;
    std::cout << " exp_power_imgs[my_image_no = 0][0] : " << exp_power_imgs[0][0] << std::endl;
//    std::cout << " exp_wsum_scale_correction_XA[ipart] : " << exp_wsum_scale_correction_XA[ipart] << " exp_wsum_scale_correction_AA[ipart] : " << exp_wsum_scale_correction_AA[ipart] << std::endl;
//    std::cout << " wsum_model.wsum_signal_product_spectra[group_id] : " << wsum_model.wsum_signal_product_spectra[group_id] << " wsum_model.wsum_reference_power_spectra[group_id] : " << wsum_model.wsum_reference_power_spectra[group_id] << std::endl;
    std::cout << " exp_min_diff2[ipart = 0] : " << exp_min_diff2[0] << std::endl;
//    std::cout << " ml_model.scale_correction[group_id] : " << mymodel.scale_correction[group_id] << std::endl;
    std::cout << " exp_significant_weight[ipart = 0] : " << exp_significant_weight[0] << std::endl;
//    std::cout << " exp_max_weight[ipart = 0] : " << exp_max_weight[0] << std::endl;

    std::cerr << " part_id : " << 0 << std::endl;
    std::cerr << " ipart : " << 0 << std::endl;
    std::cerr << " exp_min_diff2[ipart = 0] : " << exp_min_diff2[0] << std::endl;
//    std::cerr << " logsigma2 : " << logsigma2 << std::endl;
    int group_id = 0;//mydata.getGroupId(part_id, 0);
    std::cerr << " group_id : " << group_id << std::endl;
//    std::cerr << " ml_model.scale_correction[group_id] : " << mymodel.scale_correction[group_id] << std::endl;
    std::cerr << " exp_significant_weight[ipart = 0] : " << exp_significant_weight[0] << std::endl;
//    std::cerr << " exp_max_weight[ipart = 0]= " << exp_max_weight[0] << std::endl;
    std::cerr << " ml_model.sigma2_noise[group_id = 0] : " << mlModel.sigma2_noise[0][0] << std::endl;
}

}   // namespace Map2dOptimizer_new


namespace Map2dOptimizer_new {
    static bool compareScalars() {
		// TODO
		return true;
    }


static void interestingCompare(const char* name) {
	if (0 && strcmp(name,"exp_local_Fctfs_writable") == 0 && (expCfg.current_size() > 0)) {
		std::cerr << "interestingCompare " << name << std::endl;
	}
}

    // get the original data offset
template<typename T>
size_t get_original_size(T& latest,const char* name){
	interestingCompare(name);
    if (	strcmp(name, "exp_local_Minvsigma2s"   ) == 0
		||	strcmp(name, "exp_local_Fctfs"		   ) == 0
		||	strcmp(name, "exp_local_Fctfs_writable") == 0
		||	strcmp(name, "exp_Frefs_Rot_real"	   ) == 0
		||	strcmp(name, "exp_Frefs_Rot_imag"	   ) == 0
		||	strcmp(name, "exp_Fweight_Rot"         ) == 0)
        return expCfg.current_size()*(expCfg.current_size()/2+1);
    else
        return latest.size();

}

template<typename T>
size_t get_original_offset_all(T& latest,const char* name){
	interestingCompare(name);
    if (	strcmp(name, "exp_local_Minvsigma2s"   ) == 0
		||	strcmp(name, "exp_local_Fctfs"		   ) == 0
		||	strcmp(name, "exp_local_Fctfs_writable") == 0)
        return expCfg.current_size()*(expCfg.current_size()/2+1);
    else
        return latest.size();
}

template<typename T>
size_t get_exp_nr(T& latest,const char* name){
	{
        std::cerr<<"check get_exp_nr"<<std::endl;
		EXIT_ABNORMALLY;
    }
	return 0;
}

    #define TEMPLATE_compareOrCopyArrayONESingleton(T1,T2)                                              \
    bool compareOrCopyArrayONE(T1 & latest, const T2* original, int & errorCount, const char* name) {   \
        auto original_size = get_original_size(latest,name);                                            \
		bool result = true;																				\
        for (int i = 0; i < original_size; i++) {                                                       \
            if (tandemCompareSingleton(latest[i], original[i], errorCount, name)) continue;             \
			result = false;																				\
			std::cerr << " " << name << "[" << i << "]" << std::endl;									\
			if (errorCount >= 4) break;																	\
		}                                                                                               \
		return result;																					\
    }                                                                                                   // end of macro

    #define TEMPLATE_compareOrCopyArrayONEVec(T1,T2)                                                    \
    bool compareOrCopyArrayONE(T1 & latest, const T2* original, int & errorCount, const char* name) {   \
        auto size = get_exp_nr(latest,name);                                                            \
		bool result = true;																				\
        auto original_size = get_original_size(latest,name);                                            \
        for (int i = 0; i < size; i++) {                                                                \
            if (compareOrCopyArrayONE(latest[i], original+i*original_size, errorCount, name)) continue; \
			result = false;																				\
			std::cerr << " " << name << "[" << i << "][?]" << std::endl;								\
			if (errorCount >= 4) break;																	\
        }                                                                                               \
		return result;																					\
    }                                                                                                   // end of macro

#define SEP
#define ELT(CLASSNAME, BASECLASS, BASEELT, SIZE_COMPUTABLE, LENGTH) \
    TEMPLATE_compareOrCopyArrayONESingleton(CLASSNAME,BASEELT)
    RUNTIME_SIZED_ARRAY_CLASSES_0
        // Generate these sizing expressions after all the static variables are available
#undef ELT
#define ELT(CLASSNAME, BASECLASS, BASEELT, SIZE_COMPUTABLE, LENGTH) \
    TEMPLATE_compareOrCopyArrayONEVec(CLASSNAME,BASEELT)
    RUNTIME_SIZED_ARRAY_CLASSES_1
        // Generate these sizing expressions after all the static variables are available
#undef ELT
#undef SEP

    #define TEMPLATE_compareOrCopyArrayVECofATOMIC(T1,T2)												\
    bool compareOrCopyArrayVEC(VectorOfSingle<T1> & latest, const T2* original, int & errorCount, const char* name) {      \
		interestingCompare(name);																		\
		bool result = true;																				\
        for (int i = 0; i < latest.size(); i++) {														\
            auto & latestElt   = latest[i];																\
            auto & originalElt = original[i];															\
            if (tandemCompareSingleton(latestElt, originalElt, errorCount, name)) continue;				\
			result = false;																				\
			std::cerr << " " << name << "[" << i << "]" << std::endl;									\
			if (errorCount >= 4) break;																	\
        }																								\
		return result;																					\
    }																									\
	bool compareOrCopyArrayVEC(std::vector<T1> & latest, const T2* original, int & errorCount, const char* name) {			\
		VectorOfSingle<T1> latestVoS(latest);															\
		return compareOrCopyArrayVEC(latestVoS, original, errorCount, name);							\
	}																									\
	// end of macro

    #define TEMPLATE_compareOrCopyArrayVECofVEC(T1,T2)													\
    bool compareOrCopyArrayVEC(std::vector<T1> & latest, const T2* original, int & errorCount, const char* name) {      \
		interestingCompare(name);																		\
		bool result = true;																				\
        for (int i = 0; i < latest.size(); i++) {														\
            auto & latestElt   = latest[i];																\
            int original_offset = get_original_offset_all(latestElt,name);								\
            auto   originalElt = original + i*original_offset;											\
            if (compareOrCopyArrayONE(latestElt, originalElt, errorCount, name)) continue;				\
			result = false;																				\
			std::cerr << " " << name << "[" << i << "][?]" << std::endl;								\
			if (errorCount >= 4) EXIT_ABNORMALLY;														\
        }																								\
		return result;																					\
    }																									\
	// end of macro

    TEMPLATE_compareOrCopyArrayVECofATOMIC(char,char)
    TEMPLATE_compareOrCopyArrayVECofATOMIC(int,int)
    TEMPLATE_compareOrCopyArrayVECofATOMIC(double,double)
    TEMPLATE_compareOrCopyArrayVECofATOMIC(float,float)
    TEMPLATE_compareOrCopyArrayVECofATOMIC(MetaDataElem,MetaDataElem)
#define SEP
#define ELT(CLASSNAME, BASECLASS, BASEELT, SIZE_COMPUTABLE, LENGTH) \
    TEMPLATE_compareOrCopyArrayVECofVEC(CLASSNAME,BASEELT)
    RUNTIME_SIZED_ARRAY_CLASSES
        // Generate these sizing expressions after all the static variables are available
#undef ELT
#undef SEP

    template <typename T1, typename T2>
    void compareOrCopyArrayPTR(T1 & latest, const T2* original, int & errorCount, const char* name) {
        std::cerr << "Don't know how to compare " << name << std::endl;
        assert(false);
    }
    static bool compareArrays() {
		//TODO
		return true;
    }
    bool waypointCompare() {
        bool result = true;
#if 0
        if (!compareScalars()) result &= false;
        if (!compareArrays ()) result &= false;
#endif
        return result;
    }

}   // namespace Map2dOptimizer_new
