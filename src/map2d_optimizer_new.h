/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * This program is free software) you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation) either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY) without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#pragma once

#include "./map_optimizer_base_new.h"
#include "./array_vector.h"
#include <typeinfo>

namespace Map2dOptimizer_new {
    
    using namespace MapOptimizer_base_new;
    
#define exp_nr_rot exp_nr_rot_NOT_USED
#define exp_nr_psi exp_nr_psi_NOT_USED
    
#define SEP
#define ELT(T, N) extern T N;
    MAP2DOPTIMIZER_STATIC_SCALARS
#undef ELT
#undef SEP
    
#undef exp_nr_rot
#undef exp_nr_psi
    
	void set_iter(int to);

	// major phases
	//
    void setupMap2dOptimizer();
	void readImages();
	void prepare();
	void iterate();
    void writeClassesAndMetadata();
    void writeClassesAndMetadata(std::string filename_mrcs,std::string filename_star);
    void destroyMap2dOptimizer();

	// Support for tandem execution, used for debugging new versions of algorithms
	//
	TandemExecution::AlgorithmP prepareAlgorithm(TandemExecution::AlgorithmsP algorithms);
	TandemExecution::AlgorithmP iterateAlgorithm(TandemExecution::AlgorithmsP algorithms);

	bool waypointCompare();

}	// namespace Map2dOptimizer_new
