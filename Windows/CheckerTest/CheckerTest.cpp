/***************************************************************************
*
* Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
*		Dana-Farber Cancer Institute, Harvard Medical School and Peking University
*		"Brett, Bevin" Intel Corporation
*		"Brett, Bevin" After retiring from Intel Corporation
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

#include "stdafx.h"
#include "../../src/checker.h"
#include "../../src/primitives_for_each_os.h"

// Utilities
//
static size_t fails;
static std::ostream& failed() { fails++; return std::cerr; }

// Example of use

// We can declare the program structure in a single place
// and verify that the call patterns match this expection
#define CHECKER_SCOPES \
	ELT(Main     , Main     , Scope_Flags::None())										SEP \
	ELT(OuterLoop, OuterLoop, Scope_Flags::Loop())										SEP \
	ELT(InnerLoop, InnerLoop, Scope_Flags::Loop())										SEP	\
	ELT(InnerIter, InnerIter, Scope_Flags::Iter())											\
	// end of macro

#define ELT(N,W,V)  extern Scope scope_##N;
#define SEP
CHECKER_SCOPES
#undef SEP
#undef ELT

#define ELT(N,W,V)  Scope scope_##N(#N, V);
#define SEP 
CHECKER_SCOPES
#undef SEP
#undef ELT


// Testing
//
static void test01() {
	using namespace Checker;

	std::cout << "CheckerTest::test01 verifies that cells can be added to, saved, loaded, and read from a new csv" 
		<< std::endl;

	// Use an in-memory file rather than disk-based
	// 
	Instance instanceO;
	{
		auto r0 = instanceO.addRowIfMissing("r0");
		auto r1 = instanceO.addRowIfMissing("r1");
		auto c0 = instanceO.addColIfMissing("c0");
		auto c1 = instanceO.addColIfMissing("c1");
		instanceO.set(r0, c0, "3:4");
		instanceO.set(r1, c1, DoublePair(5, 6));
	}

	std::ostringstream ostringstream;
	instanceO.write(ostringstream);

	std::istringstream istringstream(ostringstream.str());

	Instance instanceI;
	instanceI.read(istringstream);

	{
		auto r0 = instanceI.addRowIfMissing("r0");
		auto r1 = instanceI.addRowIfMissing("r1");
		auto c0 = instanceI.addColIfMissing("c0");
		auto c1 = instanceI.addColIfMissing("c1");
		auto x34 = instanceI.getDbl(r0, c0);
		auto x56 = instanceI.getStr(r1, c1);
		if (x34.lo != 3 || x34.hi != 4)			 failed() << "not 3 4" << std::endl;
		if (DoublePair(5, 6) != DoublePair(x56)) failed() << "not 5 6" << std::endl;
	}

}


static void test02() {
	using namespace Checker;

	std::cout << "CheckerTest::test02 verifies that scope timings can be added to am existing csv"
		<< " and compared across runs"
		<< std::endl;

	auto checker_ftrace_fnm = currentDirectory() + "/CheckerTest_test02_ftrace_temp";	// _#.txt appended by Checker
	Checker::setFtraceFile(checker_ftrace_fnm);

	auto checker_i_fnm = currentDirectory() + "/tmp/CheckerTest02_testi.csv";
	auto checker_o_fnm = currentDirectory() + "/tmp/CheckerTest02_testo.csv";
	std::cout << "Checker::unitTest test02"
		<< " reads and writes " << checker_i_fnm << " and " << checker_o_fnm
		<< " and also " << checker_ftrace_fnm
		<< std::endl;

	double const aFewSecs = 0.1;
	Checker::configureFTraceSkipping(
		aFewSecs, 
		aFewSecs*20);

	for (size_t alternative = 0; alternative < 3; alternative++) {

		Benchmark benchmark;

		switch (alternative) {
		case 0: benchmark.setFiles("",			  checker_i_fnm); break;
		case 1: benchmark.setFiles(checker_i_fnm, checker_o_fnm); break;
		case 2: benchmark.setFiles(checker_o_fnm, ""           ); break;
		default: assert(false);
		}

		benchmark.setDefaultColumnHeader(
			benchmark.addColIfMissing("c_temp"));

		Life life_Main(scope_Main.init(__FILE__, __LINE__), nullptr);

		Life life_OuterLoop(scope_OuterLoop.init(__FILE__, __LINE__), &life_Main);
		for (int iter = 0; iter < 2; iter++) {	// This mimics the iteration of the ROME3D
			auto scope_OuterIter_name = "OuterIter_" + std::to_string(iter);
			Scope scope_OuterIter(scope_OuterIter_name.c_str(), Scope_Flags::Loop());
			Life life_OuterIter(scope_OuterIter.init(__FILE__, __LINE__), &life_OuterLoop, &benchmark);
			//...
			Life life_InnerLoop(scope_InnerLoop.init(__FILE__, __LINE__), &life_OuterIter);
			for (int i = 0; i < 10; i++) {		// This mimics some per-image loop of the ROME3D
				Life life_InnerIter(scope_InnerIter.init(__FILE__, __LINE__), &life_InnerLoop);
				AccurateTimer spinning;
				while (spinning.sinceInited() / spinning.countPerMicrosecond() < aFewSecs*0.3e6);
			}
		}
	}

	Checker::unsetFtraceFile();
}


static void test03() {
	using namespace Checker;

	std::cout << "CheckerTest::test03 verifies that load balance can be checked"
		<< std::endl;

	static const int iters = 144;
	static volatile int sum[iters];

	CHECKER_PARALLEL_COUNTED_FOR1(int, i, 0, iters) {
		sum[i] = 0;
		for (int j = 0; j < i; j++) {
			for (int k = 0; k < 1000; k++) {
				sum[i] += i*j*k;
			}
		}
	} CHECKER_PARALLEL_FOR_END
}


int main()
{
	std::cout << "CheckerTest running in " << currentDirectory() << std::endl;

	if (0) test01();
	if (0) test02();
	if (1) test03();

	if (fails) std::cerr << "FAILED" << std::endl;

	return fails ? 1 : 0;
}

