#include "../../src/peakFinder.h"

static int test()
{
	PeakFinder::Bounds	peakFinderBounds;
	auto pushBound = [&](const char* name, int bound) {
		auto lo = (bound + 1) / 2;
		auto hi = 2 * bound + 1;
		peakFinderBounds.push_back({ lo,  hi });
		std::cerr << "Bound " << name << " [" << lo << ".." << hi - 1 << "]" << std::endl;
		std::cout << "Bound " << name << " [" << lo << ".." << hi - 1 << "]" << std::endl;
	};
	pushBound("iimage_tile",     8);
	pushBound("iimage_sub_tile", 7);
	pushBound("itrans_sub_tile", 7);

	PeakFinder peakFinder(peakFinderBounds);
	// TODO learns across iterations also

	{
		PeakFinder peakFinder(peakFinderBounds);

		#pragma omp parallel for
		for (int tries = 0; tries < 10000; tries++) {
			if (peakFinder.hasFoundPeak()) continue;
			class A : public NoCopy, public PeakFinder::Assignment {
			public:
				virtual void coords(PeakFinder::ResultIndex resultIndex, PeakFinder::Coordinates const & coordinates) {
					auto x = coordinates[0]; 
					auto y = coordinates[1]; 
					auto z = coordinates[2];
					result = -float(square(x - 8) + square(y - 7) + square(z - 7));
					// std::cout << "f(x:" << x << ", y:" << y << ") = " << result << std::endl;
				}
				float result;
			} a;
			auto ri = peakFinder.askForAssignment(a);
			if (ri == PeakFinder::noAssignment) continue;
			peakFinder.reportResult(ri, a.result);
			// std::cout << a.result << ", bestResult:" << peakFinder.bestResult() << " after " << peakFinder.numberAssignmentsReported() << " trials" << std::endl;
		}

		auto msg = [&](std::ostream& os) {
			os << "bestResult:" << peakFinder.bestResult();
			auto & bestCoords = peakFinder.bestCoordinates();
			const char* sep = " found at ";
			for (auto c : bestCoords) {
				os << sep << c;
				sep = ",";
			}
			os << " in a volume " << peakFinder.volume() << " searched with " << peakFinder.numberAssignmentsReported() << " trials" << std::endl;
		};
		msg(std::cout);

	}

    return 0;
}

int main() {
	auto result = test();
	std::cerr << "test returned " << result << std::endl;
	return result;
}
