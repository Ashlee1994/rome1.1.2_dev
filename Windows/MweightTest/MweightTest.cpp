#include "../../src/exp_mweight.h"
#include "../../src/exp_mweight_old.h"


static const int numberOfRGs = 1024;	// BIG HACK
std::default_random_engine				random_generator[numberOfRGs];
std::uniform_real_distribution<float>	rand01f_distribution(0.0, 1.0);
auto rand01f = []()->float { 
	return rand01f_distribution(random_generator[omp_get_thread_num()]);
};

#define INNER_DIMENSIONS \
	ELT(idir,			exp_nr_dir,		   200,    800)		SEP \
	ELT(ipsi,			exp_nr_psi,			24,		48)		SEP \
	ELT(itrans,			exp_nr_trans,	   100,    400)		SEP \
	ELT(iover_rot,		exp_nr_over_rot,	 1, (1<<9))		SEP \
	ELT(iover_trans,	exp_nr_over_trans,   1, (1<<6))			\
	// end of macro 
#define DIMENSIONS \
	ELT(iimage,			exp_nr_images,		10,	100000)		SEP \
	ELT(iclass,			exp_nr_classes,		 1,	   100)		SEP \
	INNER_DIMENSIONS											\
	// end of macro 

#define SEP
#define ELT(N,L,LO,HI) int L = LO;
DIMENSIONS
#undef ELT
#undef SEP

void init(Exp_Mweight_new& m,
	int nr_images,
	int nr_classes,
	int nr_dir,
	int nr_psi,
	int nr_over_rot,
	int nr_trans,
	int nr_over_trans) {
	m.init(nr_images, nr_classes, nr_dir, nr_psi, nr_over_rot, nr_trans, nr_over_trans);
}
void init(Exp_Mweight_old& m,
	int nr_images,
	int nr_classes,
	int nr_dir,
	int nr_psi,
	int nr_over_rot,
	int nr_trans,
	int nr_over_trans) {
	m.init(nr_images, nr_classes, nr_dir, nr_psi, nr_over_rot*nr_trans*nr_over_trans);
}

void clear(Exp_Mweight_new & m, int iimage) { m.clear      (iimage); }
void clear(Exp_Mweight_old & m, int iimage) { m.clear_image(iimage); }

typedef std::vector< std::pair<int, double> > MSS_old;
typedef MweightsForSomeSpinsAndSlides		  MSS_new;

MSS_new & mweightsForSomeSpinsAndSlides(Exp_Mweight_new& m, int iimage, int iclass, int idir, int ipsi) {
	return m.mweightsForSomeSpinsAndSlides(iimage, iclass, idir, ipsi);
}
MSS_old & mweightsForSomeSpinsAndSlides(Exp_Mweight_old& m, int iimage, int iclass, int idir, int ipsi) {
	return m.wptr_sparse(iimage, iclass, idir, ipsi);
}

void insert(MSS_new& m, int iover_rot, int itrans, int iover_trans, double value) {
	m.insert(iover_rot, itrans, iover_trans, value);
}
void insert(MSS_old& m, int iover_rot, int itrans, int iover_trans, double value) {
	int ihidden = itrans*exp_nr_over_rot + iover_rot;
	ihidden = ihidden*exp_nr_over_trans + iover_trans;
	m.push_back(std::pair<int, double>(ihidden, value));
}

struct MSS_new_iterator {
	MSS_new* parent;
	size_t   i;
	MSS_new_iterator() : parent(nullptr), i(0) {}
	MSS_new_iterator(MSS_new* parent, size_t i) : parent(parent), i(i) {}
	MSS_new_iterator operator++(int) { return MSS_new_iterator(parent, i++); }
	bool operator!=(MSS_new_iterator& rhs) const { return i != rhs.i; }
};

inline MSS_new_iterator begin(MSS_new& m) { return MSS_new_iterator(&m, 0);        }
inline MSS_new_iterator end  (MSS_new& m) { return MSS_new_iterator(&m, m.size()); }

typedef MSS_old::iterator MSS_old_iterator;
inline MSS_old_iterator begin(MSS_old& m) { return m.begin(); }
inline MSS_old_iterator end  (MSS_old& m) { return m.end();   }

inline double value(MSS_new_iterator & i) { return i.parent->value(i.i); }
inline double value(MSS_old_iterator & i) { return i->second;            }

int numberOfHotspots = 10;		// the number of hotspots when aligning a specific class and image and rotation and translation
int hotspotRadius    =  4;		// how near a hotspot to become significant

template<class Exp_Mweight>
// typedef Exp_Mweight_new Exp_Mweight;
// typedef Exp_Mweight_old Exp_Mweight;
void test(const char* name) {
	std::cout << "Test start " << name << std::endl;
	Microseconds elapsed;
	{
		Microseconds start = timeInMicroseconds();
		Exp_Mweight exp_Mweight_coarse;
		Exp_Mweight exp_Mweight_fine;

		init(exp_Mweight_coarse,
			exp_nr_images,
			exp_nr_classes,
			exp_nr_dir,
			exp_nr_psi,
			1,
			exp_nr_trans,
			1);

		init(exp_Mweight_fine,
			exp_nr_images,
			exp_nr_classes,
			exp_nr_dir,
			exp_nr_psi,
			exp_nr_over_rot,
			exp_nr_trans,
			exp_nr_over_trans);

		// getAllSquaredDifferencesCoarse 
		std::atomic<int> hotspots = 0;
		#pragma omp parallel for collapse(3) schedule(dynamic)
		for (int iimage = 0; iimage < exp_nr_images; iimage++) {
			for (int iclass = 0; iclass < exp_nr_classes; iclass++) {
				for (int idir = 0; idir < exp_nr_dir; idir++) {
					for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) {
						auto& exp_Mweight_coarse_sub = mweightsForSomeSpinsAndSlides(exp_Mweight_coarse,iimage, iclass, idir, ipsi);
						// hack to drastically reduce the number of random numbers needed
						int remainingHotspots = numberOfHotspots;
						int remainingUntilNextHotspot = rand01f() * 2.0 * exp_nr_trans / remainingHotspots;	
						for (int itrans = 0; itrans < exp_nr_trans; itrans++) {
							int const iover_rot   = 0;
							int const iover_trans = 0;
							if (remainingUntilNextHotspot-- > 0) continue;
							insert(exp_Mweight_coarse_sub,
								iover_rot, itrans, iover_trans,
								1.0);
							hotspots++;
							remainingHotspots--;
							if (!remainingHotspots) break;
							remainingUntilNextHotspot = rand01f() * 2.0 * (exp_nr_trans - itrans) / remainingHotspots;
						}
					}
				}
			}
		}
		std::cerr << "Hotspots inserted:" << hotspots << std::endl;

		// convertSquaredDifferencesToWeights
		std::atomic<int> max_hotspots_found = 0;
		#pragma omp parallel for collapse(1)
		for (int iimage = 0; iimage < exp_nr_images; iimage++) {
			bool found(false); int found_iimage, found_iclass, found_idir, found_ipsi; double found_v;
			for (int iclass = 0; iclass < exp_nr_classes; iclass++) {
				for (int idir = 0; idir < exp_nr_dir; idir++) {
					for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) {
						auto& someWeights = mweightsForSomeSpinsAndSlides(exp_Mweight_coarse, iimage, iclass, idir, ipsi);
						for (auto i = begin(someWeights); i != end(someWeights); i++) {
							found = true;
							found_iimage = iimage; found_iclass = iclass; found_idir = idir; found_ipsi = ipsi;
							found_v = value(i);
						}
					}
				}
			}
			clear(exp_Mweight_coarse,iimage);
			if (found) {
				auto& mss = mweightsForSomeSpinsAndSlides(exp_Mweight_coarse, found_iimage, found_iclass, found_idir, found_ipsi);
				insert(mss, 0, 0, 0, found_v);
				max_hotspots_found++;
			}
		}
		std::cerr << "Max hotspots inserted:" << max_hotspots_found << std::endl;

		// getAllSquaredDifferencesFine 
		std::atomic<int> near_hotspots = 0;
		#pragma omp parallel for collapse(3) schedule(dynamic)
		for (int iimage = 0; iimage < exp_nr_images; iimage++) {
			for (int iclass = 0; iclass < exp_nr_classes; iclass++) {
				for (int idir = 0; idir < exp_nr_dir; idir++) {
					for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) {
						auto& someCoarseWeights = mweightsForSomeSpinsAndSlides(exp_Mweight_coarse, iimage, iclass, idir, ipsi);
						auto  someFineWeights = &someCoarseWeights;	someFineWeights = nullptr; // hack to get the type
						if (someCoarseWeights.size() == 0) continue;
						for (int itrans = 0; itrans < exp_nr_trans; itrans++) {
							for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++) {
								for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++) {
									if (itrans + iover_trans > hotspotRadius || iover_rot > hotspotRadius) continue;
									if (!someFineWeights) someFineWeights = &mweightsForSomeSpinsAndSlides(exp_Mweight_fine, iimage, iclass, idir, ipsi);
									insert(*someFineWeights, iover_rot, itrans, iover_trans, 1.0);
									near_hotspots++;
								}
							}
						}
					}
				}
			}
		}
		std::cerr << "Near hotspots inserted:" << near_hotspots << std::endl;

		elapsed = timeInMicroseconds() - start;
	}
	std::cout << "Test end, took " << long(elapsed) << " microsecs" << std::endl;
}

int main()
{
	// The algorithm is currently basically linear in these dimensions
	// Improving this would require multiplexing over the dir and psi 
	// - which may be possible given the way they are all cleared and then just the max pushed back in and then pulled out again later!
	//
	exp_nr_images		= 10;
	exp_nr_classes		= 10;
	exp_nr_dir			= 30;
	exp_nr_psi			= 36;

	// The algorithm should be basically constant overhead on these and linear in usage
	//
	exp_nr_trans		= 400;
	exp_nr_over_rot		= 1<<9;
	exp_nr_over_trans	= 1<<6;

	//test("?");
	test<Exp_Mweight_old>("old");
	test<Exp_Mweight_new>("new");
	return 0;
}

