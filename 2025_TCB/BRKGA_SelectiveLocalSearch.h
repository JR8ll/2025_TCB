#pragma once
#include <BRKGA.h>
#include <type_traits>

class GaDecoderJobListSched;

 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //+++++ BEWARE: This class definition requires the following modifications to the original brkgaAPI:			  +++++
 //+++++ 1) BRKGA.h -> line 123 change access modificator for attributes of class BRKGA from private to protected +++++
 //+++++ 2) Population.h -> line 25 in class Population add forward definition for friend class					  +++++
 //+++++                    "template< class Decoder, class RNG > friend class BRKGA_SelectiveLocalSearch;		  +++++
 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<class Decoder, class RNG>
class BRKGA_SelectiveLocalSearch : public BRKGA<Decoder, RNG> {
	
private: 
	static_assert(std::is_same<Decoder, GaDecoderJobListSched>::value,
		"BRKGA_SelectiveLocalSearch requires GaDecoderJobListSched");
	int nLocalSearch;

public:
	BRKGA_SelectiveLocalSearch(int n, int nPop, double pElt, double pRpM,
		double rhoe, Decoder& decoder, RNG& rng,
		int K, int MAXT)
		: BRKGA<Decoder, RNG>(n, nPop, pElt, pRpM, rhoe, decoder, rng, K, MAXT)	{
	}

	void setLocalSearchFraction(double fraction) {
		this->nLocalSearch = floor(fraction * (double)pe);
	}

	void evolve(unsigned generations = 1) {
		if (generations == 0) { throw std::range_error("Cannot evolve for 0 generations."); }
		for (unsigned i = 0; i < generations; ++i) {
			for (unsigned j = 0; j < K; ++j) {
				this->evolution(*this->current[j], *this->previous[j]);	// First evolve the population (curr, next)
				std::swap(this->current[j], this->previous[j]);		// Update (prev = curr; curr = prev == next)
			}
		}
	}

	void applyNonPersistantLocalSearch(double fraction) { // apply local search to the best <fraction> percent of the population
		if (fraction <= 0.0 || p == 0) {
			return;
		}
		if (fraction > 1.0) fraction = 1.0;

		size_t num_top = min(static_cast<size_t>(fraction * p), p);
		vector<unsigned> indices(p);
		iota(indices.begin(), indices.end(), 0);

		partial_sort(indices.begin(), indices.begin() + num_top, indices.end(), [this](unsigned a, unsigned b) {
			return current[a]->getBestFitness() < current[b]->getBestFitness();
			});

		for (size_t i = 0; i < num_top; ++i) {
			auto* gaJobListSchedDecoder = dynamic_cast<GaDecoderJobListSched*>(&current[indices[i]]);
			gaJobListSchedDecoder->applyNonPersitentLocalSearch();
		}
	}	

	//template< class Decoder, class RNG >
	inline void evolution(Population& curr, Population& next) {
		// We now will set every chromosome of 'current', iterating with 'i':
		unsigned i = 0;	// Iterate chromosome by chromosome
		unsigned j = 0;	// Iterate allele by allele

		// 2. The 'pe' best chromosomes are maintained, so we just copy these into 'current':
		while (i < pe) {
			for (j = 0; j < n; ++j) { next(i, j) = curr(curr.fitness[i].second, j); }

			next.fitness[i].first = curr.fitness[i].first;
			next.fitness[i].second = i;
			++i;
		}

		// 3. We'll mate 'p - pe - pm' pairs; initially, i = pe, so we need to iterate until i < p - pm:
		while (i < p - pm) {
			// Select an elite parent:
			const unsigned eliteParent = (refRNG.randInt(pe - 1));

			// Select a non-elite parent:
			const unsigned noneliteParent = pe + (refRNG.randInt(p - pe - 1));

			// Mate:
			for (j = 0; j < n; ++j) {
				const unsigned sourceParent = ((refRNG.rand() < rhoe) ? eliteParent : noneliteParent);

				next(i, j) = curr(curr.fitness[sourceParent].second, j);

				//next(i, j) = (refRNG.rand() < rhoe) ? curr(curr.fitness[eliteParent].second, j) :
				//		                              curr(curr.fitness[noneliteParent].second, j);
			}

			++i;
		}

		// We'll introduce 'pm' mutants:
		while (i < p) {
			for (j = 0; j < n; ++j) { next(i, j) = refRNG.rand(); }
			++i;
		}

		// Time to compute fitness, in parallel:
		#ifdef _OPENMP
			#pragma omp parallel for num_threads(MAX_THREADS)
		#endif
		for (int i = int(pe); i < int(p); ++i) {
			if (i < int(pe) + nLocalSearch) {
				next.setFitness(i, refDecoder.decodeWithLocalSearch(next.population[i]));
			}
			else {
				next.setFitness(i, refDecoder.decode(next.population[i]));
			}
			
		}



		// Now we must sort 'current' by fitness, since things might have changed:
		next.sortFitness();
	}

	
};




