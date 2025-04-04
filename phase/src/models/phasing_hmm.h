/*******************************************************************************
 * Copyright (C) 2022-2023 Simone Rubinacci
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * MIT Licence
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef _DIPLOTYPE_HMM_H
#define _DIPLOTYPE_HMM_H

#include <utils/otools.h>
#include <containers/conditioning_set.h>
#include <immintrin.h>
#include <boost/align/aligned_allocator.hpp>
#include <utils/ug.h>

template <typename T>
#if UG_AVX512
using avx_aligned_vector = std::vector<T, boost::alignment::aligned_allocator < T, 64 > >;
#else
using avx_aligned_vector = std::vector<T, boost::alignment::aligned_allocator < T, 32 > >;
#endif

#define HAP_NUMBER 8

#define VAR_PEAK_HET 0
#define VAR_PEAK_HOM -1
#define VAR_FLAT_HET -2

#define ALLELE(hap , pos) ((hap) & (1<<(pos)))

inline
float horizontal_add (const __m256& a) //@simo: implemented horizontal add with instructions in the 4B minimum
{
    __m128 vlow = _mm256_castps256_ps128(a);
    __m128 vhigh = _mm256_extractf128_ps(a, 1); // high 128
   vlow = _mm_add_ps(vlow, vhigh);     // add the low 128
   __m128 shuf = _mm_movehdup_ps(vlow);        // broadcast elements 3,1 to 2,0
   __m128 sums = _mm_add_ps(vlow, shuf);
   shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
   sums = _mm_add_ss(sums, shuf);
   return _mm_cvtss_f32(sums);
}

class phasing_hmm {
private:

	conditioning_set * C;

	//EXTERNAL DATA
	std::vector < char > VAR_TYP;
	std::vector < bool > VAR_ALT;
	std::vector < int > VAR_ABS;
	std::vector < int > VAR_REL;

	//COORDINATES & CONSTANTS
	unsigned int n_segs;
	unsigned int n_miss;

	//SEGMENTATION
	std::vector < int > segments;

	//CURSORS
	int curr_idx_locus;
	int curr_abs_locus;
	int curr_rel_locus;
	int curr_segment_index;
	int curr_segment_locus;
	int curr_missing_locus;

	//DYNAMIC ARRAYS
	float probSumT;
	avx_aligned_vector < float > prob;
	avx_aligned_vector < float > probSumK;
	avx_aligned_vector < float > probSumH;

	avx_aligned_vector < float > phasingProb;
	avx_aligned_vector < float > phasingProbSum;
	std::vector < float > phasingProbSumSum;

	avx_aligned_vector < float > imputeProb;
	avx_aligned_vector < float > imputeProbSum;
	avx_aligned_vector < float > imputeProbSumSum;
	avx_aligned_vector < float > imputeProbOf1s;
	std::vector < int > dip_sampled;

	//STATIC ARRAYS
	std::vector < float > DProbs;
	std::vector<avx_aligned_vector<float> > EMIT0;
	std::vector<avx_aligned_vector<float> > EMIT1;
	avx_aligned_vector < float > HProbs;
	float sumHProbs, sumDProbs;
	float nt, yt;

	//INLINED ROUTINES
	void INIT_PEAK_HET(int);
	void INIT_PEAK_HOM(bool);
	void INIT_FLAT_HET();

	void RUN_PEAK_HET(int);
	void RUN_PEAK_HOM(bool);
	void RUN_FLAT_HET();

	void COLLAPSE_PEAK_HET(int);
	void COLLAPSE_PEAK_HOM(bool);
	void COLLAPSE_FLAT_HET();

    void SUMK();
    bool TRANS_HAP();
    bool SAMPLE_DIP();
    void IMPUTE_FLAT_HET();

public:
	//CONSTRUCTOR/DESTRUCTOR
	phasing_hmm(conditioning_set * C);
	~phasing_hmm();

	void reallocate(const std::vector < bool > &, const std::vector < bool > &, std::vector < bool > &);
	void forward();
	void backward();
	void rephaseHaplotypes(std::vector < bool > &, std::vector < bool > &, std::vector < bool > &);
};

inline
void phasing_hmm::INIT_PEAK_HET(int curr_het)
{
	const std::array <__m256, 2 > emits = {_mm256_load_ps(&EMIT0[curr_het][0]),_mm256_load_ps(&EMIT1[curr_het][0])};
    __m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
	{
		const bool ah = C->Hvar.get(curr_rel_locus, k);
		_sum = _mm256_add_ps(_sum, emits[ah]);
		_mm256_store_ps(&prob[i], emits[ah]);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
}

inline
void phasing_hmm::INIT_PEAK_HOM(bool ag)
{
	const std::array <__m256, 2 > emits = {_mm256_set1_ps(1.0f),_mm256_set1_ps(C->ed_phs/C->ee_phs)};
    __m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
	{
		const bool ag_ah = C->Hvar.get(curr_rel_locus, k)!=ag;
		_sum = _mm256_add_ps(_sum, emits[ag_ah]);
		_mm256_store_ps(&prob[i], emits[ag_ah]);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
}

inline
void phasing_hmm::INIT_FLAT_HET() {
	fill(prob.begin(), prob.end(), 1.0f/(HAP_NUMBER * C->n_states));
	fill(probSumH.begin(), probSumH.end(), 1.0f/HAP_NUMBER);
    probSumT = 1.0f;
}

inline
void phasing_hmm::RUN_PEAK_HET(int curr_het)
{
	const __m256 _tFreq = _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
	const __m256 _nt = _mm256_set1_ps(nt / probSumT);
	const std::array <__m256, 2 > emits = {_mm256_load_ps(&EMIT0[curr_het][0]),_mm256_load_ps(&EMIT1[curr_het][0])};
    __m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
	{
		const bool ah = C->Hvar.get(curr_rel_locus, k);
		const __m256 _prob_curr = _mm256_mul_ps(_mm256_fmadd_ps(_mm256_load_ps(&prob[i]), _nt, _tFreq), emits[ah]);
		_sum = _mm256_add_ps(_sum, _prob_curr);
		_mm256_store_ps(&prob[i], _prob_curr);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
}

inline
void phasing_hmm::RUN_PEAK_HOM(bool ag)
{
#if UG_AVX512
	assert(UG_BITMATRIX);
	const __m256 _tFreq256= _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
	__m512 _tFreq = _mm512_set1_ps(0.0f);					// broadcast 0,1,2,3,4,5,6,7 -> 8,9,10,11,12,13,14,15
	_tFreq = _mm512_insertf32x8(_tFreq, _tFreq256, 0);
	_tFreq = _mm512_insertf32x8(_tFreq, _tFreq256, 1);
	const __m512 _nt = _mm512_set1_ps(nt / probSumT);

    const __m256 _mism256 = _mm256_set1_ps(C->ed_phs/C->ee_phs);
    const __m512 _mism = _mm512_set1_ps(C->ed_phs/C->ee_phs);
	__m512 _mismX1 = _mm512_set1_ps(1.0f);
	_mismX1 = _mm512_insertf32x8(_mismX1, _mism256, 0);
	__m512 _mism1X = _mm512_set1_ps(1.0f);
	_mism1X = _mm512_insertf32x8(_mism1X, _mism256, 1);

    __m512 _sum = _mm512_set1_ps(0.0f);

	unsigned char* ptr = C->Hvar.getBytePtr(curr_rel_locus, 0); // cache row pointer
	unsigned int byteAcc; // 8 bit byte accumulator
	int size = C->n_states;
	int size2 = size & ~1; // size up to last even number

	// 2-state (16 floats) loop
	int k, i;
	for(k = 0, i = 0 ; k != size2 ; k += 2, i += (2 * HAP_NUMBER))
	{
		// on mod8==0, load accumulator
		if ( !(k & 7) ) {
			byteAcc = *ptr++;
		}

		// extract 2 values value from byte accumulator
		const bool ah0 = byteAcc & 0x80;
		const bool ah1 = byteAcc & 0x40;
		byteAcc <<= 2;

		const __m512 _prob_prev = _mm512_load_ps(&prob[i]);
		__m512 _prob_curr = _mm512_fmadd_ps(_prob_prev, _nt, _tFreq);
		if ( ah0 == ah1 ) {
			if ( ag != ah0 )
				_prob_curr = _mm512_mul_ps(_prob_curr, _mism);
		} else if ( ag != ah0 ) {
			_prob_curr = _mm512_mul_ps(_prob_curr, _mismX1);
		} else {
			_prob_curr = _mm512_mul_ps(_prob_curr, _mism1X);
		}

		_sum = _mm512_add_ps(_sum, _prob_curr);
		_mm512_store_ps(&prob[i], _prob_curr);
	}

	// fold sum back into 256 bits
	__m256 _sum256 = _mm512_castps512_ps256(_sum);
	__m256 _sum_high = _mm512_extractf32x8_ps(_sum, 1);
	_sum256 = _mm256_add_ps(_sum256, _sum_high);

	// reminder loop (should be exactly 1 iteration or none)
	if ( k != size ) {
		assert((k+1) == size);
		const __m256 _tFreq = _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
		const __m256 _nt = _mm256_set1_ps(nt / probSumT);
    	const __m256 _mism = _mm256_set1_ps(C->ed_phs/C->ee_phs);

		// extract value from byte accumulator
		if ( !(k & 7) ) {
			byteAcc = *ptr++;
		}
		const bool ah = byteAcc & 0x80;

		const __m256 _prob_prev = _mm256_load_ps(&prob[i]);
		__m256 _prob_curr = _mm256_fmadd_ps(_prob_prev, _nt, _tFreq);
		if (ag!=ah) {
			_prob_curr = _mm256_mul_ps(_prob_curr, _mism);
		}
		_sum256 = _mm256_add_ps(_sum256, _prob_curr);
		_mm256_store_ps(&prob[i], _prob_curr);
	}

	_mm256_store_ps(&probSumH[0], _sum256);
	probSumT = horizontal_add(_sum256);

#else
	const __m256 _tFreq = _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
	const __m256 _nt = _mm256_set1_ps(nt / probSumT);
    const __m256 _mism = _mm256_set1_ps(C->ed_phs/C->ee_phs);
    __m256 _sum = _mm256_set1_ps(0.0f);
#if UG_BITMATRIX
	unsigned char* ptr = C->Hvar.getBytePtr(curr_rel_locus, 0); // cache row pointer
	unsigned int byteAcc; // 8 bit byte accumulator
#endif
	int size = C->n_states;
	for(int k = 0, i = 0 ; k != size ; ++k, i += HAP_NUMBER)
	{
#if UG_BITMATRIX
		// on mod8==0, load accumulator
		if ( !(k & 7) ) {
			byteAcc = *ptr++;
		}

		// extract value from byte accumulator
		const bool ah = byteAcc & 0x80;
		byteAcc <<= 1;
#else
		const bool ah = C->Hvar.get(curr_rel_locus, k);
#endif

		const __m256 _prob_prev = _mm256_load_ps(&prob[i]);
		__m256 _prob_curr = _mm256_fmadd_ps(_prob_prev, _nt, _tFreq);
		if (ag!=ah) {
			_prob_curr = _mm256_mul_ps(_prob_curr, _mism);
		}
		_sum = _mm256_add_ps(_sum, _prob_curr);
		_mm256_store_ps(&prob[i], _prob_curr);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
#endif
}

inline
void phasing_hmm::RUN_FLAT_HET()
{
	const __m256 _tFreq = _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
	const __m256 _nt = _mm256_set1_ps(nt / probSumT);
    __m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
	{
		const __m256 _prob_curr = _mm256_fmadd_ps(_mm256_load_ps(&prob[i]), _nt, _tFreq);
		_sum = _mm256_add_ps(_sum, _prob_curr);
		_mm256_store_ps(&prob[i], _prob_curr);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
}

inline
void phasing_hmm::COLLAPSE_PEAK_HET(int curr_het)
{
	const __m256 _tFreq = _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
	const __m256 _nt = _mm256_set1_ps(nt / probSumT);
	const std::array <__m256, 2 > emits = {_mm256_load_ps(&EMIT0[curr_het][0]),_mm256_load_ps(&EMIT1[curr_het][0])};
    __m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
	{
		const bool ah = C->Hvar.get(curr_rel_locus, k);
		const __m256 _prob_curr = _mm256_mul_ps(_mm256_fmadd_ps(_mm256_set1_ps(probSumK[k]), _nt, _tFreq), emits[ah]);
		_sum = _mm256_add_ps(_sum, _prob_curr);
		_mm256_store_ps(&prob[i], _prob_curr);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
}

inline
void phasing_hmm::COLLAPSE_PEAK_HOM(bool ag)
{
	const __m256 _tFreq = _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
   	const __m256 _nt = _mm256_set1_ps(nt / probSumT);
    __m256 _sum = _mm256_set1_ps(0.0f);
    const __m256 _mism = _mm256_set1_ps(C->ed_phs/C->ee_phs);

   	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
   	{
   		const bool ah = C->Hvar.get(curr_rel_locus, k);
   		__m256 _prob_curr = _mm256_fmadd_ps(_mm256_set1_ps(probSumK[k]), _nt, _tFreq);
		if (ag!=ah) _prob_curr = _mm256_mul_ps(_prob_curr, _mism);
   		_sum = _mm256_add_ps(_sum, _prob_curr);
   		_mm256_store_ps(&prob[i], _prob_curr);
   	}
   	_mm256_store_ps(&probSumH[0], _sum);
   	probSumT = horizontal_add(_sum);
}

inline
void phasing_hmm::COLLAPSE_FLAT_HET ()
{
	const __m256 _tFreq = _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
	const __m256 _nt = _mm256_set1_ps(nt / probSumT);
    __m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
	{
   		__m256 _prob_curr = _mm256_fmadd_ps(_mm256_set1_ps(probSumK[k]), _nt, _tFreq);
		_sum = _mm256_add_ps(_sum, _prob_curr);
		_mm256_store_ps(&prob[i], _prob_curr);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
}


inline
void phasing_hmm::SUMK() {
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
		probSumK[k] = horizontal_add(_mm256_load_ps(&prob[i]));
}

inline
bool phasing_hmm::TRANS_HAP()
{
	const int states_haps = C->n_states*HAP_NUMBER;
	sumHProbs = 0.0f;
	yt = C->getTransition(VAR_ABS[curr_idx_locus], VAR_ABS[curr_idx_locus+1]);
	nt = 1.0f - yt;
	const __m256 _fact2 = _mm256_set1_ps( nt / probSumT);
	int h1 = 0, j=0;
	for (; h1 < HAP_NUMBER ; h1++, j += HAP_NUMBER)
	{
		const __m256 _fact1 = _mm256_set1_ps((probSumH[h1]/probSumT) * yt / C->n_states);
		__m256 _sum = _mm256_set1_ps(0.0f);
		for(int k=0, i=0; k != C->n_states ; ++k, i += HAP_NUMBER)
		{
			const __m256 _prob0 = _mm256_fmadd_ps(_mm256_set1_ps(prob[i+h1]), _fact2, _fact1);
			_sum = _mm256_add_ps(_sum, _mm256_mul_ps(_prob0, _mm256_load_ps(&phasingProb[(curr_segment_index+1)*states_haps+i])));
		}
		_mm256_store_ps(&HProbs[j], _sum);
		sumHProbs += horizontal_add(_sum);
	}
	return (std::isnan(sumHProbs) || std::isinf(sumHProbs) || sumHProbs < std::numeric_limits<float>::min());
}

inline
bool phasing_hmm::SAMPLE_DIP() {
	sumDProbs = 0.0f;
	for (int d = 0 ; d < HAP_NUMBER ; d ++) {
		int prev_h0 = dip_sampled[curr_segment_index];
		int prev_h1 = HAP_NUMBER - dip_sampled[curr_segment_index] - 1;
		DProbs[d] = (HProbs[prev_h0 * HAP_NUMBER + d] / sumHProbs) * (HProbs[prev_h1 * HAP_NUMBER + (HAP_NUMBER - d - 1)] / sumHProbs);
		sumDProbs += DProbs[d];
	}
	if (std::isnan(sumDProbs) || std::isinf(sumDProbs) || sumDProbs < std::numeric_limits<float>::min()) return true;
	dip_sampled[curr_segment_index+1] = rng.sample(DProbs, sumDProbs);
	return false;
}

inline
void phasing_hmm::IMPUTE_FLAT_HET()
{
	const int states_haps = C->n_states*HAP_NUMBER;
	const __m256 _one = _mm256_set1_ps(1.0f);
	const __m256 _zero = _mm256_set1_ps(0.0f);
	__m256 _scaleR = _mm256_load_ps(&imputeProbSum[curr_missing_locus*HAP_NUMBER]);
	__m256 _scaleL = _mm256_load_ps(&probSumH[0]);
	_scaleR = _mm256_div_ps(_one, _scaleR);
	_scaleL = _mm256_div_ps(_one, _scaleL);
	std::array <__m256, 2 > sums = {_mm256_set1_ps(0.0f),_mm256_set1_ps(0.0f)};

	for(int k = 0, i = 0 ; k !=  C->n_states ; ++k, i += HAP_NUMBER) {
		const bool ah = C->Hvar.get(curr_rel_locus, k);
		const __m256 _p1 = _mm256_mul_ps(_mm256_load_ps(&imputeProb[curr_missing_locus*states_haps + i]), _scaleR);
		const __m256 _p2 = _mm256_mul_ps(_mm256_load_ps(&prob[i]), _scaleL);
		sums[ah] = _mm256_add_ps(sums[ah], _mm256_mul_ps(_p1,_p2));
	}
	const __m256 _norm_sum = _mm256_div_ps(sums[1], _mm256_add_ps(sums[0], sums[1]));
	const __m256 _clamp = _mm256_max_ps(_zero, _mm256_min_ps(_one, _norm_sum));
	_mm256_store_ps(&imputeProbOf1s[curr_missing_locus*HAP_NUMBER], _clamp);
}

#endif

