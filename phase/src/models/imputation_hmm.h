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

#ifndef _HAPLOTYPE_HMM_H
#define _HAPLOTYPE_HMM_H

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


class imputation_hmm {
private:
	conditioning_set * C;
	unsigned int modK;

	//DYNAMIC ARRAYS
	avx_aligned_vector < float > Emissions;
	avx_aligned_vector < float > Alpha;
	avx_aligned_vector < float > AlphaSum;
	avx_aligned_vector < float > Beta;

public:
	//CONSTRUCTOR/DESTRUCTOR
	imputation_hmm(conditioning_set *);
	~imputation_hmm();

	void resize();
	void init(const std::vector < float > &);
	void forward(std::vector < bool > &);
	void backward(const std::vector < float > &, std::vector < bool > &, std::vector < float > &);
	void computePosteriors(const std::vector < float > &, std::vector < bool > &, std::vector < float > &);

};

#endif
