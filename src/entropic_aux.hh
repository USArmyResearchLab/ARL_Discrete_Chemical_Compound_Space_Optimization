/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 3; tab-width: 3 -*- */
/*
 * This file is part of ARL Discrete Chemical Compound Space Optimization (ARL DCCSO) project.
 *
 * ARL DCSSO constitutes a work of the United States Government and is not
 * subject to domestic copyright protection under 17 USC Sec. 105.
 * Release authorized by the US Army Research Laboratory
 *
 * To the extent possible under law, the author(s) have dedicated all copyright
 * and related and neighboring rights to this software to the public domain
 * worldwide. This software is distributed without any warranty.
 *
 * You should have received a copy of the CC0 Public Domain Dedication along
 * with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
 *
 */

/*! @addtogroup DOChemS Discrete Optimization of Chemical Space */
//! \file entropic_aux.hh Auxiliary functions for *_entropic classes.

#ifndef _ENTROPIC_AUX_HH_
#define _ENTROPIC_AUX_HH_
#include <iostream>
#include <cmath>
#include <BCR_CPP_LA/linear_algebra.h>
#include <typedefs.hh>

using namespace std;
using namespace linear_algebra;

void linsolve_cg(const mat_sym_full<double>& J,
      const refvector<double>& G,
      refvector<double>& X);

refvector<double>& set_gradient(const refvector<long>& b,
      const mat_full<double>& H,
      const refvector<double>& X,
      refvector<double>& G);

mat_sym_full<double>& set_hessian(const refvector<long>& b,
      const mat_full<double>& H,
      const refvector<double>& X,
      mat_sym_full<double>& J);

ulong maximize_entropic_distance(const refvector<ulong>& A, const refvector<long>& b);

#endif // _ENTROPIC_AUX_HH_
