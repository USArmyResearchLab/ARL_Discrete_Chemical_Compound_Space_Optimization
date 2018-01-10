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
/*! \file has_gradients_hessians_data.hh */
/*! \brief Implementations of the has_gradients and has_hessians concepts
 */
/*!@{*/
#ifndef _HAS_GRADIENTS_HESSIANS_DATA
#define _HAS_GRADIENTS_HESSIANS_DATA

#include <has_gradients_hessian_data.decl>
#include <BCR_CPP_LA/refcount.h>
#include <BCR_CPP_LA/mat_sym_full.h>
#include <typedefs.hh>
#include <stdexcept>

using namespace linear_algebra;

//! Gradient computation
template<class X>
refvector<valerg> has_gradients_data<X>::gradient(ulong i) const
{
   refvector<valerg> v(get_bits());
   return gradient(i,v);
}
//! Gradient computation
template<class X>
refvector<valerg>& has_gradients_data<X>::gradient(ulong i, refvector<valerg>& v) const
{
   try {
      v.resize(get_bits());
      long maxN=get_space_size();
      long n1, n2, m, k;
      valerg a,b;
      long lbits=0;
      for(k=1;k<maxN;k*=2,lbits++)
      {
         m=((i-i%k)/k)%2;
         n1=i-m*k;
         n2=n1+k;
         a=X::compute_property(n1);
         b=X::compute_property(n2);
         v[lbits]=b-a;
      }
   }
   catch (exception& e) {
      cerr << e.what() << endl;
      throw domain_error("called by has_gradients_data<X>::gradient(ulong i, refvector<valerg>& v) const");
   }
   return v;
}
//! Hessian computation
template<class X>
mat_sym_full<valerg> has_hessians_data<X>::hessian(ulong i) const
{
   mat_sym_full<double> H(X::get_bits());
   return hessian(i,H);
}
//! Hessian computation
template<class X>
mat_sym_full<valerg>& has_hessians_data<X>::hessian(ulong i, mat_sym_full<valerg>& H) const
{
   try {
      H.resize(X::get_bits());
      long maxN=X::get_space_size();
      long n1, n2, m;
      long l, n11, n12, n21, n22;
      double p11,p12,p21,p22;
      long bitsk,bitsl,k;
      for(bitsk=0,k=1;k<maxN;k*=2,bitsk++)
      {
         m=((i-i%k)/k)%2;
         n1=i-m*k;
         n2=n1+k;
         for(bitsl=0,l=1;l<maxN;l*=2,bitsl++)
         {
            m=((n1-n1%k)/k)%2;
            n11=n1-m*k;
            n12=n1+k;
            m=((n2-n2%k)/k)%2;
            n21=n2-m*k;
            n22=n2+k;
            p11=X::compute_property(n11).property;
            p12=X::compute_property(n12).property;
            p21=X::compute_property(n21).property;
            p22=X::compute_property(n22).property;
            H[bitsk*(bitsk+1)/2+bitsl]=p11+p22-p12-p21;
         }
      }
   } catch (exception& e) {
      cerr << e.what() << endl;
      throw domain_error("called by has_hessians_data<X>::hessian(ulong i, mat_sym_full<valerg>& H) const");
   }
   return H;
}

#endif
/*! @}*/
