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

//! \file generalbaseiterator.hh

#ifndef _GENERALBASEITERATOR_HH_
#define _GENERALBASEITERATOR_HH_

#include <BCR_CPP_LA/refcount.h>
#include <typedefs.hh>
#include <chemident.hh>

//! Iterator over the potential bases
template<class X>
class general_base_iterator
{
public:
   //! Constructor. Mandatory library object.
   general_base_iterator(const X& b);
   //! Copy Constructor.
   general_base_iterator(const general_base_iterator<X>& b);
   //! Advance the iterator to the next number.
   general_base_iterator<X>& operator++(int i);
   //! Set the reference state of the library to determine the relevance of subs.
   ulong set_refstate(ulong newref);
   //! Retrieve the current refstate.
   ulong get_refstate() const;
   //! Set the current state.
   long operator=(long i);
   //! Retrieve the current state.
   long get_state() const;
   //! Retrieve the current base size.
   long operator()() const;
   //! Modulus
   long modulus() const;
   //! End of iterator?
   bool done() const
   {
      if(state<occupation.size())
         return false;
      else
         return true;
   }
   //!
   long non_empty_size() const;
protected:

private:

   long number_of_bases;
   //! Reference state
   ulong refstate;
   //! Current base
   long state;
   //! Library object which is enumerated with the bases to be iterated.
   const X& lib_object;
   //! Offsets of the bases.
   mutable refvector<long> bases;
   //! Moduli of the bases.
   mutable refvector<long> moduli;
   //! Compute offsets of groups
   void compute_bases();
   //! Compute the base for N.
   void compute_bases(long N, const refvector<long>& base_sizes);
   //! occupation of groups by current refstate.
   refvector<long> occupation;
   //! Set the occupation numbers w.r.t. number.
   void occupy(ulong number);
};

#endif // _GENERALBASEITERATOR_HH_
