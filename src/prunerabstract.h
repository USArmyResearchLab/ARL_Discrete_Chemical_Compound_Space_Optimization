/* -*- Mode: C++; indent-tabs-mode: t; c-basic-offset: 3; tab-width: 3 -*- */
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
//! \file prunerabstract.h \brief pruner_abstract implementation

#ifndef _PRUNERABSTRACT_H_
#define _PRUNERABSTRACT_H_

#include <BCR_CPP_LA/refcount.h>

//! Abstract class implementing the pruner concept.
template <class X>
class pruner_abstract : public X
{
private:
public:

   //! Default constructor
   pruner_abstract():
      X(),
      pruned_visited(),
      visited(X::visited),
      value(X::value),
      bits_computed(false),
      space_size_computed(false)
{ };

   //! Copy constructor
   pruner_abstract(const pruner_abstract<X>& a) :
      X(a),
      pruned_visited(a.pruned_visited),
      visited(X::visited),
      value(X::value),
      bits_computed(false),
      space_size_computed(false)
   {};
   pruner_abstract<X>& operator=(const pruner_abstract<X>& A)
   {
      (X&) *this = (X&) A;
      pruned_visited=A.pruned_visited;
      return *this;
   }

   //! Wrap the Library::compute_property to exclude pruned access.
   valerg compute_property(ulong N) const {
      return X::compute_property(deprune(N));
   }

   //! Adjust the space_size to reflect the absence of pruned values
   ulong get_space_size() const {
      if(!space_size_computed) {
         space_size=X::get_space_size();
         space_size-=pruned_visited.size();
         space_size_computed=true;
      }
      return space_size; }

   //! Compute the bits of the adjusted space.
   ulong get_bits() const {
      if(!bits_computed) {
         get_space_size();
         bits_computed=true;
         bits=(ulong) ceil(log((double) space_size)/log(2.0));
      }
      return bits;
   }

   //! Same as compute_property(ulong N), but absolute numbering
   valerg get_value(ulong N) const { return X::get_value(N); }

   virtual ulong prune(refvector<double>& lambda, ulong &conf1, ulong &conf2, long &config, const refvector<ulong>& visit) const=0;
   ulong prune(refvector<double>& lambda, ulong &conf1, ulong &conf2, long &config) const
   {
      return prune(lambda,conf1,conf2,config,X::visited_r);
   }
   //! Given an "absolute" reference number return the current reference number.
   virtual ulong reprune(ulong N) const=0;
   //! self-explanatory
   virtual ulong deprune(ulong N) const=0;

protected:
   //! List of pruned indices in Library.
   mutable	refvector<ulong> pruned_visited;

   //! Reference to Library::visited.
   refvector<ulong>& visited;

   //! Reference to Library::value
   refvector<valerg>& value;
   mutable bool bits_computed;
   mutable bool space_size_computed;
   mutable ulong space_size;
   mutable ulong bits;
private:
};

#endif // _PRUNERABSTRACT_H_
