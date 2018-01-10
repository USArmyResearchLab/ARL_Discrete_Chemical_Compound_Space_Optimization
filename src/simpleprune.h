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
//! \file simpleprune.h \brief simple pruning class.

#ifndef _SIMPLEPRUNE_H_
#define _SIMPLEPRUNE_H_

#include <typedefs.hh>
#include <BCR_CPP_LA/refcount.decl>
#include <prunerabstract.h>

using namespace linear_algebra;

//! Class to prune a Library.
template <class X>
class simple_prune: public pruner_abstract<X>
{
private:

protected:

   using pruner_abstract<X>::value;
   using pruner_abstract<X>::visited;
   using pruner_abstract<X>::pruned_visited;
   using X::bits_computed;
   using X::space_size_computed;

public:

   using pruner_abstract<X>::value_r;
   using pruner_abstract<X>::visited_r;
   using pruner_abstract<X>::prune;

   simple_prune() : pruner_abstract<X>() {};
   simple_prune(const simple_prune<X>& a) :
      pruner_abstract<X>(a)
      {};
   //! Prune the Library and adjust lambda. Pruned entries are in pruned_visited.
   ulong prune(refvector<double> &lambda, ulong &conf1, ulong &conf2, long &config, const refvector<ulong>& a) const;
   //! Revert an index from the pruned Library to the original Library
   ulong deprune(ulong N) const;
   //! Convert an "absolute" index to a relative index.
   ulong reprune(ulong N) const;
protected:

private:
};

#endif // _SIMPLEPRUNE_H_
