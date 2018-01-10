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
//! \file binary_entropic.hh Definition and Implementation of the binary_entropic class.

#ifndef _BINARYENTROPIC_HH_
#define _BINARYENTROPIC_HH_
#include <typedefs.hh>
#include <iostream>
#include <optimizeabstract.h>
#include <cmath>
#include <BCR_CPP_LA/linear_algebra.h>
#include <entropic_aux.hh>

using namespace std;
using namespace linear_algebra;

//! Class for enhanced sampling using an entropic measure of coverage.
template<class C>
class binary_entropic: public optimize_abstract
/*!
   This class provides a meta-optimization mechanism.
   The secondary optimizer C has to comply with the has_optimize and
   pruner concepts in order for the optimization to work.
   After each successful optimization the visited configurations are analysed
   for coverage of the library. maximize_entropic_distance() is then used to determine
   the next best starting point.
 */
{

   void print_finished_optimization(long config, ulong conf1) const
   {
      cout << "The optimized value is: " << opt_object.value_r[config].property << endl;
      cout << " Penalty: ";
      opt_object.value_r[config].penalty.display(cout);
      cout << " for compound #" << conf1 << endl;
      cout.flush();
   }

protected:
   //! Optimization object of type C. This is the underlying optimization method.
   C opt_object;
   //! This vector defines the substitution bases used in the optimization and entropic maximization ( maximize_entropic_distance())
   refvector<long> bases;

public:
   //! Number of runs.
   ulong nruns;
   //! Maximum number of steps.
   long max_steps;

   //! Copy constructor
   binary_entropic(const C& a) : opt_object(a), nruns(1), max_steps(0)
   {
      refvector<long> b(opt_object.get_bits());
      for(long i=0;i<b.dim();i++)
         b[i]=2;
      bases=b;
   };

   //! Meta-Optimize by generating maximally distant starting configurations.
   ulong optimize(ulong N) const
   /*!
      The metric used to determine the distance is implemented and explained
      in maximize_entropic_distance().
    */
   {
      try {
         long i;
         ulong number=N;
         long config=opt_object.visited_r.contains(N);

         // Actual optimization routine.
         ulong conf1;
         conf1=number;

         for(ulong runs=0;
               runs<nruns && max_steps>opt_object.visited_r.size();
               runs++)
         {
            conf1=opt_object.optimize(conf1);
            config=opt_object.visited_r.contains(opt_object.deprune(conf1));
            print_finished_optimization(config, conf1);
            // Single optimization run done.
            conf1=maximize_entropic_distance(opt_object.visited_r,bases);
         }

         // minimize penalty (assumed to be non-negative), then maximize property
         conf1=0;
         for(i=0;i<opt_object.visited_r.size();i++)
         {
            if(opt_object.value_r[i].penalty*opt_object.value_r[i].penalty<=
                  opt_object.value_r[conf1].penalty*opt_object.value_r[conf1].penalty)
               if(!(opt_object.value_r[i].penalty==opt_object.value_r[conf1].penalty) ||
                     opt_object.value_r[i].property > opt_object.value_r[conf1].property)
                  conf1=i;
         }
         conf1=opt_object.visited_r[conf1];
         return conf1;
      } catch(exception& e) {
         cerr << e.what() << endl;
         throw domain_error("binary_entropic::optimize(ulong N) const");
      }
   }

   valerg get_value(const ulong i) const
   {
      return opt_object.get_value(i);
   }

   void set_compute_property_flag(bool B) const {opt_object.set_compute_property_flag(B);}
protected:

private:

};

#endif // _BINARYENTROPIC_HH_
