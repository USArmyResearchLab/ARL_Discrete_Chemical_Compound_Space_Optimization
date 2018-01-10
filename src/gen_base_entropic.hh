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
//! \file binary_entropic.hh \brief implementation of binary entropy driven optimization.

#ifndef _GEN_BASE_ENTROPIC_HH_
#define _GEN_BASE_ENTROPIC_HH_

#include <typedefs.hh>
#include <iostream>
#include <optimizeabstract.h>
#include <cmath>
#include <BCR_CPP_LA/linear_algebra.h>
#include <entropic_aux.hh>

using namespace std;

//! Meta-Optimize by generating maximally distant starting configurations.
template <class C>
class gen_base_entropic:
      public optimize_abstract
{

protected:
   C opt_object;
   refvector<long> bases;
   bool reorder;

public:
   //! Number of runs.
   ulong nruns;
   //! Maximum number of steps.
   long max_steps;

   gen_base_entropic(const C& a, const refvector<long>& b, bool _reorder=false) :
      opt_object(a),bases(b),nruns(2), reorder(_reorder) {};

   //! Meta-Optimize by generating maximally distant starting configurations.
   ulong optimize(ulong N) const
   {
      try {
         long i;
         ulong number=N;
         long config=opt_object.lib_object_r.visited_r.contains(N);
         refvector<double> lambda(opt_object.lib_object_r.get_number_of_constraints());

         opt_object.set_id(id_r+"::opt_object");

         // Actual optimization routine.
         ulong conf1;
         conf1=number;

         for(ulong runs=0;
               runs<nruns && max_steps>opt_object.lib_object_r.visited_r.size();
               runs++)
         {
            cout << id_r << " Starting run " << runs << " with " << conf1 << endl;
            conf1=opt_object.optimize(conf1);
            config=opt_object.lib_object_r.visited_r.contains(conf1);
            cout << id_r << ":The optimized value in run " << runs << " is: " << opt_object.lib_object_r.value_r[config].property;
            cout << " Penalty: ";
            opt_object.lib_object_r.value_r[config].penalty.display();
            cout << id_r << " for compound #" << conf1 << endl;
            cout.flush();

            refvector<ulong> library(opt_object.lib_object_r.visited_r.size());
            if (reorder)
               for (long i=0; i<library.size(); i++)
                  library[i]=opt_object.lib_object_r.reprune(opt_object.lib_object_r.visited_r[i]);
            else
               library=opt_object.lib_object_r.visited_r;
            // Single optimization run done.
            conf1=maximize_entropic_distance(library,bases);
            if (reorder)
               conf1=opt_object.lib_object_r.deprune(conf1);
         }
         cout << id_r << " Total Visited configurations ";
         opt_object.lib_object_r.visited_r.display();
         cout << endl;
         cout << id_r << "Number of configurations: " << opt_object.lib_object_r.visited_r.dim() << endl;


         conf1=0;
         for(i=0;i<opt_object.lib_object_r.visited_r.size();i++)
         {
            if(opt_object.lib_object_r.value_r[i].penalty*opt_object.lib_object_r.value_r[i].penalty<=
                  opt_object.lib_object_r.value_r[conf1].penalty*opt_object.lib_object_r.value_r[conf1].penalty)
               if(!(opt_object.lib_object_r.value_r[i].penalty==opt_object.lib_object_r.value_r[conf1].penalty) ||
                     opt_object.lib_object_r.value_r[i].property > opt_object.lib_object_r.value_r[conf1].property)
                  conf1=i;
         }
         conf1=opt_object.lib_object_r.visited_r[conf1];
         return conf1;
      } catch(exception& e) {
         cerr << e.what() << endl;
         throw domain_error("gen_base_entropic::optimize(ulong N) const");
      }
   }

   valerg get_value(const ulong i) const
   {
      return opt_object.lib_object_r.get_value(i);
   }

   void set_compute_property_flag(bool B) const {opt_object.set_compute_property_flag(B);}
      protected:

      private:

};

#endif // _GEN_BASE_ENTROPIC_HH_
