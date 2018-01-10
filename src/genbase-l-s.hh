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
//! \file genbase-l-s.hh Line search optimization with general bases.

#ifndef _GENBASE_L_S_HH_
#define _GENBASE_L_S_HH_

#include <typedefs.hh>
#include <iostream>
#include <optimizeabstract.h>
#include <BCR_CPP_LA/refcount.h>

using namespace std;
using namespace linear_algebra;

//! Line search optimization class using general bases.
template <class C, class B>
class gen_base_LS: public optimize_abstract
{

   using optimize_abstract::optimize;

   const C lib_object;
   mutable B bases;

   void select_current_best(const valerg& interim, const refvector<double>& lambda, valerg& current_best_val, ulong& conf1, ulong nm, long& config) const
   {
      if (interim.property_computed && interim.property - lambda * interim.penalty >= current_best_val.property - lambda * current_best_val.penalty) {
         cout << " > ";
         cout << lib_object.deprune(conf1) << "(" << lib_object.visited_r.contains(lib_object.deprune(conf1)) << ") = " << current_best_val.property << " and penalty: ";
         current_best_val.penalty.display();

         current_best_val = interim;
         conf1 = nm;
         config = lib_object.visited_r.contains(lib_object.deprune(conf1));
      } else {
         config = lib_object.visited_r.contains(lib_object.deprune(conf1));
         cout << " < ";
         cout << lib_object.deprune(conf1) << "(" << config << ") = " << current_best_val.property << " and penalty: ";
         current_best_val.penalty.display();
      }
   }

   void update_visited_run(refvector<ulong>& visited_run, ulong nm, const valerg& interim) const
   {
      if (interim.property_computed && visited_run.contains(lib_object.deprune(nm)) < 0)
         visited_run.push_back(lib_object.deprune(nm));

      cout << id_r << "::Config: " << lib_object.deprune(nm) << "(" << lib_object.visited_r.contains(lib_object.deprune(nm)) << ")  finished with property: " << interim.property << " and penalty: ";
      interim.penalty.display();
   }

public:

   const B& bases_r;
   const C& lib_object_r;

   gen_base_LS(const C& Library):
      lib_object(Library),
      bases(lib_object), bases_r(bases),
      lib_object_r(lib_object)
   {};

   //! Optimize starting from conformation i.
   /*!
   A Lagrange multiplier method is employed to enforce
   boundary constraints. The multiplier increases as the
   optimizatoin proceeds, enforcing the constraint ever
   more rigorously.
    */
   ulong optimize(ulong N) const
   {
      try
      {
         ulong number=lib_object_r.reprune(N);
         valerg current_best_val;
         refvector<ulong> visited_run;
         long config=lib_object.visited_r.contains(N);

         current_best_val=lib_object.compute_property(number);
         while(current_best_val.energy==INFINITY && number<lib_object_r.get_space_size()-1)
         {
            number++;
            current_best_val=lib_object.compute_property(number);
         }
         if(current_best_val.property_computed)
            visited_run.push_back(number);
         refvector<double> lambda(lib_object_r.get_number_of_constraints());

         // Actual optimization routine.
         ulong conf1;
         ulong conf2;
         ulong nm;
         valerg interim;
         conf1=number;
         conf2=conf1-1;
         while (conf1!=conf2)
         {
            conf2=conf1;
            long j;
            ulong k=1;
            bases.set_refstate(lib_object.deprune(conf1));

            for(bases=0;!bases.done();bases++)
            {
               cout << "In "<< id_r << "::optimize(): " << bases.get_state() << endl;

               ulong conf3=conf1-(((conf1-conf1%bases())/bases()) % bases.modulus()) * bases();

               for(j=0;j<bases.modulus();j++)
               {
                  nm=conf3+j*bases();

                  interim=lib_object.compute_property(nm);
                  update_visited_run(visited_run, nm, interim);
                  select_current_best(interim, lambda, current_best_val, conf1, nm, config);
               }
               bases.set_refstate(conf1);
            }

            config=lib_object.visited_r.contains(lib_object.deprune(conf1));
            cout << id_r+"::optimized value is: " << lib_object.value_r[config].property;
            cout << " Penalty: ";
            lib_object.value_r[config].penalty.display();
            cout << " lambda: ";
            lambda.display();
            cout << " Result: " << (lib_object.value_r[config].property-
                  lib_object.value_r[config].penalty*lambda)
                       << " for compound #" << lib_object.deprune(conf1) << endl;
            cout.flush();

            // Single optimization run done.

            lib_object.prune(lambda,
                  conf1,
                  conf2,
                  config, visited_run);

            current_best_val=lib_object.value_r[config];
            lambda*=1.1;
            cout << id_r << "::New lambda = ";
            lambda.display();
            lib_object.get_space_size();
         }
         return lib_object.deprune(conf1);
      }
      catch (domain_error e)
      {
         cerr << e.what() << endl;
         throw domain_error("called by gen_base_LS::optimize(ulong N) const");
      }
   }

   valerg get_value(const ulong i) const
   {
      return lib_object.get_value(i);
   }
   void set_compute_property_flag(bool Bo) const {lib_object.set_compute_property_flag(Bo);}
protected:

private:

};

#endif // _GENBASE_L_S_HH_
