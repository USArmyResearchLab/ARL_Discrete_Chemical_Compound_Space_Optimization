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
//! \file genbasegdmc.hh \brief implementation of GDMC on a hypertorus

#ifndef _GENBASEGDMC_HH_
#define _GENBASEGDMC_HH_

#include <typedefs.hh>
#include <iostream>
#include <optimizeabstract.h>
#include <BCR_CPP_LA/refcount.h>

using namespace std;

template <class C>
class gen_base_gdmc: public optimize_abstract
/*!
   This class takes a Library object and performs an optimization.
   The Library is arranged in a hypertorus topology in which each molecule occupies
   a vertex on the hull of the torus. After a local optimum has been reached, a
   Monte-Carlo simulation is performed to determine a new starting point.
   Each time a new minimum is found the acceptance criteria become more strict.
 */
{

private:
   //! Optimization object
   C opt_object;
   //! Definition of the hypertorus topology.
   refvector<long> bases;

public:
   //! Temperature.
   double T;
   //! Number of tightening steps.
   ulong tight_steps;
   //! Maximum number of steps.
   ulong max_steps;

   gen_base_gdmc(const C& a, const refvector<long>& b) :
      opt_object(a),
      bases(b),
      T(0.0),
      tight_steps(1),
      max_steps(2)
   {};

   //! Gradient-directed Monte-Carlo optimization
   ulong optimize(ulong N) const
   {
      try {
         ulong i,j;
         long k;
         ulong number=N;
         long config=opt_object.lib_object_r.visited_r.contains(N);

         refvector<double> lambda(opt_object.lib_object_r.get_number_of_constraints());
         ulong steps=1;
         opt_object.set_id(id_r+"::opt_object");

         // Actual optimization routine.
         ulong conf1,conf2;
         conf1=number;

         refvector<valerg> valgrad(bases.size());

         refvector<double> grad(valgrad.dim());

         ulong conf3;
         // tight_steps counts how often lambda is ramped.
         while(tight_steps>=steps)
         {
            // To force lambda ramping the average number of computed values is limited.
            while(steps*max_steps/tight_steps>(ulong) opt_object.stacksize())
            {
               conf2=conf1;
               //This reinitializes the optimization intermediates.
               conf1=opt_object.optimize(conf1);
               cout << id_r+"::Gradient of " << conf1 << endl;
               cout.flush();
               valgrad=opt_object.gradient(opt_object.lib_object_r.reprune(conf1),valgrad);
               number=0;

               for(k=0,i=1;k<valgrad.size();k++)
               {

                  j=(conf1-conf1%i)/i % bases[i];

                  grad[k]=valgrad[k].property-lambda*valgrad[k].penalty;

                  if(valgrad[k].property>-INFINITY && valgrad[k].property<INFINITY)
                  {

                     double p=1.0/(1.0+exp(grad[k]/T));
                     double random_nr=(double) random()/(double) RAND_MAX;
                     if(random_nr > p)
                        j+=(long) (p*(double) bases[i]*0.5)%bases[i];
                     else
                        j-=(long) (p*(double) bases[i]*0.5)%bases[i];

                  }

                  number+=j*i;
                  i*=bases[i];
               }
               conf1=number;

               opt_object.lib_object_r.compute_property(conf1);
               config=opt_object.lib_object_r.visited_r.contains(conf1);
               conf3=conf1;
               conf1=opt_object.lib_object_r.deprune(conf1);

               cout << id_r+"::New starting value is: " << opt_object.lib_object_r.value_r[config].property;
               cout << " Penalty: ";
               opt_object.lib_object_r.value_r[config].penalty.display();
               cout << " lambda: ";
               lambda.display(cout);
               cout  << " Result: " << (opt_object.lib_object_r.value_r[config].property-
                     opt_object.lib_object_r.value_r[config].penalty*lambda)
								               << " for compound #" << conf1 << endl;
               cout.flush();

               // Single optimization run done.

            }
            opt_object.lib_object_r.prune(lambda,
                  conf3,
                  conf2,
                  config);

            lambda*=1.1;
            cout << id_r+"::New lambda = " << endl;
            lambda.display(cout);

            steps++;
         }
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
         throw domain_error(id_r+"::gen_base_gdmc::optimize(ulong N) const");
      }
   }

   valerg get_value(ulong i) const
   {
      return opt_object.lib_object_r.get_value(i);
   }

   void set_compute_property_flag(bool B) const {opt_object.set_compute_property_flag(B);}
protected:

private:

};

#endif // _BINARYGDMC_HH_
