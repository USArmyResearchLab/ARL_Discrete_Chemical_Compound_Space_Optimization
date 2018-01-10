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
//! \file binarygdmc.hh \brief Implementation of binary GDMC

#ifndef _BINARYGDMC_HH_
#define _BINARYGDMC_HH_

#include <typedefs.hh>
#include <iostream>
#include <optimizeabstract.h>
#include <stdexcept>

using namespace std;

//! Gradient-directed Monte Carlo performed on a hypercube
template <class C>
class binary_gdmc: public optimize_abstract
/*!
   This class takes a Library object and performs an optimization.
   The Library is arranged in a hypercube topology in which each molecule occupies
   a vertex on the hull of the cube. After a local optimum has been reached, a
   Monte-Carlo simulation is performed to determine a new starting point.
   Each time a new minimum is found the acceptance criteria become more strict.
 */
{
private:

   //! This is the library object to be used in the optimization.
   C opt_object;

public:
   //! Temperature-analogue to be used for Monte-Carlo evaluations.
   double T;
   //! Number of tightening steps. Determines the strictness of acceptance.
   ulong tight_steps;
   //! Maximum number of steps.
   ulong max_steps;

   //! Copy Constructor.
   binary_gdmc(const C& a) : opt_object(a), T(0.0), tight_steps(1), max_steps(1) {};

   //! Gradient-directed Monte-Carlo optimization
   ulong optimize(ulong N) const
   /*!
      The gradient to be used is determined by the gradient() method of the class C.
    */
   {
      try {
         ulong i,j,k;
         ulong number=N;
         long config=opt_object.visited_r.contains(N);
         refvector<double> lambda(opt_object.get_number_of_constraints());
         ulong steps=1;

         opt_object.set_id(id_r+"::opt_object");

         // Actual optimization routine.
         ulong conf1,conf2;
         conf1=number;

         refvector<valerg> valgrad(opt_object.get_bits());

         refvector<double> grad(valgrad.dim());

         ulong conf3;
         // tight_steps counts how often lambda is ramped.
         while(tight_steps>=steps)
         {
            // To force lambda ramping the average number of computed values is limited.
            while(steps*max_steps/tight_steps>(ulong) opt_object.visited_r.size())
            {
               conf2=conf1;
               //This reinitializes the optimization intermediates.
               conf1=opt_object.optimize(conf1);
               cout << id_r+"::Gradient of " << conf1 << endl;
               cout.flush();
               opt_object.get_space_size();
               opt_object.get_bits();
               conf1=opt_object.reprune(conf1);
               valgrad=opt_object.gradient(conf1,valgrad);
               number=0;
               ulong space_size=opt_object.get_space_size();
               for(k=0,i=1;i<space_size;k++)
               {
                  if(valgrad[k].property>-INFINITY && valgrad[k].property<INFINITY)
                  {
                     grad[k]=valgrad[k].property-lambda*valgrad[k].penalty;

                     double p=1.0/(1.0+exp(grad[k]/T));
                     double random_nr=(double) random()/(double) RAND_MAX;
                     if(random_nr > p)
                        j=1;
                     else j=0;
                     if(number+j*i<space_size)
                        number+=j*i;

                  }
                  else
                     if(number+((conf1 % i)%2)*i<space_size)
                        number+=((conf1 % i)%2)*i;

                  if(i<space_size/2) i*=2;
                  else i=space_size;
               }
               conf1=number;

               opt_object.compute_property(conf1);
               config=opt_object.visited_r.contains(opt_object.deprune(conf1));
               conf3=conf1;
               conf1=opt_object.deprune(conf1);

               cout << id_r+"::New starting value is: " << opt_object.value_r[config].property;
               cout << " Penalty: ";
               opt_object.value_r[config].penalty.display();
               cout << " lambda: ";
               lambda.display();
               cout << " Result: " << (opt_object.value_r[config].property-
                     opt_object.value_r[config].penalty*lambda)
                                      << " for compound #" << conf1 << endl;
               cout.flush();

               // Single optimization run done.

            }
            opt_object.prune(lambda,
                  conf3,
                  conf2,
                  config);

            lambda*=1.1;
            cout << id_r+"::New lambda = ";
            lambda.display();

            steps++;
         }
         // minimax
         conf1=0;
         for(i=0;i<(ulong) opt_object.visited_r.size();i++)
         {
            if(opt_object.value_r[i].penalty*opt_object.value_r[i].penalty<=
                  opt_object.value_r[conf1].penalty*opt_object.value_r[conf1].penalty)
               if(!(opt_object.value_r[i].penalty*opt_object.value_r[i].penalty==opt_object.value_r[conf1].penalty*opt_object.value_r[conf1].penalty) ||
                     opt_object.value_r[i].property > opt_object.value_r[conf1].property)
                  conf1=i;
         }
         conf1=opt_object.visited_r[conf1];
         return conf1;
      } catch(exception& e) {
         cerr << e.what() << endl;
         throw domain_error(id_r+"::binary_gdmc::optimize(ulong N) const");
      }
   }

   //! Return the value for molecule i.
   valerg get_value(const ulong i) const
   {
      return opt_object.get_value(i);
   }

   void set_compute_property_flag(bool B) const {opt_object.set_compute_property_flag(B);}
protected:

private:

};

#endif // _BINARYGDMC_HH_
