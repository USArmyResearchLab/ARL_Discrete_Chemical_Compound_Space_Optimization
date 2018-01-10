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
//! \file binarysteepestdescent.hh \brief Implements the binary_steepest_descent optimization class.

#ifndef _BINARYSTEEPESTDESCENT_HH_
#define _BINARYSTEEPESTDESCENT_HH_

#include <optimizeabstract.h>
#include <BCR_CPP_LA/refcount.h>
#include <has_gradients_hessian_data.hh>

using namespace std;

//! Searches a Library for an optimum using a steepest descent method.
template <class C>
class binary_steepest_descent: public optimize_abstract, public C
/*!
   As all optimization classes in the binary_* vein, the library is organized and
   searched as a hypercube on which the molecules are the vertices on the hull.
   The resultant bit-string representation of each molecule is used to perform
   optimizations.
 */
{
private:
protected:
   using C::space_size;
public:

   using C::deprune;
   using C::get_space_size;
   using C::get_bits;
   using optimize_abstract::optimize;
   using C::gradient;
   using C::set_compute_property_flag;
   using C::compute_property;
   using C::prune;

   binary_steepest_descent():
      optimize_abstract(),
      C()
   {};

   ~binary_steepest_descent() {};

   //! Optimize starting from index N.
   ulong optimize(ulong N) const
   /*!
      Each edge (i.e., bit) on the hypercube is interpreted as a search direction.
      A gradient is computed along each direction and is then searched for improvements.
      Constraints are enforced with Lagrange multipliers which are ramped after a
      local optimum is reached.
    */
   {
      try {
         ulong i,j,k;
         ulong number=N;
         valerg current_best_val;
         long config=C::visited.contains(deprune(N));
         refvector<double> lambda(C::get_number_of_constraints());

         current_best_val=compute_property(number);

         // Actual optimization routine.
         ulong conf1;
         ulong conf2;
         conf1=number;
         conf2=conf1-1;

         get_space_size();

         valerg interim;
         refvector<valerg> valgrad;
         cout << " initial gradient \n"; cout.flush();
         valgrad=gradient(conf1);
         cout << "gradient done \n"; cout.flush();
         refvector<double> grad(valgrad.dim());

         while(conf1!=conf2) {
            while(conf1!=conf2)
            {
               conf2=conf1;
               cout << "Gradient of " << conf1 << endl;
               cout.flush();
               valgrad=gradient(conf1);
               number=0;
               for(k=0,i=1;i<space_size;k++) {
                  if(valgrad[k].property>-INFINITY &&
                        valgrad[k].property<INFINITY)
                  {
                     grad[k]=valgrad[k].property-lambda*valgrad[k].penalty;
                     if(grad[k]>0)
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

               interim=compute_property(number);
               while(interim.property==-INFINITY) {
                  number++;
                  interim=compute_property(number);
               }

               conf1=number;
               config=C::visited.contains(deprune(conf1));
               cout << "Config: "
                     << number << "(" << config
                     << ")  finished with property: " << interim.property
                     << " and penalty: ";
               interim.penalty.display();
            }
            config=C::visited.contains(deprune(conf1));
            cout << "The optimized value is: " << C::value[config].property << endl;
            cout << " Penalty: ";
            C::value[config].penalty.display();
            cout << " lambda: ";
            lambda.display();
            cout << " Result: " << (C::value[config].property-
                  C::value[config].penalty*lambda)
					                  << endl;

            // Single optimization run done.

            prune(lambda,
                  conf1,
                  conf2,
                  config);

            config=C::visited.contains(deprune(conf1));
            lambda*=1.1;
            cout << " New lambda = ";
            lambda.display();
         }
         return conf1;
      } catch(exception& e) {
         cerr << e.what() << endl;
         throw domain_error("pruned_constrained_chem_opt_sd::optimize_steepest_descent(ulong N) const");
      }
   }
   //void set_compute_property_flag(bool f) {};
   void set_compute_property_flag(bool B) const {C::set_compute_property_flag(B);}
protected:

private:

};

#endif // _BINARYSTEEPESTDESCENT_HH_
