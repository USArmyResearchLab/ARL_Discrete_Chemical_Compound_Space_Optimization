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
//! \file binary_line_search.hh \brief implementation of a line search optimization

#ifndef _BINARY_LINE_SEARCH_HH
#define _BINARY_LINE_SEARCH_HH

#include <optimizeabstract.h>

using namespace std;
using namespace linear_algebra;

//! Search a pruned Library  using a linesearch algorithm using a binary hypercube
template<class C>
class binary_line_search : public optimize_abstract, public C
/*!
   In this method the Library is spanned <EM>via</EM> a hypercube.
   Each molecule in the library sits on one vertex of the cube.
   This leads to a bit representation of each molecule. The line search
   is then performed on the bits.
 */
{
private:

protected:

   using C::value;
   using C::space_size;
public:

   using C::visited_r;
   using C::compute_property;
   using C::prune;
   using C::deprune;
   using C::get_space_size;
   using C::value_r;
   using C::set_number_of_constraints;

   //! Default constructor
   binary_line_search() : optimize_abstract(), C() {};

   //! Copy constructor.
   binary_line_search(const binary_line_search& a) :
      optimize_abstract(a), C(a)
   {};

   using optimize_abstract::optimize;

   //! Default optimization with a clear history and no pruning in the first iteration.
   ulong optimize(ulong N) const
   /*!
      This is a clean slate optimization. Any pruned molecules from previous runs are
      once again viable.
      \param N gives the starting molecule's number in the library.
    */
   {
      refvector<double> l;
      refvector<ulong> visited_run;
      C::pruned_visited.clear();
      return optimize(N,l,visited_run);
   }

   //! Optimize starting from conformation N.
   /*!
     A Lagrange multiplier method is employed to enforce
     boundary constraints. The multiplier increases as the
     optimization proceeds, enforcing the constraint ever
     more rigorously.

   In essence we are following an interior point method, in which we reestablish
   feasibility every log(library size) cycles. The precise rules for pruning
   and adjusting the Lagrange multiplier can be found in prune().
   At each position in the bit string representing the reference molecule N, the molecule
   with that bit flipped from N is computed and compared. If there is an improvement
   in the objective function \f$P+\lambda \pi\f$, where P is the property, \f$\pi\f$
   is the penalty and \f$\lambda\f$ is the current value of the Lagrange multiplier,
   then the new molecule becomes the new reference.
   \see noprune::adjust_lagrange(), simple_prune::prune(), reorder_general_base::prune()
    */
   ulong optimize(ulong N,refvector<double>& lambda, refvector<ulong>& visited_run) const
   {
      try {
         ulong i,j;
         ulong number=C::reprune(N);
         valerg current_best_val;

         current_best_val=compute_property(number);
         j=deprune(number);
         if(
               visited_run.contains(j)<0 &&
               visited_r.contains(j) >= 0
         )
            visited_run.push_back(j);
         long config=visited_r.contains(j);
         if(lambda.size() != current_best_val.penalty.size()) {
            lambda.copy(value_r[config].penalty);
            lambda.zero();
         }
         // Actual optimization routine.
         ulong conf1;
         ulong conf2;
         conf1=number;
         conf2=conf1-1;

         get_space_size();

         valerg interim;

         while(conf1!=conf2)
         {
            conf2=conf1;
            for(i=1;i<space_size;)
            {
               j=((conf1-conf1 % i)/i)%2;
               number=conf1-j*i +((j+1) %2)*i;

               cout << "In "<< id_r << "::optimize(): " << i << endl;

               interim=compute_property(number);
               j=deprune(number);
               if(
                     visited_run.contains(j)<0 &&
                     visited_r.contains(j) >= 0
               )
                  visited_run.push_back(j);
               cout << id_r
                     << "::Config: "
                     << j << "(" << visited_r.contains(j)
                     << ")  finished with property: " << interim.property
                     << " and penalty: ";
               interim.penalty.display();

               if(interim.property-lambda*interim.penalty >
               current_best_val.property-lambda*current_best_val.penalty)
               {
                  cout << " > ";
                  cout << deprune(conf1) << "(" << visited_r.contains(deprune(conf1)) << ") = "
                        << current_best_val.property
                        << " and penalty: ";
                  current_best_val.penalty.display();
                  cout << "\n";

                  current_best_val=interim;
                  conf1=number;
                  config=visited_r.contains(j);
               }
               else
               {
                  config=visited_r.contains(deprune(conf1));
                  cout << " < ";
                  cout << deprune(conf1) << "(" << config << ") = "
                        << current_best_val.property
                        << " and penalty: ";
                  current_best_val.penalty.display();
                  cout << "\n";

               }
               if(i<space_size/2) i*=2;
               else i=space_size;
            }

            cout << "In " << id_r << "::optimize(): "
                  << " Penalty= ";
            value[config].penalty.display();
            cout << " lambda= ";
            lambda.display();
            cout << " Result= " << (value[config].property-
                  value[config].penalty*lambda)
                                << " for " << C::deprune(conf1)
            << endl;
            // Single optimization run done.

            cout << id_r << ": ";
            prune(lambda,
                  conf1,
                  conf2,
                  config,
                  visited_run);

            current_best_val=value[config];
            lambda*=1.1;
            cout << id_r << "::New lambda = ";
            lambda.display();
            get_space_size();
         }
         cout << id_r << "::optimized value is: " << value[config].property;

         cout << " Penalty: ";
         value[config].penalty.display();
         cout << " lambda: ";
         lambda.display();
         cout << " Result: " << (value[config].property-
               value[config].penalty*lambda)
                          << " for " << C::deprune(conf1)
         << endl;

         return deprune(conf1);
      } catch(exception& e) {
         cerr << e.what() << endl;
         throw domain_error("binary_line_search<"+id_r+">::optimize(ulong N) const");
      }
   }

   void set_compute_property_flag(bool B) const {C::set_compute_property_flag(B);}
};
#endif
