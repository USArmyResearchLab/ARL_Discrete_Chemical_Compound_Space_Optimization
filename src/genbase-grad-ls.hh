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
//! \file genbase-grad-ls.hh Line search optimization with general bases.

#ifndef _GENBASE_GRAD_LS_HH_
#define _GENBASE_GRAD_LS_HH_

#include <typedefs.hh>
#include <iostream>
#include <optimizeabstract.h>
#include <BCR_CPP_LA/refcount.h>

using namespace std;
using namespace linear_algebra;

//! Line search optimization class using general bases.
template <class C, class B>
class gen_base_grad_LS: public optimize_abstract
{
private:

   const C lib_object;
   mutable B bases;

   //! Precondition the library to induce a proper ordering.
   void precondition(ulong &conf1,refvector<ulong>& visited_run) const
   {
      long j;
      long nm;
      valerg interim;
      for(bases=0;!bases.done();bases++)
      {
         cout << "In "<< id_r << "::precondition() " << bases.get_state() << endl;

         ulong conf3=conf1-(((conf1-conf1%bases())/bases()) % bases.modulus()) * bases();

         for(j=0;j<bases.modulus();j++)
         {
            nm=conf3+j*bases();
            interim=lib_object.compute_property(nm);
            update_visited_run(visited_run, nm, interim);
            cout << endl;
         }
      }
      ulong conf2=conf1;
      long config=0;
      refvector<double> l(lib_object.get_number_of_constraints());
      cout << "In "<< id_r << "::precondition() " << "Started with: " << conf1 << endl;
      nm=lib_object.deprune(conf1);
      lib_object.prune(l,
            conf1,
            conf2,
            config,
            visited_run);
      conf1=lib_object.reprune(nm);
      cout << "Ended with " << conf1 << endl;
      return;
   }

   void select_current_best(const valerg& interim, const refvector<double>& lambda, valerg& current_best_val, ulong& conf1, ulong np, long& config) const
   {
      if ((interim.property_computed && interim.property - lambda * interim.penalty > current_best_val.property - lambda * current_best_val.penalty) ||
          lib_object_r.is_badval(current_best_val) ) {
         cout << " > ";
         cout << lib_object.deprune(conf1) << "(" << lib_object.visited_r.contains(lib_object.deprune(conf1)) << ") = " << current_best_val.property << " and penalty: ";
         current_best_val.penalty.display();
         cout << endl;
         conf1 = np;
         current_best_val = interim;
         config = lib_object.visited_r.contains(lib_object.deprune(conf1));
      } else {
         config = lib_object.visited_r.contains(lib_object.deprune(conf1));
         cout << " <= ";
         cout << lib_object.deprune(conf1) << "(" << config << ") = " << current_best_val.property << " and penalty: ";
         current_best_val.penalty.display();
         cout << endl;
      }
   }

   void update_visited_run(refvector<ulong>& visited_run, ulong np, const valerg& interim) const
   {
      if (visited_run.contains(lib_object.deprune(np)) < 0)
         visited_run.push_back(lib_object.deprune(np));

      cout << id_r << "::Config: " << lib_object.deprune(np) << "(" << lib_object.visited_r.contains(lib_object.deprune(np)) << ")  finished with property: " << interim.property << " and penalty: ";
      interim.penalty.display();
   }

   void sweep_direction(ulong &conf1, ulong conf3, refvector<ulong> visited_run, const refvector<double>& lambda, valerg &current_best_val, long &config) const
   {
      ulong np;
      ulong nm;
      valerg interimp;
      valerg interimm;
      valerg old;
      long j;
      long dumbcounter = 0;
      while (conf1 != conf3 && dumbcounter < bases.modulus()) {
         conf3 = conf1;
         dumbcounter++;
         old=current_best_val;
         j = ((conf1 - conf1 % bases()) / bases()) % bases.modulus();
         cout << "In " << id_r << "::optimize(): " << bases.get_state() << endl;
         np = conf1 + (-j + (j + 1) % bases.modulus()) * bases();
         nm = conf1 + (-j + (j + bases.modulus() - 1) % bases.modulus()) * bases();
         interimp = lib_object.compute_property(np);
         update_visited_run(visited_run, np, interimp);
         select_current_best(interimp, lambda, current_best_val, conf1, np, config);
         interimm = lib_object.compute_property(nm);
         update_visited_run(visited_run, nm, interimm);
         select_current_best(interimm, lambda, current_best_val, conf1, nm, config);
         refvector<double> l=(interimp.penalty-interimm.penalty);
         l*=(interimp.property-interimm.property);
         l-=old.penalty*(4.0*(interimp.property+interimm.property-2.0*old.property-
               lambda*(interimp.penalty+interimm.penalty-old.penalty*2.0)));
         l*=1.0/((interimp.penalty-interimm.penalty)*(interimp.penalty-interimm.penalty));
         cout << "In " << id_r << "::optimize():lambda*: ";
         l.display();
         cout << endl;

      }
   }

public:

   bool precondition_flag;
   using optimize_abstract::optimize;

   const B& bases_r;
   const C& lib_object_r;

   gen_base_grad_LS(const C& Library):
      lib_object(Library),
      bases(lib_object),
      precondition_flag(false),
      bases_r(bases),
      lib_object_r(lib_object)
   {};


   //! Copy constructor
   gen_base_grad_LS(const gen_base_grad_LS<C,B>& a):
      lib_object(a.lib_object),
      bases(a.bases),
      precondition_flag(a.precondition_flag),
      bases_r(bases),
      lib_object_r(lib_object)
   {};

   //! Gradient computation
   refvector<valerg> gradient(const ulong conf1) const
	            {
      refvector<valerg> r(bases.non_empty_size());
      return gradient(conf1,r);
	            }
   //! Gradient computation
   refvector<valerg>& gradient(const ulong conf1, refvector<valerg>& r) const
	            {
      long i,j;
      const long saved_state=bases.get_state();
      const ulong saved_refstate=bases.get_refstate();
      bases.set_refstate(conf1);

      ulong np;
      ulong nm;

      for(bases=0,i=0;!bases.done();bases++,i++)
      {
         j=((conf1 - conf1 % bases())/bases()) % bases.modulus();

         np=conf1+(-j+(j+1) % bases.modulus())*bases();
         nm=conf1+(-j+(j+bases.modulus()-1) % bases.modulus())*bases();

         r[i]=lib_object.compute_property(np);
         r[i]-=lib_object.compute_property(nm);
      }
      bases.set_refstate(saved_refstate);
      bases=saved_state;
      return r;
	            }

   long stacksize() const
   {
      return lib_object_r.visited_r.size();
   }

   //! Optimize starting from conformation i.
   /*!
   A Lagrange multiplier method is employed to enforce
   boundary constraints. The multiplier increases as the
   optimization proceeds, enforcing the constraint ever
   more rigorously.
    */
   ulong optimize(ulong N) const
   {
      try
      {
         ulong number=lib_object_r.reprune(N);
         valerg current_best_val;
         long config=lib_object.visited_r.contains(N);
         refvector<ulong> visited_run;

         current_best_val=lib_object.compute_property(number);
         visited_run.push_back(lib_object.deprune(number));
         while(lib_object_r.is_badval(current_best_val) && number<lib_object_r.get_space_size()-1)
         {
            number++;
            current_best_val=lib_object.compute_property(number);
            ulong renumb=lib_object.deprune(number);
            if(visited_run.contains(renumb)<0)
                visited_run.push_back(renumb);
         }
         
         if(precondition_flag)
            precondition(number,visited_run);
         refvector<double> lambda(lib_object_r.get_number_of_constraints());

         // Actual optimization routine.
         ulong conf1;
         ulong conf2;
         conf1=number;
         conf2=conf1-1;
         ulong cycle=0;

         while (conf1!=conf2)
         {
            conf2=conf1;
            bases.set_refstate(lib_object.deprune(conf1));
            cout << id_r << " Starting cycle " << cycle << endl;

            for(bases=0;!bases.done();bases++)
            {
               ulong conf3=conf1-1;
               sweep_direction(conf1, conf3, visited_run, lambda, current_best_val, config);
               bases.set_refstate(lib_object.deprune(conf1));
            }

            config=lib_object.visited_r.contains(lib_object.deprune(conf1));
            {
               cout << id_r+"::optimized value in cycle " << cycle << " is: " << lib_object.value_r[config].property;
               cout << " Penalty: ";
               lib_object.value_r[config].penalty.display();
               cout << " lambda: ";
               lambda.display();
               cout << " Result: " << (lib_object.value_r[config].property-
                     lib_object.value_r[config].penalty*lambda)
									            << " for compound #" << lib_object.deprune(conf1) << endl;
               cout.flush();
               cycle++;
            }
            // Single optimization run done.

            lib_object.prune(lambda,
                  conf1,
                  conf2,
                  config,
                  visited_run);

            current_best_val=lib_object.value_r[config];
            lambda*=1.1;
            cout << id_r << "::New lambda = ";
            lambda.display();
            cout << endl;
            lib_object.get_space_size();
         }
         cout << id_r << "::Visited Run =  ";
         visited_run.display();
         cout << " Number of compounds = " << visited_run.dim() << " in " << cycle << " cycles" << endl;
         return lib_object.deprune(conf1);
      }
      catch (domain_error e)
      {
         cerr << e.what() << endl;
         throw domain_error("called by gen_base_grad_LS::optimize(ulong N) const");
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

#endif // _GENBASE_GRAD_LS_HH_
