/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*- */
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
//! \file noprune.h \brief Header of class noprune

#ifndef _NOPRUNE_H_
#define _NOPRUNE_H_

#include <prunerabstract.h>
#include <sorting_functions.hh>

//! Dummy pruning class. No actual pruning done.
template<class X>
class noprune : public pruner_abstract<X>
{
public:

   using pruner_abstract<X>::visited_r;
   using pruner_abstract<X>::value_r;
   using pruner_abstract<X>::prune;
   using pruner_abstract<X>::get_badval;

   bool minimax;

   noprune() : pruner_abstract<X>(), minimax(false) {};

   noprune(const noprune<X>& a) :
      pruner_abstract<X>(a), minimax(a.minimax)
      {};

   noprune<X>& operator=(const noprune<X>& A)
   {
      (pruner_abstract<X>&) *this=(pruner_abstract<X>&) A;
      minimax = A.minimax;
      return *this;
   }
   //! Compute the increase of lambda and assess current best value.
   void adjust_lagrange(refvector<double> &oldlambda, ulong & conf1, ulong & conf2, long &config, const refvector<ulong>& visited_run) const
   {
      config=visited_r.contains(deprune(conf1));
      const ulong conf3=deprune(conf1); // save original
      string serr="reorder_general_base<X>::adjust_lagrange(double& lambda,ulong& conf1,ulong& conf2,long& config)";
      if(config<0)
         throw domain_error(serr+": conf1 not contained in visited) const");
      long min_max=0;
      //NOTE: We assume penalty >= 0.0
      const refvector<double>& newlambda=value_r[config].penalty;
      refvector<double> newlambda2=value_r[config].penalty;
      double lambda=0.0;
      if(!minimax) {
         // Lagrange multiplier evaluation.
         // minimize over all positive lagrange estimators
         // (p_0-p_i)/(pi_0-pi_i) >= 0, i>0
         long i,j;
         // smallest penalty and largest property.
         try {
            conf1=visited_r[config];

            for(i=0;i< visited_run.size() &&
            (
                  (value_r[config].penalty-
                        value_r[visited_r.contains(visited_run[i])].penalty)*newlambda<1e-16 ||
                        value_r[visited_r.contains(visited_run[i])].penalty == get_badval().penalty
            )
            ;i++)
            {
               j=visited_r.contains(visited_run[i]);
               // min
               if (value_r[min_max].penalty*newlambda>=
                     value_r[j].penalty*newlambda)
                  // max
                  if (!(value_r[min_max].penalty*newlambda==
                        value_r[j].penalty*newlambda &&
                        value_r[min_max].property-value_r[min_max].penalty*oldlambda>=
                        value_r[j].property-value_r[j].penalty*oldlambda)
                  )
                     min_max=j;
            }
            if(i< visited_run.size())
            {
               j=visited_r.contains(visited_run[i]);
               if(value_r[config].penalty*newlambda>value_r[j].penalty*newlambda)
                  lambda=
                        fabs((value_r[config].property-value_r[config].penalty*oldlambda-
                              (value_r[j].property-value_r[j].penalty*oldlambda))
                              /
                              ((value_r[config].penalty-value_r[j].penalty)*newlambda)
                        );
               else lambda*=1.1;
               conf1=visited_r[j];
               i++;
            }
            for(; i< visited_run.size(); i++) {
               j=visited_r.contains(visited_run[i]);
               // min
               if (value_r[min_max].penalty*newlambda>=
                     value_r[j].penalty*newlambda)
                  // max
                  if (!(value_r[min_max].penalty*newlambda==
                        value_r[j].penalty*newlambda &&
                        value_r[min_max].property-value_r[min_max].penalty*oldlambda>
                  value_r[j].property-value_r[j].penalty*oldlambda)
                  )
                     min_max=j;

               if(fabs(value_r[config].property-value_r[config].penalty*oldlambda-value_r[j].property+value_r[j].penalty*oldlambda)
                     <
                     ((value_r[config].penalty-value_r[j].penalty)*newlambda)*lambda
               )
               {
                  lambda=
                        fabs(
                              (value_r[config].property-value_r[config].penalty*oldlambda-value_r[j].property+value_r[j].penalty*oldlambda)
                              /
                              ((value_r[config].penalty-value_r[j].penalty)*newlambda)
                        );
                  conf1=visited_r[j];
               }
            }
            oldlambda+=newlambda*lambda;
         } catch(domain_error &e) {
            cerr << e.what() << endl;
            throw domain_error("called by "+serr);
         }
      }// end if (!minimax)
      else {
         long i,j;
         conf1=conf3;
         refvector<double> min=value_r[config].penalty;
         double max=value_r[config].property;
         long newconf=conf1;
         min_max=config;

         //minimize the penalty
         for(i=0;i<visited_run.size();i++)
         {
            j=visited_r.contains(visited_run[i]);
            // If at least one component is lower
            if(lesseq<double>(value_r[j].penalty,min))
               // Only if ALL components are lower or the property value increases
               if(!lesseq<double>(min,value_r[j].penalty) || value_r[j].property>max) {
                  min.copy(value_r[j].penalty);
                  max = value_r[j].property;
                  conf1 = visited_r[j];
                  min_max=j;
               }
         }
         // adjust lambda such that min_max is optimal
         /*! We demand that \f$ P(x)-\lambda^T\pi(x)\ge P(y)-\lambda^T\pi(y)\forall y\f$.
          * We let \f$\lambda_{n+1}=\lambda_n+\lambda^\prime \Delta \lambda_n\f$.
          * Hence
          * \f[
          * P(x)-P(y)-\lambda_n^T\Delta\pi\ge \lambda^\prime\Delta\lambda_n^T\Delta \pi
          * \f]
          * Since \f$y\f$ is the final result of a run, the lhs is less than or equal to 0.
          * Let \f$\Delta\lambda_n=|-\Delta \pi|\f$, where \f$|f|\f$ sets negative components to zero.
          * Then \f$\lambda^\prime = \frac{P(x)-P(y)-\lambda_n^T\Delta\pi}{|-\Delta \pi|^T\Delta \pi}\f$.
          */
         newlambda2-=min;
         if((max-value_r[config].property)-(min*oldlambda-value_r[config].penalty*oldlambda)<-1e-16)
         {
            for(i=0;i<newlambda2.dim();i++)
               if(newlambda2[i]<0.0) newlambda2[i]=0.0;
            lambda=(max-value_r[config].property-(min*oldlambda-value_r[config].penalty*oldlambda))/(min*newlambda2-value_r[config].penalty*newlambda2);
         }
         oldlambda+=newlambda2*lambda;
      }
      // Ensure that conf1 and conf2 point to reasonable answers.
      // If the previous iteration's result has not changed
      if(deprune(conf2)==conf1) {
         // but the best guess did change during optimization, but before readjustment
         if(conf3!=conf1) conf2=conf3; // set previous iteration's result to interim position
         else conf2=visited_r[min_max]; // otherwise set to minimax
      }
      else conf2=deprune(conf2);//keep conf2 as before

      config=visited_r.contains(conf1);
   }


   ulong prune(refvector<double>& lambda, ulong &conf1, ulong &conf2, long &config, const refvector<ulong>& a) const {
      adjust_lagrange(lambda, conf1, conf2, config, a);
      return conf1;
   }

   ulong deprune(ulong N) const { return N; }
   ulong reprune(ulong N) const { return N;}
protected:
private:

};

#endif // _NOPRUNE_H_
