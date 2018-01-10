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

#ifndef _REORDERGENERALBASE_HH_
#define _REORDERGENERALBASE_HH_

#include <prunerabstract.h>
#include <typedefs.hh>
#include <BCR_CPP_LA/refcount.h>
#include <sorting_functions.hh>

/*!
This class takes an array of bases and orders the bases so as
to produce a smooth optimization process. Each digit per base
is weighted by the average computed values for compounds with
the appropriate digit. There is no actual pruning done!
 */
template <class X>
class reorder_general_base: public pruner_abstract<X>
{
   //! Records the order of digits in each base.
   mutable refvector < refvector < long > > base_order;
   //! Records the average values associated with each digit and base.
   mutable refvector < refvector<valerg> > base_averages;
   //! Records the size of visited prior to the latest iteration.
   mutable ulong oldvisitedsize;

public:

   using pruner_abstract<X>::visited_r;
   using pruner_abstract<X>::value_r;
   using pruner_abstract<X>::prune;
   using pruner_abstract<X>::get_badval;

   bool minimax;
   bool at_max;
   bool at_current;

   //! Copy Constructor
   reorder_general_base(reorder_general_base<X>& L):
      pruner_abstract<X>(L),
      base_order(L.base_order),
      base_averages(L.base_averages),
      oldvisitedsize(L.oldvisitedsize),
      minimax(L.minimax),
      at_max(L.at_max),
      at_current(L.at_current)
      {};
   //! Copy Constructor
   reorder_general_base(const reorder_general_base<X>& L):
      pruner_abstract<X>(L),
      base_order(L.base_order),
      base_averages(L.base_averages),
      oldvisitedsize(L.oldvisitedsize),
      minimax(L.minimax),
      at_max(L.at_max),
      at_current(L.at_current)
      {};

   //! Constructor with a vector of bases.
   reorder_general_base(const refvector < long >& ibases):
      pruner_abstract<X>(),
      base_order(ibases.size ()),
      base_averages(ibases.size()),
      oldvisitedsize(0),
      minimax(false),
      at_max(false),
      at_current(false)
      {
      for (long i = 0; i < ibases.size (); i++)
      {
         refvector<long> r(ibases[i]);
         refvector<valerg> v(ibases[i]);
         base_averages[i]=v;
         for(long j=0;j<ibases[i];j++)
            r[j]=j;
         base_order[i]=r;
      }
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

   //! Adjust the Lagrange multiplier and prune the library.
   ulong prune(refvector<double> &lambda, ulong & conf1, ulong & conf2, long &config, const refvector<ulong>& visited_run) const
   {
      string serr="reorder_general_base<X>::prune(refvector<double> &lambda, ulong & conf1, ulong & conf2, long &config, const refvector<ulong>& visited_run) const";

      refvector<refvector<double> > val_averages;
      try {
         // conf1 and conf2 go in pruned and come out depruned
         adjust_lagrange(lambda, conf1, conf2, config, visited_run);

         // Lambda has been discerned and conf1 is the index in the global library


         // Compute the orders.
         long i,j;

         // first update averages

         refvector< refvector<valerg> > new_averages;
         new_averages=(const refvector< refvector<valerg> >) base_averages;
         refvector<refvector<ulong> > base_visited;
         for(i=0; i<base_averages.size();i++)
         {
            refvector<double> v(base_averages[i].size());
            val_averages.push_back(v);
            refvector<ulong> a(base_averages[i].size());
            base_visited.push_back(a);
         }

         for(i=0;i<visited_run.size();i++)
         {
            ulong N=visited_run[i];
            long m=visited_r.contains(N);
            for(j=0;j<base_order.size();j++)
            {
               long k=N % base_order[j].size();
               N=(N-k)/base_order[j].size();
               bool infinities=(value_r[m].property == INFINITY ||
                                value_r[m].property == -INFINITY);
               for(long dim=0;dim<value_r[m].penalty.size();dim++)
                  infinities = infinities ||
                              value_r[m].penalty[dim] == INFINITY ||
                              value_r[m].penalty[dim] == -INFINITY;
               if(infinities)
                  val_averages[j][k]-=M_PI*0.5;
               else
                  val_averages[j][k]+=atan(value_r[m].property-lambda*value_r[m].penalty);
               base_visited[j][k]++;
            }
         }
         long N=conf1;
         for(i=0;i<val_averages.size();i++) {
            long k=N % base_order[i].size();
            N=(N-k)/base_order[i].size();
            for(j=0;j<val_averages[i].size();j++)
               if(base_visited[i][j]==0) {
                  // this ensures that unvisited positions are injected into the search path
                  cout << serr << ": " << i << "," << j << " has no representative."
                        << " Current best bit: " << k << endl;
                  if(at_current)
                     val_averages[i][j]=val_averages[i][k];
                  if(at_max)
                     val_averages[i][j]=0.5*M_PI*(double) base_visited[i][k];
                  // otherwise it stays at zero.
                  base_visited[i][j]=base_visited[i][k];
               }
            for(j=0;j<val_averages[i].size();j++)
               // average value a certain base digit has
               val_averages[i][j]*=1.0/(double) base_visited[i][j];
         }

         new_averages*=
               (double) 1/(double) (visited_run.size() - oldvisitedsize);
         for(i=0;i<base_order.size();i++)
         {
            double average=0.0;
            refvector<ulong> index;
            index=sort_ascending(val_averages[i]);
            for(j=0;j<index.size()-1;j+=2) {
               base_order[i][j/2]=index[j];
               base_order[i][index.size()-1-j/2]=index[j+1];
            }
            if(j<index.size())
               base_order[i][j/2]=index[j];
         }
         oldvisitedsize=visited_run.size();

         config=visited_r.contains(conf1);
         conf1=reprune(conf1);
         conf2=reprune(conf2);

      } catch(domain_error &e) {
         cerr << e.what() << endl;
         throw domain_error("called by "+serr);
      }
#ifdef DEBUG
      cout << "The new order is:\n";
      for (long i=0;i<base_order.size();i++)
         for(long j=0;j<base_order[i].size();j++)
            cout << "Base " << i << ": Substituent " << j << " -> " << base_order[i][j]
                 << " with average "
                 << val_averages[i][base_order[i][j]]
                 << endl;
#endif
      return conf1;
   }

   //! Reverse-map N to the global index.
   ulong deprune (ulong N) const
   {
      long i;
      ulong j;
      ulong M=0;
      j=1;
      for(i=0;i<base_order.size();i++)
      {
         ulong k;
         k=N % base_order[i].size();
         N=(N-k)/base_order[i].size();
         M+=(ulong) base_order[i][k]*j;
         j*=(ulong) base_order[i].size();
      }
      return M;
   }

   ulong reprune(ulong N) const
   {
      ulong M=0;
      ulong k=1;
      ulong conf1=N;
      long i,j;
      for(i=0;i<base_order.size();i++)
      {
         j=conf1 % base_order[i].size();
         conf1=(conf1-j)/base_order[i].size();
         M+=base_order[i].contains(j)*k;
         k*=(ulong) base_order[i].size();
      }
      conf1=M;
      return conf1;
   }
protected:

private:

};

#endif // _REORDERGENERALBASE_HH_
