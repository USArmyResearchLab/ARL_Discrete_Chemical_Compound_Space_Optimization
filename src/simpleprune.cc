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
//! \file simpleprune.cc \brief Implementation file of simple_prune

#include <BCR_CPP_LA/refcount.h>
#include <simpleprune.h>
#include <sorting_functions.hh>

using namespace linear_algebra;

//! Compute the new lagrange multiplier and prune the library of inferior values.
template <class X>
ulong simple_prune<X>::prune(refvector<double>& oldlambda,
      ulong& conf1,
      ulong& conf2,
      long& config,
      const refvector<ulong>& visited_run) const
      /*!
       * Inferior values have larger penalties than the current best.
       * Lambda is adjusted such that \f$ Prop_{conf1}-\lambda^\prime Pen_{conf1} < Prop_j-\lambda^\prime Pen_j \f$ for some j and \f$\lambda^\prime>\lambda\f$
       */
      {
   try {
      // Lagrange multiplier evaluation.
      long i,j;
      // smallest penalty and largest property.
      long min_max=0;
      ulong conf3=deprune(conf1);
      config=visited.contains(deprune(conf1));
      if(config<0)
         throw domain_error("simple_prune<X>::prune(refvector<double>& lambda,ulong& conf1,ulong& conf2,long& config: conf1 not contained in visited_run) const");
      conf1=visited[config];
      const refvector<double>& newlambda=value_r[config].penalty;
      double lambda=0.0;

      for(i=0;i< visited_run.size() && (
            value[config].penalty*newlambda-
            value[visited_r.contains(visited_run[i])].penalty*newlambda<1e-16 ||
            value[visited_r.contains(visited_run[i])].penalty*newlambda == INFINITY
      )
      ;i++)
      {
         j=visited_r.contains(visited_run[i]);
         // min
         if (value[min_max].penalty*newlambda>=
               value[j].penalty*newlambda)
            // max
            if (!(value[min_max].penalty*newlambda==
                  value[j].penalty*newlambda &&
                  value[min_max].property-value[min_max].penalty*oldlambda>=
                  value[j].property-value[j].penalty*oldlambda)
            )
               min_max=j;
      }

      if(i< visited_run.size())
      {
         j=visited_r.contains(visited_run[i]);
         if(value[config].penalty*newlambda>value[j].penalty*newlambda)
            lambda=
                  fabs((value[config].property-value[config].penalty*oldlambda-
                        value[j].property-value[j].penalty*oldlambda)
                        /
                        ((value[config].penalty-value[j].penalty)*newlambda)
                  );
         else lambda*=1.1; 
         conf1=visited_r[j];
      }
      for(;
            i< visited_run.size(); i++)
      {
         j=visited_r.contains(visited_run[i]);
         // min
         if (value[min_max].penalty*newlambda>=
               value[j].penalty*newlambda)
            // max
            if (!(value[min_max].penalty*newlambda==
                  value[j].penalty*newlambda &&
                  value[min_max].property-value[min_max].penalty*oldlambda>
            value[j].property-value[j].penalty*oldlambda)
            )
               min_max=j;

         if(fabs(value[config].property-value[config].penalty*oldlambda-
               value[j].property-value[j].penalty*oldlambda)
               <
               ((value[config].penalty-value[j].penalty)*newlambda)*lambda
         )
         {
            lambda=
                  fabs((value[config].property-value[config].penalty*oldlambda-
                        value[j].property-value[j].penalty*oldlambda)
                        /
                        ((value[config].penalty-value[j].penalty)*newlambda)
                  );
            conf1=visited_run[i];
         }
         oldlambda+=newlambda*lambda;
      }
      if(deprune(conf2)==conf1 && (visited_r[min_max]!=conf1 || conf3!=conf1)) {
         conf2=visited_r[min_max];
      }
      else conf2=deprune(conf2);

      config=visited_r.contains(conf1);

      // Lambda has been discerned and conf1 is the index in the global library


      // Compute the pruned library

      refvector<ulong> dummy;
      for(i=0;i<visited_run.size();i++)
         if(value[config].penalty*oldlambda<value[visited_r.contains(visited_run[i])].penalty*oldlambda)
            dummy.push_back(visited_run[i]);

      refvector<unsigned long> dummy_index;

      dummy_index=sort_ascending(dummy);

      cout << "pruned indices:\n";
      dummy.display();

      pruned_visited->resize(dummy.size());
      for(i=0;i<dummy.size();i++)
         pruned_visited[i]=dummy[dummy_index[i]];

      // Pruned indices determined
      conf3=conf1;

      // Compute pruned index conf1, conf2
      for(i=pruned_visited.size()-1;i>=0;i--)
         if(conf1>=pruned_visited[i])
            conf1--;
      for(i=pruned_visited.size()-1;i>=0;i--)
         if(conf2>=pruned_visited[i])
            conf2--;

      if(deprune(conf1)!=conf3)
         throw domain_error("prune something wrong in deprune/prune");
      config=visited_r.contains(conf3);

      space_size_computed=false;
      bits_computed=false;

      return conf1;
   } catch(domain_error e) {
      cerr << e.what() << endl;
      throw domain_error("called by ulong simple_prune<X>::prune(double& lambda,\
      ulong& conf1,\
      ulong& conf2,\
      long& config,\
      const refvector<ulong>& visited_run) const");
   }
      }

template <class X>
ulong simple_prune<X>::deprune(ulong N) const
{
   ulong i=N;
   for(long j=0;j<pruned_visited.size();j++)
      if(i>=pruned_visited[j])
      {
         i++;
      }
   return i;
}


template <class X>
ulong simple_prune<X>::reprune(ulong N) const
{
   ulong i;
   ulong conf1=N;
   for(i=pruned_visited.size();i>=1;i--)
      if(conf1>=pruned_visited[i-1])
         conf1--;
   return conf1;
}
