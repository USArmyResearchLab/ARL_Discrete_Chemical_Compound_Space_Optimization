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
//! \file sorting_functions.hh Various sorting functions

#ifndef _SORTING_FUNCTIONS__H
#define _SORTING_FUNCTIONS__H

#include <BCR_CPP_LA/refcount.h>
//#include <sorting_functions.hh>

using namespace linear_algebra;

//! Determine whether all components of a are less than or equal to components of b elementwise.
template<class T> static bool lesseq(const refvector<T>& a, const refvector<T>& b, double tol=1e-16)
{
#ifdef DEBUG
   if(a.dim() != b.dim())
      throw domain_error("less for vectors are a mismatch in dimensions");
#endif
   try {
      for(long i=0;i<a.dim();i++)
         if(a[i]>b[i]+tol)
            return false;
   } catch(exception& e) {
      cerr << e.what();
      throw domain_error("called from lesseq");
   }
   return true;
}

//! Do a qsort on E with immediate storage back into index. ensure index has the right dimension. 
template<class T1> static refvector<unsigned long>& sort_ascending(const refvector<T1>& E,
      unsigned long start,
      unsigned long end,
      refvector<unsigned long>& index)
{
   try {
      if(end-start<1) return index;
      unsigned long current=index[start];
      {
         unsigned long less=0;
         unsigned long more=0;
         refvector<unsigned long> split(end-start);
         unsigned long i,l;
         for(i=start+1;i<end;i++)
            if(E[index[i]]<E[index[start]]) {
               split[less]=index[i];
               less++;
            } else {
               split[end-start-more-1]=index[i];
               more++;
            }

         for(i=start, l=0;i<start+less;i++,l++)
            index[i]=split[l];
         index[i]=current;
         current=i;
         i++;
         for(l=0;i<end;i++,l++)
            index[i]=split[less+1+l];
      }
      index=sort_ascending(E,start,current,index);
      return sort_ascending(E,current+1,end,index);
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("sort_ascending");
   }
}

template<class T1> static refvector<unsigned long> sort_ascending(const refvector<T1>& E)
{
   try {
      refvector<unsigned long> index(E.size());
      for(unsigned long i=0;i<(ulong) E.size();i++)
         index[i]=i;
      return sort_ascending(E,(unsigned long) 0,(unsigned long) E.size(),index);
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("sort_ascending");
   }
}


#endif /* _SORTING_FUNCTIONS__H */


