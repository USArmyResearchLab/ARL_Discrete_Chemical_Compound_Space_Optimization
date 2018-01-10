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
//! \file generalbaseiterator.cc

#include <chem_opt.hh>
#include <chemident.hh>
#include <generalbaseiterator.hh>

using namespace linear_algebra;

//! Set occupations according to number.
template<>
void general_base_iterator<chem_opt>::occupy(ulong number)
{
   long i,j,m, groupsize,g;

   try {
      g=0;
      for(i=0;i< lib_object.Substituent_Groups_r.size();i++) {
         const ChemIdent& Group((lib_object.Substituent_Groups_r[i]));
         for(j=0;j< Group.allowed_Substituents_r.size();j++,g++) {
            groupsize=Group.allowed_Substituents_r[j].size();
            m=number % groupsize;
            number=(number -m)/ groupsize;
            occupation[g]=m;
         }
      }
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::occupy");
   }
   return;
}
//! Set occupations according to number.
template<>
void general_base_iterator<zmat_opt>::occupy(ulong number)
//! TODO: Implement non-trivial version
{
   return;
}

//! Compute bases based on group Group
template<>
void general_base_iterator<chem_opt>::compute_bases(long Group, const refvector<long>& base_sizes)
{
   long i,M;
   for(i=0;i<lib_object.Substituent_Groups_r[Group].Connector_r.size();i++)
   {
      M=occupation[base_sizes[Group]+i];
      if(bases[base_sizes[Group]+i]==0) {
         bases[base_sizes[Group]+i]=1;
         compute_bases(lib_object.
               Substituent_Groups_r[Group].
               allowed_Substituents_r[i][M],
               base_sizes);
      }
   }
}

template<>
void general_base_iterator<chem_opt>::compute_bases()
{
   long i,j;
   bases.zero();
   long k=1;
   long m=0;

   refvector<long> base_sizes(lib_object.Substituent_Groups_r.size());
   base_sizes[0]=0;

   for(i=1;i<base_sizes.size();i++)
      base_sizes[i]=lib_object.Substituent_Groups_r[i-1].allowed_Substituents_r.size()+base_sizes[i-1];
   number_of_bases=0;
   compute_bases(0,base_sizes);
   for(i=0; i< lib_object.Substituent_Groups_r.size(); i++)
      for(j=0; j<lib_object.Substituent_Groups_r[i].allowed_Substituents_r.size(); j++,m++)
      {
         if(bases[m]==1) {
            bases[m]=k;
            number_of_bases++;
         }
         if(lib_object.Substituent_Groups_r[i].allowed_Substituents_r[j].size()!=0)
            k*=lib_object.Substituent_Groups_r[i].allowed_Substituents_r[j].size();
      }
   return;
}

//! Copy Constructor.
template<>
general_base_iterator<chem_opt>::general_base_iterator(const general_base_iterator<chem_opt>& b):
number_of_bases(b.number_of_bases),
refstate(b.refstate),
state(b.state),
lib_object(b.lib_object),
bases(b.bases),
moduli(b.moduli),
occupation(b.occupation)
{};

//! Constructor.  Mandatory library object.
template<>
general_base_iterator<chem_opt>::general_base_iterator(const chem_opt& b):
lib_object(b)
{
   state=0;
   refstate=0;
   long i,j,m;
   m=0;
   for(i=0;i<b.Substituent_Groups_r.size();i++)
      for(j=0;j<b.Substituent_Groups_r[i].allowed_Substituents_r.size();j++)
         m++;
   refvector<long> bases1(m);
   refvector<long> moduli1(m);
   refvector<long> occupation1(m);
   bases=bases1;
   moduli=moduli1;
   m=0;
   for(i=0;i<b.Substituent_Groups_r.size();i++)
      for(j=0;j<b.Substituent_Groups_r[i].allowed_Substituents_r.size();j++,m++)
         moduli[m]=b.Substituent_Groups_r[i].allowed_Substituents_r[j].size();
   occupation=occupation1;
   occupy(0);
   compute_bases();
}

template<>
long general_base_iterator<chem_opt>::modulus() const
{
   if(state<moduli.size())
      return moduli[state];
   else
      throw domain_error("general_base_iterator<chem_opt>::modulus() const:state out of range");
}

template<>
general_base_iterator<chem_opt>& general_base_iterator<chem_opt>::operator++(int i)
{
   if(!done()) {
      state++;
      while(state<bases.size() &&
            bases[state]==0)
         state++;
   }
   return *this;
}

template<>
ulong general_base_iterator<chem_opt>::get_refstate() const
{
   return refstate;
}

template<>
ulong general_base_iterator<chem_opt>::set_refstate(ulong newref)
{
   if(refstate!=newref)
   {
      refstate=newref;
      occupy(refstate);
      compute_bases();
   }
   return refstate;
}

template<>
long general_base_iterator<chem_opt>::operator=(long i)
{
   if(i!=state && i<bases.size()) {
      state=i;
      while(state<bases.size() &&
            bases[state]==0)
         state++;
   }
   return state;
}

template<>
long general_base_iterator<chem_opt>::get_state() const
{
   return state;
}

template<>
long general_base_iterator<chem_opt>::operator()() const
{
   if(state<bases.size())
      return bases[state];
   else
      throw domain_error("long general_base_iterator<chem_opt>::operator()(): invalid state to call method.");
}

template<>
long general_base_iterator<chem_opt>::non_empty_size() const
{
   return number_of_bases;
}
