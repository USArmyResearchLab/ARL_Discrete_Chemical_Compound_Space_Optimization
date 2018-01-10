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
/*! \file typedefs.h \brief contains typedef information */

#ifndef _TYPEDEFS_H
#define _TYPEDEFS_H

#include <cmath>
#include <limits>
#include <BCR_CPP_LA/refcount.decl>


typedef unsigned long ulong;

//! Collects property value, constraint violations, and energy.
typedef struct { double property; linear_algebra::refvector<double> penalty; double energy; bool property_computed; bool energy_computed;} valerg;

//! Pair of configurational and conformational index, respectively.
typedef struct { int configuration; int conformation;} double_index;

//! Comparison operator.
inline bool operator==(const double_index& a, const double_index& b) {
   if(a.configuration!=b.configuration || a.conformation!=b.conformation) return false;
   else return true;
}
//! Comparison operator.
inline bool operator==(const valerg& a, const valerg& b) {
   if(a.property!=b.property || a.energy!=b.energy || a.property_computed!=b.property_computed|| a.energy_computed!=b.energy_computed) return false;
   else return true;
}

//! Subtract two valerg.
inline valerg operator-(const valerg& a, const valerg& b)
{
   valerg r;
   r.property=a.property-b.property;
   r.energy=a.energy-b.energy;
   r.penalty=a.penalty-b.penalty;
   r.property_computed=a.property_computed && b.property_computed;
   r.energy_computed=a.energy_computed && b.energy_computed;
   return r;
}
//! Subtract two valerg.
inline valerg operator-=(valerg& a, const valerg& b)
      {
   a.property-=b.property;
   a.energy-=b.energy;
   a.penalty-=b.penalty;
   a.property_computed=a.property_computed && b.property_computed;
   a.energy_computed=a.energy_computed && b.energy_computed;
   return a;
      }
//! Subtract two valerg.
inline valerg operator+=(valerg& a, const valerg& b)
      {
   a.property+=b.property;
   a.energy+=b.energy;
   a.penalty+=b.penalty;
   a.property_computed=a.property_computed && b.property_computed;
   a.energy_computed=a.energy_computed && b.energy_computed;
   return a;
      }

//! Add two valerg.
inline valerg operator+(const valerg& a, const valerg& b)
{
   valerg r;
   r.property=a.property+b.property;
   r.energy=a.energy+b.energy;
   r.penalty=a.penalty+b.penalty;
   r.property_computed=a.property_computed && b.property_computed;
   r.energy_computed=a.energy_computed && b.energy_computed;
   return r;
}

//! Multiply
inline valerg operator*(const valerg& a, const double b)
{
   valerg r;
   r=a;
   r.property*=b;
   r.penalty*=b;
   r.energy*=b;
   r.energy_computed=a.energy_computed;
   r.property_computed=a.property_computed;
   return r;
}
//! Multiply
inline valerg& operator*=(valerg& r, const double b)
      {
   r.property*=b;
   r.penalty*=b;
   r.energy*=b;
   return r;
      }

using namespace linear_algebra;

inline refvector<double> atan(const refvector<double>& l)
      {
   refvector<double> r(l);
   for (long i=0;i<r.size();i++)
      r[i]=atan(r[i]);
   return r;
      }
inline refvector<double> tan(const refvector<double>& l)
      {
   refvector<double> r(l);
   for (long i=0;i<r.size();i++)
      r[i]=tan(r[i]);
   return r;
      }

inline valerg atan(const valerg& a)
{
   valerg r=a;
   r.property=atan(a.property);
   r.penalty=atan(a.penalty);
   r.energy=atan(a.energy);
   return r;
}
static const refvector<double> BADPENALTY=refvector<double>();

#endif
