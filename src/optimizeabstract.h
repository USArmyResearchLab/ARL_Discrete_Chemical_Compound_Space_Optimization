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
//! \file optimizeabstract.h \brief implement abstract class optimize_abstract

#ifndef _OPTIMIZEABSTRACT_H_
#define _OPTIMIZEABSTRACT_H_

//! Abstract class that implements the has_optimize concept. Only good for generic referencing of optimizations.
class optimize_abstract
{

private:

   //! Identifier string;
   mutable string id;
public:
   const string& id_r;

   virtual ~optimize_abstract() {};
   //! Default optimize.
   ulong optimize() const { return optimize(0); }
   virtual ulong optimize(const ulong N) const=0;
   optimize_abstract() : id(), id_r(id) {};
   optimize_abstract(const optimize_abstract& a) :
      id(a.id_r), id_r(id)
   {};
   void set_id(const string& ID) const
   {
      id=ID;
   }
   virtual void set_compute_property_flag(bool B) const=0;

protected:
};

#endif // _OPTIMIZEABSTRACT_H_
