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
//! \file chem_opt_class.hh chem optimization class definition.

#ifndef _chem_OPT_CLASS_HH
#define _chem_OPT_CLASS_HH

#include <BCR_CPP_LA/refcount.h>
#include <typedefs.hh>
#include <chemgroup.hh>
#include <zmat.hh>
#include <Library_data.hh>

//! Chemical optimization class.
class chem_opt: public ChemGroup, public Library_data
/*!
	The class implements a Library of molecules that can be enumerated.
 */
{
public:

   //! Default constructor.
   chem_opt();

   //! Default constructor.
   chem_opt(const chem_opt& a);

   //! Construct from a ChemGroup.
   chem_opt(const ChemGroup& a);

   //! Construction assigning maximum number of threads.
   chem_opt(ulong nmax);
   //! Assignment operator.
   chem_opt& operator=(const ChemGroup& a);

   //! Output of the chemical pattern.
   void output() const;

   //! Compute the property.
   valerg compute_property(ulong i) const;

   //! Compute the size of the optimization space.
   ulong get_space_size() const;

   //! Compute the size of the optimization space of a specific group.
   ulong get_space_size(const long Group) const;

   //! Compute the number of bits to address the optimization space.
   ulong get_bits() const;

};

#endif
