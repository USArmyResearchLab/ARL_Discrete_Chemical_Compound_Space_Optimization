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
//! \file zmat_opt.hh Z-matrix optimization class definition.

#ifndef _ZMAT_OPT_HH
#define _ZMAT_OPT_HH

#include <typedefs.hh>
#include <zmat.hh>
#include <Library_data.hh>
#include <noprune.h>
/*!
 zmat_opt describes a Library using class zmat.
 It has two associated properties as expressed by compute_energy(ulong i)
 and compute_property(ulong i). Set compute_property_flag in order to compute the
 property and energy of a given conformation.
 */
class zmat_opt:
      public Library_data
{
private:
   mutable zmat Z;
public:
   //! Read-only access to Z.
   const zmat &Z_r;


   //! @name Constructors
   //!@{

   //! Default constructor
   zmat_opt();
   //! Copy constructor.
   zmat_opt(const zmat_opt& A);
   //! Initialize with a zmatrix
   zmat_opt(const zmat& A);

   //!@}

   //! Assignment operator for a zmat.
   const zmat& operator=(const zmat& A);

   //! Assignment operator for a zmat_opt.
   zmat_opt& operator=(const zmat_opt& A);

   //! Compute the property. (Here energy is important.)
   valerg compute_property(ulong i) const;

   //! Compute the property. (Here energy is important.)
   /*! This method is necessary for initialisation purposes. */
   valerg compute_property(ulong i, zmat& A) const;

   //! Compute the energy. (Here energy is important.)
   valerg compute_energy(ulong i) const;

   //! Compute the energy. (Here energy is important.)
   /*! This method is necessary for initialisation purposes. */
   valerg compute_energy(ulong i, zmat& A) const;

   //! Compute the size of the optimization space.
   ulong get_space_size() const;

   //! Compute the number of bits to address the optimization space.
   ulong get_bits() const;

   //! Find a converged starting geometry.
   bool pre_opt(ulong N) const;
};

#endif
