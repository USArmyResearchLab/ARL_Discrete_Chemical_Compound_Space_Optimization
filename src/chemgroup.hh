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
/*! \file chemgroup.hh \brief Header file for ChemGroup class */

#ifndef _CHEMGROUP_HH
#define _CHEMGROUP_HH

#include <BCR_CPP_LA/refcount.decl>
#include <zmat.hh>
#include <iostream>
#include <sstream>
#include <typedefs.hh>
#include <chemident.hh>

using namespace linear_algebra;

typedef unsigned long ulong;


/*!\brief
  This class describes a group of substitution sites and their respective
  substitution options as well as computation history.
 */
/*!
  It assumes the existence of a
  valerg calc_property(const zmat& A, const stringstream& id) function to compute 
  properties. 
  In essence this is a collection of ChemIdent objects.
	The first substitution group is used as the root for generating any Z-matrices
	based on the specific chosen substitutions.

	The general strategy is to split each substitution group (ChemIdent) from the
	overall structure. Hence, a ChemIdent identifies a molecular framework/scaffold
	and possible substitution sites on this. Singletons don't have any substitution sites
	and cap substitutions.
	The ChemGroup collects these and gives meaning to specific substitutions at each defined 
	site.

  See class ChemIdent for the substitution rules.

 */
class ChemGroup {
private:
   /*!
    This array holds the possible ChemIdent substitutions for each
    substitution site.
    */
   refvector<ChemIdent> Substituent_Groups;

   //! Return the number of substitutions possible for a specific group.
   ulong enumerate(long Group) const;

public:

   bool error_free() const;

   //! Build the Z-matrix using occupation i.
   zmat& build_zmat(long i,
         ulong& number,
         const zmat_connector& e,
         zmat& A,
         zmat_connector& y) const;

   //! Build the Z-matrix using occupation i.
   zmat& build_zmat(long i,
         const zmat_connector& e,
         zmat& A,
         zmat_connector& y) const;

   //! Set occupations according to number.
   void occupy(ulong number) const;

   /*!
    This double array holds the possible substitutions for each
    substitution site.
    (READ ONLY)
    */
   const refvector<ChemIdent> &Substituent_Groups_r;

   //! @name Constructors
   //!
   //! @{

   //! Empty constructor
   ChemGroup();

   //! Copy constructor.
   ChemGroup(const ChemGroup& a);

   //! Construction from a stringstream.
   ChemGroup(stringstream& in);

   //! Construction from a file.
   ChemGroup(istream& in);
   //! @}

   //! Assignment operator.
   ChemGroup& operator=(const ChemGroup& a);

   //! Comparison operator.
   bool operator==(const ChemGroup& a) const;

   //! @name Substituent manipulations
   //!
   //! @{

   //! Add a substituent.
   long add_substituent(const ChemIdent& a);

   //! Add a list of substituents.
   refvector<long> add_substituents(const refvector<ChemIdent>& a);

   //! Add a substituent to a site.
   void add_substituent(long i, long j, long k);

   //! Add a list of substituents to a site.
   void add_substituents(long i, long k, const refvector<long>& j);

   //! Randomize the order of substituents
   void randomize();
   //! @}

   //! Display to cout what's going on.
   void output() const;

   //! Enumerate and print all possible combinations
   void enumerate() const;

};
#endif
