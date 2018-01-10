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
/*! \file chemident.hh \brief Header file for ChemIdent class */

#ifndef _CHEMIDENT_HH
#define _CHEMIDENT_HH

#include <BCR_CPP_LA/linear_algebra.h>
#include <zmat.hh>
#include <zmat_opt.hh>
#include <iostream>
#include <sstream>
#include <typedefs.hh>

using namespace linear_algebra;


/*!\brief
  This class describes a group of substitution sites on a Z-matrix and their respective
  substitution options.
 */
/*!
  The representation of a chemical substitution pattern is subdivided into a 
  list of substitution sites which are potentially linked (return_connector).
  In the case of linkage the groups are considered to be linked in a chain,
  i.e., only to the predecessor. The return_connector has three distinct references:
  - -6 to -4 refer to the original connector of the following group. (default)
  - -3 to -1 refer to values of the connector passed down to the current group.
  - >-1 refers to values in the local Z-matrix.

  For the purpose of preoptimizing starting structures, which in the construction
  routine may be overlapped, the zmat_connector::opt_val and zmat_entry::opt_val structures
  direct the optimization of a specific entry. By default all dihedral angles between adjoined
  Z-matrices are optimized. This may cause problems in ring-structures. To change this the 
  fix_Substituent_dihedrals() call cancels this behavior. Additionally it is possible to select
  dihedrals for conformational searching via ChemIdent::dihedrals.
  zmat_connector supplies zmat_connector::angle and zmat_connector::angle_val for conformational 
  optimizations with respect to the connectors.

 */
class ChemIdent {
private:

   /*!
    This double array holds the possible substitutions for each
    substitution site.
    */
   refvector<refvector<long> > allowed_Substituents;

   /*!
    This array holds current substituents at each site
    */
   mutable refvector<long> occupation;

   /*!
    This number holds the possible number of combinations
    for the substitution sites.
    */
   mutable long Space_Size;
   /*!
    This Z-matrix holds the connectivity data for this group.
    Respective substitution sites will be added via their own
    build_zmat() calls. For unambiguous building of Z-matrices
    it may be necessary to include dummy atoms. These are easily
    omitted upon conversion to cartesian coordinates.
    */
   zmat Z;
   /*!
    This array of entries holds the connectivity data
    on the substituents within Z and the complete framework.
    */
   refvector<zmat_connector> Connector;

   /*!
    This Z-matrix entry is returned in order to connect subsequent
    groups to this chemgroup.
    */
   zmat_connector return_connector;

public:
   /*!
    This double array holds the possible substitutions for each
    substitution site. (READ ONLY)
    */
   const refvector<refvector<long> > &allowed_Substituents_r;

   /*!
    This array holds current substituents at each site. (READ ONLY)
    */
   const refvector<long> &occupation_r;

   /*!
    This Z-matrix holds the connectivity data for this group.
    Respective substitution sites will be added via their own
    build_zmat() calls. For unambiguous building of Z-matrices
    it may be necessary to include dummy atoms. These are easily
    omitted upon conversion to cartesian coordinates.
    (READ ONLY)
    */
   const zmat &Z_r;

   /*!
    This array of entries holds the connectivity data
    on the substituents within Z and the complete framework.
    (READ ONLY)
    */
   const refvector<zmat_connector> &Connector_r;

   /*!
    This Z-matrix entry is returned in order to connect subsequent
    groups to this chemgroup.
    */
   const zmat_connector& return_connector_r;

   //! @name Constructors
   //!
   //! @{

   //! Empty constructor
   ChemIdent();

   //! Primary constructor.
   ChemIdent(const string& Name);

   //! Primary constructor from a Z-matrix.
   ChemIdent(const zmat& A);

   //! Copy constructor.
   ChemIdent(const ChemIdent& a);

   //! Construction from a stringstream.
   ChemIdent(stringstream& in);

   //! Construction from a file.
   ChemIdent(istream& in);
   //! @}

   //! Assignment operator.
   ChemIdent& operator=(const ChemIdent& a);

   //! Comparison operator.
   bool operator==(const ChemIdent& a) const;

   //! Set the Z-matrix.
   ChemIdent& set_Z(const zmat& Z);

   //! Set the return connector.
   ChemIdent& set_return_connector(const zmat_connector& A);

   //! Add another dihedral for conformational searching.
   ChemIdent& add_to_dihedrals(int zentry, int maxn, int n);

   //! Set the value of a dihedral.
   ChemIdent& set_dihedrals(int dihedrals, int n);

   //! @name Site manipulations
   //!
   //! @{

   //! Increase the number of substitution sites and add a list of substituents for that site.
   ChemIdent& add_substitution_site(const refvector<long>& a,
         const zmat_connector& e);

   //! Increase the number of substitution sites and add a substituent number to that site.
   ChemIdent& add_substitution_site(long a,
         const zmat_connector& e);

   //! Randomize the order of substitutions
   void randomize();

   //! @}

   //! @name Substituent manipulations
   //!
   //! @{

   //! Add a substituent.
   ChemIdent& add_substituents(long i,
         const refvector<long>& a);

   //! Add a substituent.
   ChemIdent& add_substituent(long i, long j);

   //! @}

   //! Update a connector to fit, e.g., after combining two Z-matrices.
   static zmat_entry& update_connector(const zmat_connector& a,
         const zmat_connector& e,
         const long A,
         zmat_connector& x);

   //! Set occupation to generate the zmat.
   void occupy(long i, long j) const;

   //! Computes the size of the library of molecules described by this instance.
   long compute_space_size() const;

   //! Fix the substituent dihedrals for connection purposes.
   void fix_Substituent_dihedrals();

   //! Display to cout what's going on.
   void output() const;
};
#endif

