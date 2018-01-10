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
/*!@{*/
/*! \file Library_data.decl \brief Abstract Library class for
  discrete spaces.
 */


#ifndef _LIBRARY_DATA_DECL
#define _LIBRARY_DATA_DECL

#include <BCR_CPP_LA/refcount.decl>
#include <typedefs.hh>
#include <string>

using namespace linear_algebra;
using namespace std;

//! Abstract class for discrete optimizations.
class Library_data
/*!
   The library is considered an access class to a static Library. Therefore,
   Many of its members are mutable since the the static library is too large to
   be stored completely. Instead, entries are computed on the fly and memoized.
 */
{
private:
   long number_of_constraints;
   valerg badval;
   valerg set_badval() const
   {
      valerg val;
      val.energy=INFINITY;
      val.property=-INFINITY;
      val.property_computed=false;
      val.energy_computed=false;
      val.penalty=refvector<double>(number_of_constraints);
      for(long i=0;i<number_of_constraints;i++)
         val.penalty[i]=INFINITY;
      return val;
   }
protected:

   //! Visited numbers
   mutable refvector<ulong> visited;

   //! Values of visited numbers
   mutable refvector<valerg> value;

   mutable ulong space_size;
   mutable bool space_size_computed;
   mutable ulong bits;
   mutable bool bits_computed;

   mutable string Name;

   //! Class must have a way to compute values.
   virtual valerg compute_property(ulong i) const=0;

   mutable bool compute_property_flag;

public:
   //! Access to Name (Read-Only)
   const string &Name_r;

   //! Visited numbers  (Read-Only)
   const refvector<ulong> &visited_r;

   //! Values of visited numbers (Read-Only)
   const refvector<valerg> &value_r;

   //! Default Constructor
   Library_data():
      number_of_constraints(0),
      badval(set_badval()),
      visited(),
      value(),
      space_size(0),
      space_size_computed(false),
      bits(0),
      bits_computed(false),
      Name(""),
      compute_property_flag(false),
      Name_r(Name),
      visited_r(visited),
      value_r(value)
   {};

   //! Copy constructor
   Library_data(const Library_data& a):
      number_of_constraints(a.number_of_constraints),
      badval(set_badval()),
      visited(a.visited),
      value(a.value),
      space_size(0),
      space_size_computed(false),
      bits(0),
      bits_computed(false),
      Name(a.Name_r),
      compute_property_flag(a.compute_property_flag),
      Name_r(Name),
      visited_r(visited),
      value_r(value)
   {};

   //! Destructor
   virtual ~Library_data() {};

   //! Assignment operator.
   Library_data& operator=(const Library_data& d);

   //! Return the number of constraints.
   long get_number_of_constraints() const { return number_of_constraints;}

   //! Sets the number of constraints. The return value indicates success.
   bool set_number_of_constraints(long n) { number_of_constraints = n; return true;}


   //! Compute the size of the optimization space.
   virtual ulong get_space_size() const=0;

   //! Compute the number of bits to address the optimization space.
   virtual ulong get_bits() const=0;

   //! Retrieve a value for a number.
   valerg get_value(ulong i) const
   {
      return compute_property(i);
   }
   const valerg& get_badval() const
   {
      return badval;
   }

   //! Check whether a value is bad
   bool is_badval(const valerg& val) const
   {
      bool Bad=false;
      Bad = Bad || val.energy==INFINITY;
      Bad = Bad || val.property==-INFINITY;
      Bad = Bad || !val.property_computed;
      Bad = Bad || !val.energy_computed;
      for(long i=0;i<val.penalty.size();i++)
         Bad = Bad || val.penalty[i]==INFINITY;
      return Bad;
   }
   //! Sets the name for the purpose of writing files etc.
   void set_Name(const string& A) const;

   //! Set whether only energies(=false) or also properties are computed (=true)
   void set_compute_property_flag(bool B) const { compute_property_flag=B; }
   /*! Changes the behavior of compute_property.
	 - true: compute property and energy
	 - false: compute energy only.
    */
};
#endif
/*! @}*/
