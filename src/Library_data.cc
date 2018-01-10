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
/*! @{*/

/*! \file Library_data.cc \brief Implementations of the Library class for 
  discrete spaces.
 */


#ifndef _LIBRARY_DATA_CC
#define _LIBRARY_DATA_CC

#include <BCR_CPP_LA/refcount.h>
#include <Library_data.hh>

using namespace linear_algebra;

//! Assignment operator.
Library_data& Library_data::operator=(const Library_data& d)
{
   try {
      space_size_computed=d.space_size_computed;
      space_size=d.space_size;
      bits_computed=d.bits_computed;
      bits=d.bits;
      visited=d.visited;
      value=d.value;
      Name=d.Name;
      return *this;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("Library_data::operator=");
   }
}

//! Set the name
void Library_data::set_Name(const string& s) const
{
   Name=s;
}
#endif
/*! @}*/
