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
//! \file chem_opt.cc Implementation of chemical optimization class.

#include <BCR_CPP_LA/refcount.h>
#include <chemgroup.hh>
#include <zmat_opt.hh>
#include <chem_opt.hh>
#include <cmath>
#include <binary_line_search.hh>

using namespace std;
using namespace linear_algebra;

#define CLASS template <class P> chem_opt

//! Assignment operator.
chem_opt& chem_opt::operator=(const ChemGroup& a)
/*!
   This should really never be invoked on anything that has already been started.
 */
{
   ChemGroup::operator=(a);
   return *this;
}

//! Constructor of optimization object from a ChemGroup.
chem_opt::chem_opt(const ChemGroup& a):
               ChemGroup(a)
{};

//! Default constructor
chem_opt::chem_opt(): ChemGroup(), Library_data()
{};

//! Copy constructor
chem_opt::chem_opt(const chem_opt& a):
            ChemGroup(a),
            Library_data(a)
{};

void chem_opt::output() const
{
   ChemGroup::output();
}

/*!
   Computes the Z-matrix of molecule i and then performs a conformational search
   and finally computes and returns the computed property value (as well as
   constraint violations).
 */
valerg chem_opt::compute_property(const ulong i) const
{
   try {
      long j=visited.contains(i);

      if(j>-1) return value[j];

      zmat_connector dummy1,dummy2;
      dummy1.set_opt_val(0,0,false);
      dummy1.set_opt_val(0,1,false);
      dummy1.set_opt_val(0,2,false);

      zmat Z;

      occupy(i);
      build_zmat(0,
            dummy1,
            Z,
            dummy2);
      stringstream s;
      s << Name << i << "_";

      binary_line_search<noprune<zmat_opt> >	opt_object;
      opt_object.set_number_of_constraints(get_number_of_constraints());

      (zmat_opt&) opt_object=Z;

      cout << "Conformational Analysis of " << i << endl;

      opt_object.zmat_opt::set_Name(s.str());
      opt_object.set_id(opt_object.Name_r+"::binary_line_search<noprune<zmat_opt> >::Conformational Analysis");
      opt_object.zmat_opt::set_compute_property_flag(false);

      if(!opt_object.pre_opt(0))
      {
         cout << opt_object.id_r << " of " << i << " failed\n";
         // Added here that we still push this back.
         return get_badval();
      }

      ulong config=opt_object.optimize(0);
      cout << opt_object.id_r << " of " << i << " done!\n";
      cout.flush();
      opt_object.set_compute_property_flag(true);
      valerg val=opt_object.compute_property(config);

      visited.push_back(i);
      value.push_back(val);
      return val;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("chem_opt::compute_property(const ulong i) const");
   }
}


//! Compute the size of the optimization space assuming no constraints.
ulong chem_opt::get_space_size(const long Group) const
{
   try {
      long i,j;
      ulong my_space_size=1;
      for(j=0;j<Substituent_Groups_r[Group].allowed_Substituents_r.size();j++) {
         ulong interim_size=0;
         for(i=0;i<Substituent_Groups_r[Group].allowed_Substituents_r[j].size(); i++)
            interim_size+=get_space_size(Substituent_Groups_r[Group].allowed_Substituents_r[j][i]);
         my_space_size*=interim_size;
      }

      return my_space_size;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("chem_opt::get_space_size(const long Group) const");
   }
}

//! space_size as defined by the constraints.
ulong chem_opt::get_space_size() const
{
   try {

      if(space_size_computed) return space_size;

      long i,j;
      space_size_computed=true;
      space_size=1;
      for(i=0;i<Substituent_Groups_r.size();i++)
         for(j=0;j<Substituent_Groups_r[i].allowed_Substituents_r.size();j++)
            if(Substituent_Groups_r[i].allowed_Substituents_r[j].size()>0)
               space_size*=Substituent_Groups_r[i].allowed_Substituents_r[j].size();

      return space_size;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("chem_opt::get_space_size() const");
   }
}

//! Compute the number of bits to address the optimization space.
ulong chem_opt::get_bits() const
{
   try {
      get_space_size();

      if(bits_computed) return bits;

      bits_computed=true;
      bits=(ulong) ceil(log((double) space_size)/log(2.0));
      return bits;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("chem_opt::get_bits() const");
   }
}
