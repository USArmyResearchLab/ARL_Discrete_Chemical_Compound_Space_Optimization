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
/*! \file chemgroup.cc \brief Implementation file of class ChemGroup */

#include <BCR_CPP_LA/refcount.h>
#include <chemgroup.hh>
#include <zmat.hh>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>

using namespace std;
using namespace linear_algebra;

//! This allows an external application to provide a filter of the library
extern bool filter(const ChemGroup& A);

//! Empty construction.
ChemGroup::ChemGroup():
              Substituent_Groups(),
              Substituent_Groups_r(Substituent_Groups)
{};

//! Construct a ChemGroup from a stringstream. \see ChemGroup::ChemGroup(istream& in)
ChemGroup::ChemGroup(stringstream& s):
              Substituent_Groups(),
              Substituent_Groups_r(Substituent_Groups)
{
   string serr="ChemGroup(stringstream& s):";
   string callederr=serr+"called";
   char c;
   string interim;
   try {
      interim="";
      try {
         while(s.peek()!=')' && s.good())
         {
            ChemIdent x(s);
            Substituent_Groups.push_back(x);
         }
      } catch(domain_error e) { cerr << e.what(); throw domain_error(callederr);}
      s >> c;// )
      if(c!=')')
         throw domain_error(serr+"no closing bracket on ChemGroup");
      if(!error_free())
         throw domain_error(serr+"Subgroups refer to illegal values.");
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error(serr);
   }
}

//! Check whether there are any errors in the definition.
bool ChemGroup::error_free() const
{
   long i,j,k;
   try {
      for(i=0; i<Substituent_Groups_r.size();i++)
         for(j=0;j<Substituent_Groups_r[i].allowed_Substituents_r.size();j++)
            for(k=0;k<Substituent_Groups_r[i].allowed_Substituents_r[j].size();k++)
               if(Substituent_Groups_r[i].allowed_Substituents_r[j][k]-Substituent_Groups_r.size()>=0)
                  return false;
      return true;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::error_free");
   }
}

//! Initialize from an input stream, e.g., a file.
/*!
  The input file format is as follows:\n
  ChemGroup( ChemIdent1 ChemIdent2 ... )\n
  \see ChemIdent::ChemIdent(istream& in) for ChemIdent input format and examples.
 */
ChemGroup::ChemGroup(istream& in):
              Substituent_Groups(),
              Substituent_Groups_r(Substituent_Groups)
{
   try {
      string interimstring;
      stringstream s;
      stringstream o;
      while(!in.eof())
      {
         in >> skipws >> interimstring;
         s << interimstring;
         // This is to allow a single comment
         if (s.peek() == '#') {
            char c = s.get();
            c = s.get();
            while(c!='\n'&& c!='#' && s.good())
               c = s.get();
         }
         o << s.str();
      }
      cout << o.str() << endl;
      ChemGroup a(o);
      *this=a;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::ChemGroup(istream& in)");
   }
}

//! Copy constructor.
ChemGroup::ChemGroup(const ChemGroup& a):
              Substituent_Groups(a.Substituent_Groups),
              Substituent_Groups_r(Substituent_Groups)
{};

/*! Assignment operator. */
ChemGroup& ChemGroup::operator=(const ChemGroup& a)
{
   try {
      Substituent_Groups.copy(a.Substituent_Groups);
      return *this;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::operator=");
   }
}

/*! Comparison operator. */
bool ChemGroup::operator==(const ChemGroup& a) const
{
   try {
      if(
            Substituent_Groups==a.Substituent_Groups
      )
         return true;
      else return false;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::operator==");
   }
            }

//! Add a substitution group.
long ChemGroup::add_substituent(const ChemIdent& a)
/*! 
This appends another ChemIdent object to the list of possible substution groups.
 */
{
   try {
      long l=Substituent_Groups.size();
      Substituent_Groups.push_back(a);
      return l;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::add_substituent(const ChemIdent& a)");
   }
}

//! Add a substituent.
void ChemGroup::add_substituent(long i, long j, long k)
/*!
   A new substituent is added at site j of group i. k references the global group
   that is taken from ChemGroup and linked to the substition site.
   Hence, j has to be a valid site on Substituent_Groups[i] and k as well as i
   are an index of Substituent_Groups.
 */
{
   try {
      ulong l=Substituent_Groups.size();
      if(i-l>0)
         throw domain_error("ChemGroup::add_substituent(long i, long j, long k): Group i does not exist");
      if(k-l>0)
         throw domain_error("ChemGroup::add_substituent(long i, long j, long k): Group k does not exist");
      if(Substituent_Groups[i].allowed_Substituents_r.size()-j<=0)
         throw domain_error("ChemGroup::add_substituent(long i, long j, long k): Connector j does not exist");
      if(!Substituent_Groups[i].allowed_Substituents_r[j].contains(k) && i!=k)
         Substituent_Groups[i].add_substituent(j,k);
      return;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::add_substituent");
   }
}

//! Add a vector of substituents.
refvector<long> ChemGroup::add_substituents(const refvector<ChemIdent>& a)
/*!
   This appends a list of substituents as add_substituent() does for a single group.
 */
{
   try {
      const long l=Substituent_Groups.size();
      Substituent_Groups.resize(l+a.size());
      refvector<long> v(a.size());
      for(long j=0;j<a.size();j++)
      {
         Substituent_Groups[l+j]=a[j];
         v[j]=l+j;
      }
      return v;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::add_substituents(const refvector<ChemIdent>& a)");
   }
}

//! Add a list of substituents to an identifier.
void ChemGroup::add_substituents(long i, long m, const refvector<long>& j)
/*!
   This adds a number of substituents to site m on group i. 
   Similar to add_substituent().
 */
{
   try {
      long l=Substituent_Groups.size();
      long k;
      if(i-l>0)
         throw domain_error("ChemGroup::add_substituent(long i, long m, refvector<long>& j): Group i does not exist");

      if(Substituent_Groups[i].allowed_Substituents_r.size()-m<=0)
         throw domain_error("ChemGroup::add_substituent(long i, long m, refvector<long>& j): Connector m does not exist");
      for(k=0;k<j.size();k++) {
         if(j[k]-l>0)
            throw domain_error("ChemGroup::add_substituent(long i, long j): Gourp j[k] does not exist");
         if(!Substituent_Groups[i].allowed_Substituents_r[m].contains(j[k]) && j[k]!=i)
            Substituent_Groups[i].add_substituent(m,j[k]);
      }
      return;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::add_substituents(long i, long m, const refvector<long>& j)");
   }
}

//! outputs a number of spaces
static void spaces(long n)
{
   long i;
   for(i=0;i<n;i++)
      cout << " ";
}

//! Display the information to cout.
void ChemGroup::output() const
{
   long i;
   try {
      cout << "(Subgroups(";
      for(i=0;i<Substituent_Groups.size();i++)
         Substituent_Groups[i].output();
      cout << "))";
      return;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::output()");
   }
}

//! Set occupations according to number.
void ChemGroup::occupy(ulong number) const
/*!
   Each ChemIdent has default occupations for its sites.
   This sets these to reflect the molecule referenced by number.
 */
{
   long i,j,m, groupsize;

   try {
      for(i=0;i< Substituent_Groups_r.size();i++) {
         const ChemIdent& Group((Substituent_Groups[i]));
         for(j=0;j< Group.allowed_Substituents_r.size();j++) {
            groupsize=Group.allowed_Substituents_r[j].size();
            m=number % groupsize;
            number=(number -m)/ groupsize;
            Group.occupy(j,m);
         }
      }
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::occupy");
   }
   return;
}

/*!
  This routine generates the local Z-matrix and combines it with the 
  global Z-matrix for molecule number. Due to connectivity concerns it may be necessary to define 
  appropriate dummy atoms in zmat::Z. A suggestion would be a dummy atom for each 
  substitution group to define the topology.
  The combinatorics are different in this version, as it frees all constraints imposed
  by the general structure.
   \param Group index of the group to be added to the Z-matrix
   \param number index of the molecule to be built.
   \param e defines the connection of Group to the Z-matrix
   \param A Z-matrix to attach Group topology. It will also be returned. 
   Note the non-const reference.
   \param y defines the return connector.
 */
zmat& ChemGroup::build_zmat(long Group,
      ulong& number,
      const zmat_connector& e,
      zmat& A,
      zmat_connector& y) const
{
   try {
      long i;
      zmat_connector x,z;
      y=x;
      y.set_opt_val(0,0,false);
      y.set_opt_val(0,1,false);
      y.set_opt_val(0,2,false);
      const long add=A.list.size()+A.offset_r;
      A.add_zmat(Substituent_Groups[Group].Z_r,e);
      ulong M,m;
      ulong groupsize;

      for(i=0;i<Substituent_Groups_r[Group].Connector_r.size();i++)
      {
         groupsize=Substituent_Groups_r[Group].allowed_Substituents_r[i].size();
         m=number % groupsize;
         number=(number -m)/ groupsize;

         zmat_connector::update_connector(Substituent_Groups_r[Group].Connector_r[i],e,add,x);
         zmat_connector::update_connector(y,x,0,z);

         M=Substituent_Groups_r[Group].occupation_r[i];
         M=m;
         build_zmat(Substituent_Groups_r[Group].allowed_Substituents_r[i][M],
               number,
               z,
               A,
               y);
      }
      zmat_connector::update_connector(Substituent_Groups_r[Group].return_connector_r,e,add,x);
      zmat_connector::update_connector(y,x,0,y);
      return A;
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::build_zmat");
   }
}

/*!
  This routine generates the local Z-matrix and combines it with the 
  global Z-matrix using the default values as provided by occupation. Due to connectivity concerns it may be necessary to define 
  appropriate dummy atoms in zmat::Z. A suggestion would be a dummy atom for each 
  substitution group to define the topology.

   \param Group index of the group to be added to the Z-matrix
   \param e defines the connection of Group to the Z-matrix
   \param A Z-matrix to attach Group topology. It will also be returned. 
   Note the non-const reference.
   \param y defines the return connector.
 */
zmat& ChemGroup::build_zmat(long Group,
      const zmat_connector& e,
      zmat& A,
      zmat_connector& y) const
{
   if(Group-Substituent_Groups_r.size()<0)
      try {
         long i;
         zmat_connector x,z;
         y=x;
         y.set_opt_val(0,0,false);
         y.set_opt_val(0,1,false);
         y.set_opt_val(0,2,false);
         const long add=A.list.size()+A.offset_r;
         A.add_zmat(Substituent_Groups[Group].Z_r,e);
         long M;

         for(i=0;i<Substituent_Groups_r[Group].Connector_r.size();i++)
         {
            zmat_connector::update_connector(Substituent_Groups_r[Group].Connector_r[i],e,add,x);
            zmat_connector::update_connector(y,x,0,z);

            M=Substituent_Groups_r[Group].occupation_r[i];
            build_zmat(Substituent_Groups_r[Group].allowed_Substituents_r[i][M],
                  z,
                  A,
                  y);
         }
         zmat_connector::update_connector(Substituent_Groups_r[Group].return_connector_r,e,add,x);
         zmat_connector::update_connector(y,x,0,y);
         return A;
      } catch(exception& e) {
         cerr << e.what() << endl;
         throw domain_error("ChemGroup::build_zmat");
      }
      else
         throw domain_error("ChemGroup::build_zmat: Group out of range");
}

//! Order substituent choices randomly.
/*!
   This is not the same as reordering done by
   reorder_general_base.
 */
void ChemGroup::randomize() {
   for( long i =0; i<Substituent_Groups.size();i++)
      Substituent_Groups[i].randomize();
}


//! Enumerate all compounds, build their Z-matrices and print together.
void ChemGroup::enumerate() const
{
   ulong space_size;
   try {
      long i,j;
      space_size=enumerate(0);
      cout << "Space size: " << space_size << endl;
      for(ulong M=0; M<space_size;M++) {
         zmat_connector dummy1,dummy2;
         zmat Z=zmat();
         ulong N=M;
         Z=build_zmat(0,N,dummy1,Z,dummy2);
         cout << "Molecule Number: " << M << "\n";
         stringstream ss;
         cout << Z.zmat_to_string(0,ss).str();
         cout << endl;
      }
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::enumerate() const");
   }

   return;

}

//! Return the number of options for a given substitution group.
ulong ChemGroup::enumerate(long Group) const
{
   ulong space_size=1;
   for(long i=0; i<Substituent_Groups_r[Group].allowed_Substituents_r.size(); i++)
   {
      ulong site_size=0;
      for(long j=0; j<Substituent_Groups_r[Group].allowed_Substituents_r[i].size();j++)
         site_size+=enumerate(Substituent_Groups_r[Group].allowed_Substituents_r[i][j]);
      space_size*=site_size;
   }
   return space_size;
}
