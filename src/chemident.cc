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
/*! \file chemident.cc \brief Implementation file of class ChemIdent */

#include <BCR_CPP_LA/refcount.h>
#include <chemident.hh>
#include <zmat.hh>
#include <zmat_opt.hh>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;
using namespace linear_algebra;

static const double zeros[3]={0.0,0.0,0.0};
static const int    conn[3]={-6,-5,-4};
static const zmat_connector default_return_connector(zmat_entry("_",zeros,conn));

//! Empty construction.
ChemIdent::ChemIdent():
  allowed_Substituents(),
  occupation(), Space_Size(0), Z(),
  return_connector(default_return_connector),
  allowed_Substituents_r(allowed_Substituents), 
  occupation_r(occupation),
  Z_r(Z),
  Connector_r(Connector),
  return_connector_r(return_connector)
{};

//! Default construction. \example ChemGroup ChemIdent("H")
ChemIdent::ChemIdent(const string& Name):
  allowed_Substituents(),
  occupation(), Space_Size(0), Z(),
  return_connector(default_return_connector),
  allowed_Substituents_r(allowed_Substituents),
  occupation_r(occupation),
  Z_r(Z),
  Connector_r(Connector),
  return_connector_r(return_connector)
{
  try {
    Z.add_entry(zmat_entry(Name));
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error("ChemIdent::ChemIdent(const string& Name)");
  }
}

//! Default construction from Z-matrix.
ChemIdent::ChemIdent(const zmat& A):
  allowed_Substituents(), 
  occupation(), Space_Size(0), Z(A),
  return_connector(default_return_connector),
  allowed_Substituents_r(allowed_Substituents),
  occupation_r(occupation),
  Z_r(Z),
  Connector_r(Connector), return_connector_r(return_connector)
{}

//! Construct a ChemGroup from a stringstream.
ChemIdent::ChemIdent(stringstream& s):
  allowed_Substituents(), 
  occupation(), Space_Size(0), Z(),
  return_connector(default_return_connector),
  allowed_Substituents_r(allowed_Substituents),
  occupation_r(occupation),
  Z_r(Z),
  Connector_r(Connector), 
  return_connector_r(return_connector)
/*!
  The input file format is as follows:\n
  ( 
   Z(zmat)
   \n
  ReturnConnector (zmat_connector)\n
  Connector ( zmat_connector
   zmat_connector
  ... )\n
  allowed_Substituents( ( group number 1, group number 2 , ... ) ... )
  )\n
   Under Z a zmat object is read in. The appropriate format can be found under 
   zmat::zmat(stringstream& s).
   ReturnConnector reads in a zmat_connector using 
   zmat_connector::zmat_connector(stringstream& s).
   Connector reads in a list of zmat_connector objects, each adhering to the 
   format found under zmat_connector::zmat_connector(stringstream& s).
  \see \link carbazoles.inp carbazoles.inp\endlink and \link vanilla-rings.inp vanilla-rings.inp\endlink in the examples section.

*/
{
  string serr="ChemIdent(stringstream& s):";
  string callederr=serr+"called";
  try {
    char c;
    s >> c;
    if(c!='(')
      throw domain_error(serr);

    // Read the Z-Matrix.
    s >> c;
    if(c!='Z')
      throw domain_error(serr+"Incorrect file format: Z-Matrix improperly defined.");
    try {
      zmat A(s);
      Z=A;
    } catch(domain_error e) {
       cerr << e.what(); throw domain_error(callederr);
    }
    string interim;
    //! Read the connector.
    interim="";
    getline(s,interim,'(');
    if(interim!="ReturnConnector")
      throw domain_error(serr+"Incorrect file format: at keyword 'ReturnConnector': "+interim);
    try {
      if(s.peek()!=')')
	return_connector=zmat_connector(s);
    } catch(domain_error e) { cerr << e.what(); throw domain_error(callederr);}

    s >> c; // )

    //! Read the connector.
    interim="";
    getline(s,interim,'(');
    if(interim!="Connector")
      throw domain_error(serr+"Incorrect file format: at keyword 'Connector'");
    try {
      while(s.peek()!=')' && s.good())
	Connector.push_back(zmat_connector(s));
    } catch(domain_error e) { cerr << e.what(); throw domain_error(callederr);}

    s >> c; // )

    interim="";
    getline(s,interim,'(');
    if(interim!="allowed_groups")
      throw domain_error(serr+"Incorrect file format: at keyword 'allowed_groups'");
    while(s.peek()!=')' && s.good())
      {
	s >> c;// (
	refvector<long> S;
	while(c!=')' && s.good())
	  {
	    long x;
	    s >> x;
	    S.push_back(x);
	    s >> c; // , or )
	  }
	allowed_Substituents.push_back(S);
	//      s >> c;// )
      }
    if(allowed_Substituents.size()!=Connector.size())
      throw domain_error(serr+"Connector size and allowed_groups size do not match up.");
    s >> c >> c;
    occupation.resize(Connector.size());
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error("ChemIdent::ChemIdent(stringstream& s)");
  }
}


//!  \example carbazoles.inp 

//!  \example vanilla-rings.inp 

//! Initialize from an input stream, e.g., a file.  \see ChemIdent::ChemIdent(stringstream& in)
ChemIdent::ChemIdent(istream& in):
  allowed_Substituents(),
  occupation(), Space_Size(0),
  Z(),
  return_connector(default_return_connector),
  allowed_Substituents_r(allowed_Substituents),
  occupation_r(occupation),
  Z_r(Z),
  Connector_r(Connector),
  return_connector_r(return_connector)
{
  try {
    string s;
    stringstream o;
    while(!in.eof())
      {
	in >> skipws >> s;
	o << s;
      }
    cout << o.str() << endl;
    ChemIdent a(o);
    *this=a;
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error("ChemIdent::ChemIdent(istream& in)");
  }
}

//! Copy constructor including the parent. 
ChemIdent::ChemIdent(const ChemIdent& a):
  allowed_Substituents(a.allowed_Substituents),
  occupation(a.occupation), Space_Size(a.Space_Size),
  Z(a.Z),
  Connector(a.Connector),
  return_connector(a.return_connector),
  allowed_Substituents_r(allowed_Substituents),
  occupation_r(occupation),
  Z_r(Z),
  Connector_r(Connector),
  return_connector_r(return_connector)
{};

/*! Assignment operator. */
ChemIdent& ChemIdent::operator=(const ChemIdent& a)
{
  try {
    allowed_Substituents.copy(a.allowed_Substituents);
    Z=a.Z;
    Connector=a.Connector;
    return_connector=a.return_connector;
    occupation.copy(a.occupation);
    return *this;
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error("ChemIdent::operator=");
  }
}

/*! Comparison operator. */
bool ChemIdent::operator==(const ChemIdent& a) const
{
  try {
    if(
       Z==a.Z &&
       Connector==a.Connector &&
       return_connector==a.return_connector &&
       allowed_Substituents==a.allowed_Substituents
       )
      return true;
    else return false;
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error("ChemIdent::operator==");
  }
}

//! Set the Z-matrix.
ChemIdent& ChemIdent::set_Z(const zmat& A)
{
  try {
    Z=A;
    return *this;
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error("ChemIdent::set_Z");
  }
}

//! Set the return connector.
ChemIdent& ChemIdent::set_return_connector(const zmat_connector& A)
{
  try {
    return_connector=A;
    return *this;
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error("ChemIdent::set_return_connector");
  }
}


//! Add another dihedral for conformational searching. (use with caution)
ChemIdent& ChemIdent::add_to_dihedrals(int zentry, int maxn, int n)
/*!
   \param zentry dihedral of Z-matrix entry to be optimized.
   \param maxn   division fo 360 degrees.
   \param n      current instantiation (<maxn).
*/
{
  string serr="ChemIdent::add_to_dihedrals():";
  try {
    if(zentry-Z_r.list.size()>0)
      throw domain_error(serr+"zentry out of range");
    if(n>=maxn)
      throw domain_error(serr+"n exceeds maxn");
    long i;
    for(i=1;i<maxn;i++)
      Z.add_increment(zentry,2,(double) (((i+n)%(maxn-1))+1)*360.0/(double) maxn);

    return *this;
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error(serr);
  }
}

//! Increase the number of substitution sites and add a list of substituents for that site.
ChemIdent& ChemIdent::add_substitution_site(const refvector<long>& a,
					    const zmat_connector& e)
{
  try {
    Connector.push_back(e);
    allowed_Substituents.push_back(a);
    return *this;
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error("ChemIdent::add_substitution_site(refvector a,e)");
  }
}

//! Increase the number of substitution sites and add a substituent number for that site.
ChemIdent& ChemIdent::add_substitution_site(long a,
					    const zmat_connector& e)
{
  try {
    Connector.push_back(e);
    refvector<long> r;
    r.push_back(a);
    allowed_Substituents.push_back(r);
    return *this;
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error("ChemIdent::add_substitution_site(a,e)");
  }
}

//! Add a vector of substituents.
ChemIdent& ChemIdent::add_substituents(long i,
				       const refvector<long>& a)
{
  try {
    if(i-Connector_r.size()>0)
      throw domain_error("ChemIdent::add_substituent : i exceeds site list");
    const ulong k=allowed_Substituents[i].size();
    allowed_Substituents[i].resize(k+a.size());
    for(long j=0;j<a.size();j++)
      {
	allowed_Substituents[i][k+j]=a[j];
      }
    return *this;
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error("ChemIdent::add_substituents(i,a)");
  }
}

//! Add a substituent.
ChemIdent& ChemIdent::add_substituent(long i,
				      long j)
{
  try {
    if(i-occupation_r.size()>0)
      throw domain_error("ChemIdent::add_substituent : i exceeds site list");
    if(!allowed_Substituents_r[i].contains(j))
      allowed_Substituents[i].push_back(j);
    return *this;
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error("ChemIdent::add_substituent(i,j)");
  }
}

void ChemIdent::randomize()
{
   long randomnumber=0;
   for (long i =0; i<allowed_Substituents.size();i++)
   {
      refvector<long> neworder;
      while(allowed_Substituents[i].size()>0)
      {
         randomnumber=random() % allowed_Substituents[i].size();
         neworder.push_back(allowed_Substituents[i][randomnumber]);
         allowed_Substituents[i].erase(randomnumber);
      }
      allowed_Substituents[i]=neworder;
   }
}

static void spaces(long n)
{
  long i;
  for(i=0;i<n;i++)
    cout << " ";
}

//! Display the information to cout.
void ChemIdent::output() const
{
  try {
    long i,j;
    stringstream s;
    cout << "(Z(";
    for(i=0;i<Z.list.size();i++)
      cout << "(" << Z.list[i].Name_r << ","
	   << Z.list[i].connect_r[0] << "," << Z.list[i].variable_r[0] << ","
	   << Z.list[i].connect_r[1] << "," << Z.list[i].variable_r[1] << ","
	   << Z.list[i].connect_r[2] << "," << Z.list[i].variable_r[2] << ")";
    cout << ")ReturnConnector(";
    return_connector_r.output(s);
    cout << s.str();
    cout << ")Connector(";
    for(j=0;j<Connector_r.size();j++)
      {
	s.str("");
	Connector_r[j].output(s);
	cout << s.str();
      }
    cout << ")allowed_groups(";
    for(i=0;i<allowed_Substituents_r.size();i++)
      {
	cout << "(";
	for(j=0;j<allowed_Substituents_r[i].size()-1;j++)
	  cout << allowed_Substituents_r[i][j] << ",";
	if(allowed_Substituents_r[i].size()>0)
	  cout << allowed_Substituents_r[i][j];
	cout << ")";
      }
    cout << "))";
    return;
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error("ChemIdent::output");
  }
}

//! Set default occupation of group i to m.
void ChemIdent::occupy(long i, long m) const
{
  try {
    if(i-allowed_Substituents_r.size()>=0)
      throw domain_error("ChemIdent::occupy(long i, long m): i too large");

    if(m-allowed_Substituents_r[i].size()>=0)
      throw domain_error("ChemIdent::occupy(long i, long m): m too large");

    occupation[i]=m;
  
  } catch(exception& e) {
    cerr << e.what() << endl;
    throw domain_error("ChemIdent::occupy(i,m)");
  }
}
