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
//! \file zmat.cc \brief Define zmat and zmat_entry.

#include <sstream>
#include <zmat.hh>
#include <BCR_CPP_LA/linear_algebra.h>

using namespace linear_algebra;

const Bool zmat_entry::default_opt_val[3]={Bool(false),Bool(false), Bool(false)};


void zmat_entry::update_variables(const zmat_entry& b)
{
   for(long i=0;i<3;i++)
      variable[i]=b.variable[i];
}


//! Copy assignment.
zmat_entry& zmat_entry::operator=(const zmat_entry& a) {
   Name=a.Name;
   variable=a.variable;
   connect=a.connect;
   opt_val=a.opt_val;
   increment=a.increment;
   return (*this);
}

//! Comparison operator.
bool zmat_entry::operator==(const zmat_entry& a) const {
   if(Name==a.Name &&
         variable==a.variable &&
         increment==a.increment &&
         connect==a.connect)
      return true;
   else return false;
}

//! Default constructor.
zmat_entry::zmat_entry(): Name(""), variable(3), increment(3),
      connect(3), opt_val(3,default_opt_val),
      Name_r(Name),
      variable_r(variable), increment_r(increment),
      connect_r(connect), opt_val_r(opt_val)
{
   for(int i=0;i<3;i++)
      connect[i]=-3+i;
}

//! Copy constructor.
zmat_entry::zmat_entry(const zmat_entry& a): 
	      Name(a.Name), variable(a.variable), increment(a.increment),
	      connect(a.connect), opt_val(a.opt_val),Name_r(Name),
	      variable_r(variable), increment_r(increment),
	      connect_r(connect), opt_val_r(opt_val) {};

//! Named default construction.
zmat_entry::zmat_entry(const string& N): Name(N), variable(3), increment(3),
      connect(3), opt_val(3,default_opt_val),Name_r(Name),
      variable_r(variable), increment_r(increment),
      connect_r(connect), opt_val_r(opt_val)
{
   for(int i=0;i<3;i++)
      connect[i]=-3+i;
}

//! Construction with full initialization.
zmat_entry::zmat_entry(const string& N, const double* v, const int* c):
         Name(N), variable(3), increment(3),
         connect(3), opt_val(3,default_opt_val),Name_r(Name),
         variable_r(variable), increment_r(increment),
         connect_r(connect), opt_val_r(opt_val)
{
   for(int i=0;i<3;i++) {
      variable[i]=v[i];
      connect[i]=c[i];
   }
}

//! Construction with variable initialization.
zmat_entry::zmat_entry(const string& N, const double* v):
         Name(N), variable(3), increment(3),
         connect(3), opt_val(3,default_opt_val),Name_r(Name),
         variable_r(variable), increment_r(increment),
         connect_r(connect), opt_val_r(opt_val)
{
   Name=N;
   connect.resize(3);
   variable.resize(3);
   for(int i=0;i<3;i++) {
      variable[i]=v[i];
      connect[i]=-3+i;
   }
}

//! Construction with variable initialization.
zmat_entry::zmat_entry(const string& N, const refvector<double>& v):
         Name(N), variable(v), increment(3),
         connect(3), opt_val(3,default_opt_val),Name_r(Name),
         variable_r(variable), increment_r(increment),
         connect_r(connect), opt_val_r(opt_val)
{
   if(v.dim() != 3) 
      throw domain_error("initialization of zmat_entry is wrong.");
}

zmat_entry::zmat_entry(const string& N,
      const refvector<double>& v,
      const refvector<long>& c):
         Name(N), variable(v), increment(3),
         connect(c), opt_val(3,default_opt_val),Name_r(Name),
         variable_r(variable), increment_r(increment),
         connect_r(connect), opt_val_r(opt_val)
{
   if(v.dim() != 3 || c.dim() !=3) 
      throw domain_error("initialization of zmat_entry is wrong. Name:"+N+" ");
}

zmat_entry::zmat_entry(stringstream& s):
         variable(3), increment(3), connect(3), opt_val(3,default_opt_val),Name_r(Name),
         variable_r(variable), increment_r(increment),
         connect_r(connect), opt_val_r(opt_val)
{
   static string serr="zmat_entry(stringstream): incorrect file format";
   char c;
   double v;
   s >> c; // This should be (
   getline(s,Name,',');
   s >> connect[0] >> c;
   if(c=='(') {
      s >> c;
      if(c!='0')
         opt_val[0]=true;
      s >> c;
      s >> c;
      if(c!=',')
         throw domain_error(serr+"No proper delimiter for optimization flag");
   }
   s >> variable[0] >> c;
   if(c!='(' && c!=',')
      throw domain_error(serr+"No '(' or ',' for first variable.");
   if(c!=',') {
      if(s.peek()!=')')
         while (c!=')' && s.good())
         {
            s >> v;
            increment[0].push_back(v);
            s >> c;
         }
      else
         s >> c;
      s >> c; // ,
   }

   if(!s.good())
      throw domain_error(serr+"No ')' or ',' for first variable before EOF.");
   s >> connect[1] >> c;
   if(c=='(') {
      s >> c;
      if(c!='0')
         opt_val[1]=true;
      s >> c;
      s >> c;
      if(c!=',')
         throw domain_error(serr+"No proper delimiter for optimization flag");
   }
   s >> variable[1] >> c;
   if(c!='(' && c!=',')
      throw domain_error(serr+"No '(' for second variable.");
   if(c!=',') {
      if(s.peek()!=')')
         while (c!=')' && s.good())
         {
            s >> v;
            increment[1].push_back(v);
            s >> c;
         }
      else
         s >> c;
      s >> c;
   }

   if(!s.good())
      throw domain_error(serr+"No ')' for second variable before EOF.");
   s >> connect[2] >> c;
   if(c=='(') {
      s >> c;
      if(c!='0')
         opt_val[2]=true;
      s >> c;
      s >> c;
      if(c!=',')
         throw domain_error(serr+"No proper delimiter for optimization flag");
   }
   s >> variable[2] >> c;
   if(c!='('&& c!=')')
      throw domain_error(serr+"No '(' or ')' for third variable.");
   if(c!=')') {
      if(s.peek()!=')')
         while (c!=')' && s.good())
         {
            s >> v;
            increment[2].push_back(v);
            s >> c;
         }
      else
         s >> c;
      s >> c;
   }

   if(!s.good())
      throw domain_error(serr+"No ')' for third variable before EOF.");
   if(!s.good())
      throw domain_error(serr+"string bad too early:" + Name+ " ");
}

stringstream& zmat_entry::output(stringstream& s) const
{

   s << "(" << Name;
   for(long j=0;j<2;j++) {
      s << ","
            << connect[j] << ","<< variable[j] << "(";
      for(long k=0;k<increment[j].size();k++)
         s << increment[j][k] << " ";
      s << ")";
   }
   s << ")" << endl;
   return s;
}

stringstream& zmat_connector::output(stringstream& s) const 
{
   long i;
   s << "((";

   s << centers_r[0] << ","
         << centers_r[1] << ","
         << centers_r[2] << ")";
   for(i=0;i<3;i++) {
      s << "(";
      s << modifiers_r[i][0] << ","
            << modifiers_r[i][1] << ","
            << modifiers_r[i][2] << ")";
   }
   for(i=0;i<3;i++) {
      s << "(";
      s << opt_val_r[i][0] << ","
            << opt_val_r[i][1] << ","
            << opt_val_r[i][2] << ")";
   }
   s << "(";
   for(i=0;i<angles.size();i++)
      s << angles[i] << " ";
   s << ")";
   s << ")";
   return s;
}

//! Default constructor.
zmat_connector::zmat_connector():
         centers(3), modifiers(3,3),
         opt_val(zmat_connector::default_opt_val()),
         angles(0),
         centers_r(centers), modifiers_r(modifiers),
         opt_val_r(opt_val), angles_r(angles)
{
   centers[0]=-3;
   centers[1]=-2;
   centers[2]=-1;
}

//! Construct from a zmat_entry.
zmat_connector::zmat_connector(const zmat_entry& a):
         centers(3), modifiers(3,3),
         opt_val(3,3),
         angles(0),
         centers_r(centers), modifiers_r(modifiers),
         opt_val_r(opt_val),
         angles_r(angles)
{
   centers.copy(a.connect_r);
   opt_val=zmat_connector::default_opt_val();
   for(int i=0;i<3;i++)
      modifiers[i][i]=a.variable_r[i];
}

//! Construct from a zmat_entry and a set of optimization flags.
zmat_connector::zmat_connector(const zmat_entry& a, const mat_full<Bool>& nopt_val):
         centers(3), modifiers(3,3),
         opt_val(nopt_val),
         angles(0),
         centers_r(centers), modifiers_r(modifiers),
         opt_val_r(opt_val),
         angles_r(angles)
{
   centers.copy(a.connect_r);
   for(int i=0;i<3;i++)
      modifiers[i][i]=a.variable_r[i];
}


//! Construct from an input string stream.
zmat_connector::zmat_connector(stringstream& s):
         centers(3), modifiers(3,3),
         opt_val(default_opt_val()),
         angles(0),
         centers_r(centers), modifiers_r(modifiers),
         opt_val_r(opt_val),
         angles_r(angles)
{
   static string serr="zmat_connector(stringstream): incorrect file format ";
   char c;
   s >> c; // This should be (
   s >> c;
   if(c==')' && s.good()) {
      *this=zmat_connector();
      return;
   }
   if(!s.good() || c != '(')
      throw domain_error(serr+"Bad opening, no '(' ");
   string interim;
   stringstream s2;

   // Get the connecting centers
   getline(s,interim,')');
   s2 << interim;
   s2 >> centers[0] >> c;
   if(!s2.good() || c != ',')
      throw domain_error(serr+interim);
   s2 >> centers[1] >> c;
   if(!s2.good() || c != ',')
      throw domain_error(serr+interim);
   s2 >> centers[2];
   if(!s2.good() && !s2.eof())
      throw domain_error(serr+interim);

   // Read in the connection modifiers.
   for(int i=0;i<3; i++) {
      s >> c;
      if(!s.good() || c != '(')
         throw domain_error(serr);
      s2.clear();
      s2.str("");	
      getline(s,interim,')');
      s2 << interim;
      s2 >> modifiers[i][0] >> c;
      if(!s2.good() || c != ',')
         throw domain_error(serr+" "+interim);
      s2 >> modifiers[i][1] >> c;
      if(!s2.good() || c != ',')
         throw domain_error(serr+" "+interim);
      s2 >> modifiers[i][2];
      if(!s2.good() && !s2.eof())
         throw domain_error(serr+" "+interim);
   }

   // Read in the connection optimization values.
   for(int i=0;i<3; i++) {
      s >> c;
      if(!s.good() || c != '(')
         throw domain_error(serr);
      s2.clear();
      s2.str("");	
      getline(s,interim,')');
      s2 << interim;
      bool b0;
      s2 >> b0 >> c;
      opt_val[i][0]=b0;
      if(!s2.good() || c != ',')
         throw domain_error(serr+" "+interim);
      s2 >> b0 >> c;
      opt_val[i][1]=b0;
      if(!s2.good() || c != ',')
         throw domain_error(serr+" "+interim);
      s2 >> b0;
      opt_val[i][2]=b0;
      if(!s2.good() && !s2.eof())
         throw domain_error(serr+" "+interim);
   }
   // Read in the angle and occupation.
   s >> c;
   if(!s.good() || c != '(')
      throw domain_error(serr+"Missing ( in angles occupation ");
   double angle;
   if(s.peek()!=')')
      while(s.good() && c!=')') {
         s >> angle;
         angles.push_back(angle);
         s >> c;
      }
   else
      s >> c;

   if(!s.good() || c != ')')
      throw domain_error(serr+"Missing ) in angles occupation");
   s>> c;
   if((!s.good() && !s.eof())|| c!=')')
      throw domain_error(serr + "Bad closing no ')' ");
}

//! Construct from another zmat_connector.
zmat_connector::zmat_connector(const zmat_connector& A):
         centers(A.centers), modifiers(A.modifiers),
         opt_val(A.opt_val),
         angles(A.angles),
         centers_r(centers), modifiers_r(modifiers),
         opt_val_r(opt_val),
         angles_r(angles)
{};

//! Assignment operator.
zmat_connector& zmat_connector::operator=(const zmat_connector& A)
{
   centers.copy(A.centers);
   modifiers.copy(A.modifiers);
   opt_val.copy(A.opt_val);
   angles.copy(A.angles);
   return *this;
}

//! Comparison operator.
bool zmat_connector::operator==(const zmat_connector& A) const
      {
   if(centers==A.centers &&
         modifiers==A.modifiers &&
         A.angles==angles &&
         opt_val==A.opt_val)
      return true;
   else
      return false;
      }

void zmat_connector::set_opt_val(long i, long j, Bool v)
{
   if(i<0 || i>2 || j<0 || j>2)
      throw domain_error("zmat_connector::set_opt_val out of range");
   opt_val.set(i,j,v);
   return;
}
//! Update a connector to fit after two Z-matrices have been combined.
zmat_connector& zmat_connector::update_connector(const zmat_connector& a, 
      const zmat_connector& e,
      const long add,
      zmat_connector& x)
{
   x=a;
   for(int j=0;j<3;j++) {
      if(x.centers[j]<-3)
         x.centers[j]+=3;
      else if(x.centers[j]<0)
      {
         x.modifiers[j]+=e.modifiers[x.centers[j]+3];
         for(int i=0;i<3;i++) {
            bool interimbool=
                  x.opt_val[j][i] || e.opt_val[x.centers[j]+3][i];
            x.opt_val[j][i]=interimbool;
         }
         if(j==2 && x.centers[j]==-1) {
            x.angles.copy(e.angles);
         }
         x.centers[j]=e.centers[x.centers[j]+3];
      }
      else
         x.centers[j]+=add;
   }
   return x;
}

const Bool* zmat_connector::default_y() {
   static const Bool default_y_[3]={Bool(false), Bool(false), Bool(false)};
   return default_y_;
}

const refvector<Bool>& zmat_connector::default_yrv()
{
   static const refvector<Bool> default_yrv_(3,default_y());
   return default_yrv_;
}

const refvector<Bool>* zmat_connector::default_yrvv() {
   static const refvector<Bool> default_yrvv_[3]={
         zmat_connector::default_yrv(),
         zmat_connector::default_yrv(),
         zmat_connector::default_yrv()};
   return default_yrvv_;
}
const refvector<refvector<Bool> >& zmat_connector::default_yrvrv() {
   static const refvector<refvector<Bool> > default_yrvrv_(3,zmat_connector::default_yrvv());
   return default_yrvrv_;
}

const mat_full<Bool>& zmat_connector::default_opt_val() {
   static const mat_full<Bool> default_opt_val_(3,3,zmat_connector::default_yrvrv());
   return default_opt_val_;
}

//! Update a connector to fit after two Z-matrices have been combined.
zmat_connector& zmat_connector::update_connector(const zmat_connector& a, 
      const long add,
      zmat_connector& x)
{
   x=a;
   for(int j=0;j<3;j++) {
      if(x.centers[j]>=0)
         x.centers[j]+=add;
   }
   return x;
}
