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
//! \file zmat.hh \brief Define zmat and zmat_entry.
#ifndef _ZMAT_HH
#define _ZMAT_HH

#include <sstream>
#include <BCR_CPP_LA/linear_algebra.h>

using namespace linear_algebra;

class zmat;

//! Wrapper for bool required for refvector<Bool>, because std::vector<bool> does not support direct access the same way.
class Bool
{
private:
   bool ok_;
public:
   Bool():ok_(true) {}
   Bool(const bool b):ok_(b) {}
   Bool(const Bool& b):ok_(b.ok_) {}
   Bool& operator=(const bool b) {ok_=b; return *this;}
   Bool& operator=(const Bool& b) {ok_=b.ok_; return *this;}
   operator bool() const { return ok_; }
   Bool operator==(const Bool a) const {return Bool(ok_==a.ok_);}
   Bool operator&&(const Bool a) const {return Bool(ok_&&a.ok_);}
   Bool operator!=(const Bool a) const {return Bool(ok_!=a.ok_);}
   Bool operator||(const Bool a) const {return Bool(ok_||a.ok_);}
   Bool operator==(const bool a) const {return (ok_==a);}
   Bool operator&&(const bool a) const {return (ok_&&a);}
   Bool operator!=(const bool a) const {return (ok_!=a);}
   Bool operator||(const bool a) const {return (ok_||a);}
   bool operator!()  const {return (!ok_);}
};

//! Class of Z-matrix entries.
class zmat_entry {
private:

   friend class zmat;

   //! Atom name
   string Name;
   //! Array of variables
   refvector<double> variable;
   //! Potential alternate values for each variable.
   refvector<refvector<double> > increment;
   //! Connectivity data
   refvector<long> connect;
   //! Optimization flags for each variable.
   refvector<Bool> opt_val;

public:

   //! Atom name
   const string &Name_r;
   //! Array of variables
   const refvector<double> &variable_r;
   //! Potential alternate values for each variable.
   const refvector<refvector<double> > &increment_r;
   //! Connectivity data
   const refvector<long> &connect_r;
   //! Optimization flags for each variable.
   const refvector<Bool> &opt_val_r;
   static const Bool default_opt_val[3];

   //! Copy assignment.
   zmat_entry& operator=(const zmat_entry& a);
   //! Comparison operator.
   bool operator==(const zmat_entry& a) const;

   //! @name Constructors
   //! @{

   //! Default constructor.
   zmat_entry();
   //! Copy constructor.
   zmat_entry(const zmat_entry& a);
   //! Named default construction.
   zmat_entry(const string& N);
   //! Construction with full initialization.
   zmat_entry(const string& N, const double* v, const int* c);
   //! Construction with variable initialization.
   zmat_entry(const string& N, const double* v);
   //! Construction with variable initialization.
   zmat_entry(const string& N, const refvector<double>& v);
   //! Construction from explicit values,
   zmat_entry(const string& N,
         const refvector<double>& v,
         const refvector<long>& c);
   //! Construction from a stringstream
   /*!
     Each entry in the string has the followin form:
     (Name, Integer(f),Length(double,...,double),Integer(f),Angle(double,...,double),Integer(f),Dihedral(double,...,double))

     - Name is a string.
     - f is a flag of 0 or 1.
     The inner parentheses are only needed if specific values are given.
    */
   zmat_entry(stringstream& s);

   //! @}

   //! Output the entry to a stringstream.
   stringstream& output(stringstream& s) const;

   //! Update the variables in the zmat_entry without touching connectivity or increments.
   void update_variables(const zmat_entry& b);


};

//! Class of connectors between z-matrices
class zmat_connector
/*!
A zmat_connector has three distinct references:
  - -6 to -4 refer to the original connector of the following group. (default)
  - -3 to -1 refer to values of the connector passed down to the current group.
  - >-1 refers to values in the local Z-matrix.
 */
{
private:
   //! Translation vector of -3, -2, -1 to alternate centers.
   refvector<long> centers;
   //! Modification matrix. Depending on the use of each center modifications may be different.
   mat_full<double> modifiers;
   //! Matrix of flags to set optimization.
   mat_full<Bool> opt_val;

   //! \f$\deg\f$/angle for conformational purposes.
   refvector<double> angles;
public:

   static const Bool* default_y();
   static const refvector<Bool>& default_yrv();
   static const refvector<Bool>* default_yrvv();
   static const refvector<refvector<Bool> >& default_yrvrv();
   static const mat_full<Bool>& default_opt_val();

   //! Translation vector of -3, -2, -1 to alternate centers. (READ ONLY)
   const refvector<long> &centers_r;
   //! Modification matrix. Depending on the use of each center modifications may be different. (READ ONLY)
   const mat_full<double> &modifiers_r;
   //! Matrix of flags to set optimization.
   const mat_full<Bool> &opt_val_r;

   //! \f$\deg\f$/angle for conformational purposes. (READ ONLY)
   const refvector<double> &angles_r;

   //! @name Constructors
   //! @{

   //! Default constructor.
   zmat_connector();
   //! Construct from a zmat_entry.
   zmat_connector(const zmat_entry& a);
   //! Construct from a zmat_entry and a set of optimization flags.
   zmat_connector(const zmat_entry& a, const mat_full<Bool>& nopt_val);
   //! Construct from an input string stream.
   /*!
     The connector information is generated from:\n
     ( C1,C2,C3)\n
     (V11,V12,V13)\n
     (V21,V22,V23)\n
     (V31,V32,V33)\n
     (O11,O12,O13)\n
     (O11,O12,O13)\n
     (O11,O12,O13)\n
     (Oa,Ob)\n
     Ci are integers referencing connectors,
     Vij are variables depending on the position Ci in a zmat_entry
     Oij is a flag 0/1 for optimization of this variable
    */
   zmat_connector(stringstream& s);

   //! Construct from another zmat_connector.
   zmat_connector(const zmat_connector& A);
   //!@}

   //! Assignment operator.
   zmat_connector& operator=(const zmat_connector& A);
   //! Comparison operator.
   bool operator==(const zmat_connector& A) const;

   //! Output the connector to a stringstream.
   stringstream& output(stringstream& s) const;

   //! Set the optimization flag .
   void set_opt_val(long i, long j, Bool v);
   //! Add an angle  to the
   void add_angle(double i);

   //! Update a connector to fit after two Z-matrices have been combined.
   static zmat_connector& update_connector(const zmat_connector& a,
         const zmat_connector& e,
         const long add,
         zmat_connector& x);

   //! Update a connector to fit after two Z-matrices have been combined.
   static zmat_connector& update_connector(const zmat_connector& a,
         const long add,
         zmat_connector& x);
};

//! Class for Z-matrices. Just a wrapper class around a vector of entries.
class zmat
/*!
  The Z-matrix may have negative connectivity entries in the first upper triangle
  of definition. This is for the purpose of combining Z-matrices. These entries are
  ignored when the zmat is output. The zmat entries are always checked for consistency
  upon construction.
 */
{
private:
   //! List of zmat_entries.
   refvector<zmat_entry> list2;
   //! Auxiliary offset for connecting matrices.
   long offset;

public:
   //! read-only access to list2
   const refvector<zmat_entry> &list;

   //! read-only access to offset.
   const long &offset_r;


   //! @name Constructors
   //! @{

   //! Default constructor.
   zmat(): list2(), offset(0), list(list2), offset_r(offset) {};

   //! Constructor which uses a non-zero offset.
   zmat(const long o):list2(), offset(o), list(list2), offset_r(offset) {};

   //! Copy constructor.
   zmat(const zmat& a):
      list2(a.list2),offset(a.offset), list(list2),
      offset_r(offset) {};

   //! Construction from a vector of zmat_entries.
   zmat(const refvector<zmat_entry>& z):
      list2(),
      offset(0),
      list(list2),
      offset_r(offset)
   {
      for(long i=0;i<z.size();i++)
         add_entry(z[i]);
   }

   //! Construct from a stringstream.
   zmat(stringstream& s):
      list2(),
      offset(0),
      list(list2),
      offset_r(offset)
   /*!
        The format goes (zmat_entry ... zmat_entry)
    */
   {
      static string serr="zmat(stringstream&): incorrect file format ";
      char c;
      s >> c;
      if(c!='(')
         throw domain_error(serr+"no opening (");
      zmat_entry x;
      while(s.peek() != ')' && !s.eof())
      {
         zmat_entry x(s);
         add_entry(x);
      }
      s >> c;
   }

   //! @}

   //! Assignment operator.
   zmat& operator=(const zmat& a)
   {
      list2=a.list2;
      offset=a.offset;
      return *this;
   }

   //! Comparison operator.
   bool operator==(const zmat& a) const
         {
      if(list2 ==a.list2 &&
            offset==a.offset)
         return true;
      else return false;
         }

   //! Set a value.
   zmat& set_val(int i, int j, double val)
   {
      if(i>=(long) list2.size())
         throw domain_error("zmat::set_val():i entry does not exist."+i);
      if(j<0 || j>2)
         throw domain_error("zmat::set_val(): j only 0 to 2 allowed.");
      list2[i].variable[j]=val;
      return *this;
   }

   //! Add to a value.
   zmat& add_val(int i, int j, double val)
   {
      if(i>=(long) list2.size())
         throw domain_error("zmat::set_val():i entry does not exist."+i);
      if(j<0 || j>2)
         throw domain_error("zmat::set_val():j only 0 to 2 allowed.");
      list2[i].variable[j]+=val;
      return *this;
   }

   //! Add an entry to the Z-matrix.
   zmat& add_entry(const zmat_entry& e)
   {
      for(int i=0;i<3;i++)
         if(e.connect[i]>=list.dim()+offset || e.connect[i]<-3)
         {
            stringstream serr;
            serr << "zmat::add_entry: Connector does not exist";
            e.output(serr);
            throw domain_error(serr.str());
         }
      if(e.connect[0]==e.connect[1] ||
            e.connect[1]==e.connect[2] ||
            e.connect[2]==e.connect[0])
      {
         stringstream serr;
         serr << "zmat::add_entry: Double connector reference in zmat_entry";
         e.output(serr);
         throw domain_error(serr.str());
      }

      list2.push_back(e);
      return (*this);
   }

   //! Output the z-matrix definition.
   void output(ostream& out) const
   {
      long i,j,k;
      out << "Z(\n";
      for(i=0;i<list.size();i++)
      {
         out << "(" << list[i].Name << ",";
         for(j=0;j<3;j++) {
            out << list[i].connect[j] << "," << list[i].variable[j] << "(";
            for(k=0;k<list[i].increment[j].size();k++) {
               out << list[i].increment[j][k];
               if(k<list[i].increment[j].size()-1)
                  out << ",";
            }
            out << ")";
            if(j<2)
               out << ",";
         }
         out << ")\n";
      }
      out << ")\n";
   }

   //! Add another value to entry i in position j.
   zmat& add_increment(long i,long j, double a)
   {
      if(i<list.size() && j< list[i].increment.size())
         list2[i].increment[j].push_back(a);
      else
         throw domain_error("zmat::add_increment: incorrect values.");
      return *this;
   }

   //! Combine two Z-matrices. B.offset must be zero!
   zmat& add_zmat(const zmat& B)
   /*!
        It is assumed that all B entries refer to the last three entries of (*this) or B itself.
    */
   {
      const long add=list.size()+offset;
      if(B.offset!=0)
         throw domain_error("zmat::add_zmat(const zmat&) has encountered incompatible offsets.");
      long i;
      zmat_entry x;
      for(i=0;i<B.list.size();i++)
      {
         x=B.list[i];
         x.connect[0]+=add;
         x.connect[1]+=add;
         x.connect[2]+=add;
         add_entry(x);
      }
      return (*this);
   }
   //! Set the optimization flag for zmat_entry i in position j
   void set_opt_val(long i, long j, Bool val)
   {
      if(i>list.size())
         throw domain_error("zmat::set_opt_val: i out of range");

      if(j>3 || j<0)
         throw domain_error("zmat::set_opt_val: j out of range");
      list2[i].opt_val[j]=val;
   }

   //! Combine two Z-matrices. B.offset must be zero!
   zmat& add_zmat(const zmat& B, const zmat_connector& e)
   /*!
        The first entry is defined to be exactly e with the corresponding Name, i.e.,
        [ B.list[0].Name -1 e.connect[0] -2 e.connect[1] -3 e.connect[2] ].
    */
   {
      if(B.offset!=0)
         throw domain_error("zmat::add_zmat(const zmat&, const zmat_entry&) has encountered incompatible offsets.");
      const long add=list.size()+offset;
      long i;
      zmat_entry x;

      for(i=0;i<B.list.size();i++)
      {
         x=B.list[i];
         for(int j=0;j<3;j++)
            if(x.connect[j]<0){
               x.variable[j]+=e.modifiers_r[x.connect[j]+3][j];
               if(j==2 && x.connect[j]==-1)
                  x.increment[j].concat(e.angles_r);
               const Bool bull=e.opt_val_r[x.connect[j]+3][j];
               const Bool bull2=((const zmat_entry&) x).opt_val[j];
               Bool result=(bull || bull2);
               result= (bool) result && !(bull && bull2);
               x.opt_val.set(j,result);
               x.connect[j]=e.centers_r[x.connect[j]+3];
            }
            else
               x.connect[j]+=add;
         add_entry(x);
      }
      return (*this);
   }

   //! Concat two Z-matrices with different offsets.
   zmat& concat_zmat(const zmat& B)
   {
      if(0!=B.offset-offset-list.size())
         throw domain_error("zmat::concat_zmat(const zmat&, const zmat_entry&): incompatible offsets.");
      for(long i=0;i<B.list.size();i++)
         add_entry(B.list[i]);
      return *this;
   }

   //! This has the effect of inserting blank space between the original offset matrix and the new matrix as well
   //! as updating the connectors.
   zmat& correct_zmat(const zmat_entry& x, const long newoff)
   {
      long i,j;
      for(i=0;i<list.size();i++)
      {
         for(j=0;j<3;j++)
            if(list2[i].connect[j]<0)
            {
               list2[i].variable[j]+=x.variable[list2[i].connect[j]+3];
               list2[i].connect[j]=x.connect[list2[i].connect[j]+3];
            }
            else if (list2[i].connect[j]-offset>=0)
               list2[i].connect[j]+=newoff-offset;
      }
      offset=newoff;
      return *this;
   }

   //! Count the number of constants.
   long count_constants() const
   {
      int nconstants=0;
      if((*this).list.dim()>1)
      {
         long i=1;
         for(long j=0;j<1;j++)
            if(!((*this).list[i].opt_val[j] || list[i].increment[j].size() >0)) nconstants++;
      }
      if((*this).list.dim()>2)
      {
         long i=2;
         for(long j=0;j<2;j++)
            if(!((*this).list[i].opt_val[j] || list[i].increment[j].size() >0)) nconstants++;
      }

      for(long i=3;i<(*this).list.size(); i++)
      {
         for(long j=0;j<3;j++)
            if(!((*this).list[i].opt_val[j] || list[i].increment[j].size() >0)) nconstants++;
      }
      return nconstants;
   }

   //! Count the number of variables.
   long count_variables() const
   {
      int nvars=0;
      if((*this).list.dim()>1)
      {
         long i=1;
         for(long j=0;j<1;j++)
            if((*this).list[i].opt_val[j]|| list[i].increment[j].size() >0) nvars++;
      }
      if((*this).list.dim()>2)
      {
         long i=2;
         for(long j=0;j<2;j++)
            if((*this).list[i].opt_val[j]|| list[i].increment[j].size() >0) nvars++;
      }

      for(long i=3;i<(*this).list.size(); i++)
      {
         for(long j=0;j<3;j++)
            if((*this).list[i].opt_val[j]|| list[i].increment[j].size() >0) nvars++;
      }
      return nvars;
   }

   //! Set the constants and variables of a Z-matrix from two refvectors.
   void set_constants_variables(const refvector<double>& consts,
         const refvector<double>& vars)
   {
      int nvars=0;
      int nconsts=0;

      if((*this).list.dim()>1)
      {
         long i=1;
         for(long j=0;j<1;j++)
            if((*this).list[i].opt_val[j]|| list[i].increment[j].size() >0)
               list2[i].variable[j]=vars[nvars++];
            else
               list2[i].variable[j]=consts[nconsts++];
      }

      if((*this).list.dim()>2)
      {
         long i=2;
         for(long j=0;j<2;j++)
            if((*this).list[i].opt_val[j]|| list[i].increment[j].size() >0)
               list2[i].variable[j]=vars[nvars++];
            else
               list2[i].variable[j]=consts[nconsts++];
      }

      for(long i=3;i<(*this).list.size(); i++)
         for(long j=0;j<3;j++) {
            if((*this).list[i].opt_val[j]|| list[i].increment[j].size() >0)
               list2[i].variable[j]=vars[nvars++];
            else
               list2[i].variable[j]=consts[nconsts++];
            if(nconsts-consts.size()>0) {
               cerr << "zmat::set_constants_variables(const refvector<double>& consts,"<< endl;
               cerr << "                              const refvector<double>& vars): nconsts discrepancy" << endl;
            }
            if(nvars-vars.size()>0){
               cerr << "zmat::set_constants_variables(const refvector<double>& consts,"<< endl;
               cerr << "                               const refvector<double>& vars): nvars discrepancy" << endl;
            }
         }
   }

   //! Output a Z-matrix to a stringstream.
   stringstream& zmat_to_string(long N, stringstream& output) const
   /*!
        \param N signifies the conformation of the combinatorial product of options.
    */
   {
      long m;
      output.str("");
      output.precision(2);
      output.setf(ios::fixed);

      int ndihedrals=0;
      int nconstants=0;
      refvector<double> dihedrals;
      refvector<double> constants;
      if((*this).list2.dim()>0)
         output << (*this).list[0].Name << endl;
      if((*this).list.dim()>1) {
         long i=1;
         output << (*this).list[i].Name << " ";
         for(long j=0;j<1;j++)
            if((*this).list[i].opt_val[j] ||
                  list[i].increment[j].size() >0) {
               output << (*this).list[i].connect[j]+1
                     << " dih" << ndihedrals++ << endl;
               dihedrals.push_back((*this).list[i].variable[j]);
               if(list[i].increment[j].size() >0)
               {
                  m=N % (list[i].increment[j].size()+1);
                  N=(N-m)/(list[i].increment[j].size()+1);
                  if(m>0)
                     dihedrals[ndihedrals-1]+=list[i].increment[j][m-1];
               }
            }
            else {
               output << (*this).list[i].connect[j]+1
                     << " c" << nconstants++  << endl;
               constants.push_back((*this).list[i].variable[j]);
            }
      }
      if((*this).list.dim()>2) {
         long i=2;
         output << (*this).list[i].Name << " ";
         for(long j=0;j<2;j++)
            if((*this).list[i].opt_val[j] ||
                  list[i].increment[j].size() >0) {
               output << (*this).list[i].connect[j]+1
                     << " dih" << ndihedrals++ << " ";
               dihedrals.push_back((*this).list[i].variable[j]);
               if(list[i].increment[j].size() >0)
               {
                  m=N % (list[i].increment[j].size()+1);
                  N=(N-m)/(list[i].increment[j].size()+1);
                  if(m>0)
                     dihedrals[ndihedrals-1]+=list[i].increment[j][m-1];
               }
            }
            else {
               output << (*this).list[i].connect[j]+1
                     << " c" << nconstants++ << " ";
               constants.push_back((*this).list[i].variable[j]);
            }
         output << endl;
      }

      for(long i=3;i<(*this).list.size(); i++) {
         output << (*this).list[i].Name << " ";
         for(long j=0;j<3;j++)
            if((*this).list[i].opt_val[j] ||
                  list[i].increment[j].size() >0) {
               output << (*this).list[i].connect[j]+1
                     << " dih" << ndihedrals++ << " ";
               dihedrals.push_back((*this).list[i].variable[j]);
               if(list[i].increment[j].size() >0)
               {
                  m=N % (list[i].increment[j].size()+1);
                  N=(N-m)/(list[i].increment[j].size()+1);
                  if(m>0)
                     dihedrals[ndihedrals-1]+=list[i].increment[j][m-1];
               }
            }
            else {
               output << (*this).list[i].connect[j]+1
                     << " c" << nconstants++  << " ";
               constants.push_back((*this).list[i].variable[j]);
            }
         output << endl;
      }

      if(ndihedrals>0)
      {
         output << endl;
         for(long i=0;i<ndihedrals;i++)
            output << "dih" << i << " " << dihedrals[i] << endl;
      }

      if(nconstants>0)
      {
         output << endl;
         for(long i=0;i<nconstants;i++)
            output << "c" << i << " " << constants[i] << endl;
      }
      return output;
   }

   //! Update the variables in the matrix without touching connectivity or increments.
   void update_variables(const zmat& B)
   {
      if(B.list.size()!=list.size())
         throw domain_error("zmat::update_variables(const zmat& B): trying to update with incompatible zmat");
      long i;
      for(i=0;i<list2.size();i++)
      {
         list2[i].update_variables(B.list[i]);
      }
   }
};
#endif
