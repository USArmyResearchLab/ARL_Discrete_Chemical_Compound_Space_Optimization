/* -*- Mode: C++; indent-tabs-mode: nil; c++-basic-offset: 3; tab-width: 3 -*- */
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
//! \file zmat_opt.cc Z-matrix optimization class implementation.

#include <typedefs.hh>
#include <zmat.hh>
#include <zmat_opt.hh>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

extern valerg calc_property(const zmat& A, const string& out, const string& id, zmat& returnA, long nconstraints);

//! Setup the external computations for the energy and execute them.
/*! requires \verbatim ./energy_script\endverbatim. */
static valerg calc_energy(const zmat& A, const string& out, const string& id, zmat& returnA, long nconstraints)
{
   string s;
   valerg value;
   value.property=-INFINITY;
   value.energy=INFINITY;
   value.energy_computed=false;
   value.penalty=refvector<double>(nconstraints);
   {
      s=id+".zmat";
      ofstream output_file(s.c_str());
      output_file << out << endl;
      output_file.close();
   }
   s="./energy_run "+id+"\n";

   int r=system(s.c_str());
   if(r==0)
   {
      {
         s=id+".energy";
         ifstream gfile3 (s.c_str());
         gfile3 >> value.energy;
         value.energy_computed=true;
         gfile3.close();
      }
      value.penalty=refvector<double>(nconstraints);
      {
         long i;
         s=id+".rconsts"; // optimized constants.
         long nconsts=A.count_constants();
         returnA=A;
         refvector<double> consts(nconsts);
         ifstream gfile3 (s.c_str());
         for(i=0;gfile3.good() && i<nconsts;i++)
            gfile3 >> consts[i];
         gfile3.close();
         if(i<nconsts) return value;

         s=id+".rvars"; // optimized variables.
         long nvars=A.count_variables();
         refvector<double> vars(nvars);
         gfile3.open(s.c_str());
         for(i=0;gfile3.good() && i<nvars;i++)
            gfile3 >> vars[i];
         gfile3.close();
         if(i<nvars) return value;
         returnA.set_constants_variables(consts,vars);
      }
   }
   return value;
}

//! Default constructor
zmat_opt::zmat_opt():
        Z(), Z_r(Z)
{
   visited.resize(0);
   value.resize(0);
   space_size_computed=false;
   bits_computed=false;
   Library_data::Name="";
   compute_property_flag=false;
}
//! Copy constructor
zmat_opt::zmat_opt(const zmat_opt& A):
        Z(A.Z_r),Z_r(Z)
{
   visited=A.visited;
   value=A.value;
   Library_data::Name="";
   space_size_computed=false;
   bits_computed=false;
   compute_property_flag=A.compute_property_flag;
}

//! Assignment
zmat_opt& zmat_opt::operator=(const zmat_opt& A)
{
   Z=A.Z_r;
   visited=A.visited;
   value=A.value;
   Library_data::Name="";
   space_size_computed=false;
   bits_computed=false;
   compute_property_flag=A.compute_property_flag;
   return *this;
}

zmat_opt::zmat_opt(const zmat& A) :
        Z(A),Z_r(Z)
{
   visited.resize(0);
   value.resize(0);
   Library_data::Name="";
   space_size_computed=false;
   bits_computed=false;
   compute_property_flag=false;
}

const zmat& zmat_opt::operator=(const zmat& A)
{
   Z=A;
   Library_data::Name="";
   space_size_computed=false;
   bits_computed=false;
   visited.clear();
   value.clear();
   compute_property_flag=false;
   return Z_r;
}

//! Compute the energy. (Here energy is important.)
valerg zmat_opt::compute_energy(const ulong i) const
{
   long j=visited.contains(i);
   if(j>=0)
   {
      valerg val=value[j];
      val.property_computed=false;
      val.property=val.energy*(double) -1;
      val.penalty.zero();
      return val;
   }

   static stringstream s;
   s.str("");
   s << Name << i << "_"
         << visited.size();

   stringstream out;
   valerg val;
   zmat D;

   val=calc_energy(Z_r, Z.zmat_to_string(i,out).str(),s.str(),D,get_number_of_constraints());
   val.property=val.energy;
   val.property*=(double) -1;
   val.property_computed=false;

   visited.push_back(i);
   value.push_back(val);
   return val;
}

//! Compute the energy. (Here energy is important.)
valerg zmat_opt::compute_energy(const ulong i, zmat& A) const
{
   long j=visited.contains(i);
   if(j>=0)
   {
      valerg val=value[j];
      val.property_computed=false;
      val.property=val.energy*(double) -1;
      val.penalty.zero();
      return val;
   }

   static stringstream s;
   s.str("");
   s << Name << i << "_"
         << visited.size();

   stringstream out;
   valerg val;

   val=calc_energy(Z_r, Z_r.zmat_to_string(i,out).str(),s.str(),A,get_number_of_constraints());
   val.property=val.energy;
   val.property_computed=false;
   val.property*=(double) -1;

   visited.push_back(i);
   value.push_back(val);

   return val;
}


//! Compute the property. (Here energy is important.)
valerg zmat_opt::compute_property(const ulong i) const
{   if(!compute_property_flag) {
      return compute_energy(i);
   }
   long j=visited.contains(i);
   if(j>=0 && value[j].property_computed) return value[j];
   if(j<0) {
      compute_energy(i);
      j=visited.contains(i);
   }

   static stringstream s;
   s.str("");
   s << Name << i << "_"
         << j;
   stringstream out;
   valerg val;
   zmat A;
   val=calc_property(Z_r, Z_r.zmat_to_string(i,out).str(),s.str(),A,get_number_of_constraints());

   if(j>=0) {
      value[j].property=val.property;
      value[j].property_computed=val.property_computed;
      value[j].penalty=val.penalty;
   }
   else {
      visited.push_back(i);
      value.push_back(val);
   }

   return val;
}

//! Compute the property. (Here energy is important.)
valerg zmat_opt::compute_property(const ulong i, zmat& A) const
{
   if(!compute_property_flag) return compute_energy(i,A);
   long j=visited.contains(i);
   if(j>=0 && value[j].property_computed) return value[j];
   if(j<0) {
      compute_energy(i,A);
      j=visited.contains(i);
   }

   static stringstream s;
   s.str("");
   s << Name << i << "_"
         << j;
   stringstream out;
   valerg val;

   val=calc_property(Z_r, Z_r.zmat_to_string(i,out).str(),s.str(),A,get_number_of_constraints());

   if(j>=0) {
      value[j].property=val.property;
      value[j].property_computed=val.property_computed;
   }
   else {
      visited.push_back(i);
      value.push_back(val);
   }

   return val;
}

//! Compute the size of the optimization space.
ulong zmat_opt::get_space_size() const
{
   if(space_size_computed) return space_size;

   space_size_computed=true;
   long i,j;
   ulong k=1;
   for(i=2;i<Z_r.list.size();i++)
      for(j=0;j<3;j++)
         k*=Z_r.list[i].increment_r[j].size()+1;
   space_size=k;
   return k;
}

//! Compute the number of bits to address the optimization space.
ulong zmat_opt::get_bits() const
{
   if(bits_computed) return bits;


   bits_computed=true;
   bits=(ulong) ceil(log((double) get_space_size())/log(2.0));
   return bits;
}

bool zmat_opt::pre_opt(ulong N) const
{
   ulong number=N;
   valerg current_best_val;

   current_best_val.energy=INFINITY;
   current_best_val.property=-INFINITY;

   // Create a starting structure
   zmat A;
   string oldName;
   oldName=Name;
   Name+="s";
   while(current_best_val.energy==INFINITY && number<get_space_size())
   {
      current_best_val=compute_energy(number,A);
      number++;
   }

   Name=oldName;
   if(current_best_val.energy!=INFINITY) {
      stringstream sid;
      sid << "./move_script "
            << Name << "s" << number-1 << "_" << visited.size()-1 << " "
            << Name << "0_0" << endl;
      number=0;
      system(sid.str().c_str());
      Z.update_variables(A);
      visited.clear();
      value.clear();
      visited.push_back(0);
      value.push_back(current_best_val);
   }
   else return false;
   return true;
}
