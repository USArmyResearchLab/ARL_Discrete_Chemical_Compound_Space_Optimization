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

/*! \addtogroup DOChemS Discrete Optimization of Chemical Space */

//! \file DiscreteCCSOpt.cc Main routine and setup optimization.
/*!
   This file also contains the definition of the actual computations to be fulfilled.
 */
/*!
\mainpage Introduction
\section usage_sec Usage
These modules include a general substitution class based on recursion (ChemIdent).
In the main driver, the function calc_property() executes <EM> ./property_script</EM>.
The secondary driver executes an energy computation via <EM>./energy_script</EM>.

To test the program, copy all files from the test directory to your working directory and invoke
\verbatim /path_to_executable/DiscreteCCSOpt -f ftc-input.inp \endverbatim
The result should be a simple optimization of the hyperpolarizibility.

calc_property() will generate a <EM>.zmat</EM> file which will hold the Z-matrix of the current molecule. 
The prefix will be the current number of visited configurations followed by the current conformation and overall structures visited.
See calc_property() for details.
<EM>energy_script</EM> must return <EM>.rconsts, .rvars, .energy</EM> files on a successful run.
- <EM>.rconsts</EM> must contain the optimized constants of the Z-matrix from the geometry optimization.
- <EM>.rvars</EM>  must contain the optimized variables of the Z-matrix from the geometry optimization.
- <EM>.result</EM> holds the result of the objective function.
- <EM>.energy</EM> holds the energy of the final geometry.

<EM>property_script</EM> must return <EM>.result</EM> and potentially <EM>.penalty</EM> files on a successful run.
- <EM>.result</EM> holds the result of the objective function.
- <EM>.penalty</EM> contains a penalty value(s).

To access the program flow go to main().
For information on the representation of molecules see class ChemIdent.
For information on the input file syntax see ChemIdent::ChemIdent(stringstream& s).

\subsection cmd_opt Command Line Options
\param -f <filename> use file with filename to read the input. If not given, the filename will be taken
from standard input interactively.
\param --start_compound <number>
\param --sc <number> Use molecule associated with number as starting compound. If not given, the starting number
will be taken from standard input interactively.
\param --checkinput Checks whether the input files follow correct syntax and exits.
\param -p Switches pruning or heuristics on. See noprune, simple_prune, reorder_general_base.
\param --sub-method <method>
\param -sm <method>
\newline details the sub-method to be used.
Currently this includes <EM>BLS (binary_line_search), GBLS (gen_base_LS), GBGLS (gen_base_grad_LS), SD or sd or  steepest_descent (binary_steepest_descent)</EM>
\param -m <method> \newline
Supplies the method to be used. Options include all listed under sub-methods and <EM>
GBEN (gen_base_entropic), GDMC (binary_gdmc, gen_base_gdmc)</EM>. The default is BLS.
GDMC only accepts GBGLS, and BLS, whereas GBEN also accepts GBLS in addition.
\param --pre For GBGLS, indicate that a full sweep of each search direction is to be performed prior to optimization.
This is only useful if -p is also set.
\param --atmax places substituents which have not been encountered previously in a GBGLS run with pruning at the top of the
search circle.
\param --atcurrent treats substituents which have not been encountered previously in a GBGLS run with pruning equivalent to
the active substituent at the time of reordering.
\param --minimax enforces adjustment of Lagrange multipliers \f$\lambda\f$ associated with constraints to solve the minimax problem
over the set of molecules which have already been computed.
\param -T <temperature> provides a fictitious temperature for GDMC calculations.
\param -TS <steps> provides the number of steps to tighten selection criteria for GDMC or number of runs for GBEN.
\param -MS <steps> Maximum number of evaluations used in GBEN and GDMC.
\param --random-order Order substituents randomly for each site.
\param --gben-reorder
\param --gbenr requires that GBEN with GBGLS use the order established in the last run to assess the diversity metric.
\param --enumerate Prints every Z-matrix and its corresponding library number

\section install_sec Installation
- Untar and unzip the tarball. The source will be extracted into <CODE>DiscreteCCSOpt/</CODE>.
- Change direcory to DiscreteCCSOpt/Build and edit Doxyfile to change the documentation path to your liking.
- make 
- For Documentation execute: \verbatim make docs\endverbatim
- For the debugable version execute: \verbatim make DiscreteCCSOpt.d \endverbatim

MakeVars and Makefile contain the compilation information which you may want to change.
The executable will be named <EM>DiscreteCCSOpt</EM>.
 */

#include <BCR_CPP_LA/refcount.h>
#include <chemgroup.hh>
#include <chem_opt.hh>
#include <Library_data.hh>
#include <zmat.hh>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <noprune.h>
#include <simpleprune.h>
#include <simpleprune.cc>
#include <optimizeabstract.h>
#include <binary_line_search.hh>
#include <binarysteepestdescent.hh>
#include <binarygdmc.hh>
#include <genbase-l-s.hh>
#include <genbase-grad-ls.hh>
#include <reordergeneralbase.hh>
#include <generalbaseiterator.hh>
#include <genbasegdmc.hh>
#include <gen_base_entropic.hh>
#include <binary_entropic.hh>
#include "parse.hh"

using namespace std;
using namespace linear_algebra;

//! Setup the external computations and execute them.
valerg calc_property(const zmat& A, const string& out, const string& id, zmat& returnA, long nconstraints)
/*!
 * Ensure that even failed jobs return the exact number of penalties. Otherwise, access errors will occur.
 */
{
   string s;
   valerg value;
   value.property=-INFINITY;
   value.energy=INFINITY;
   value.penalty=refvector<double>(nconstraints);
   for(long i=0;i<nconstraints;i++) value.penalty[i]=INFINITY;
   value.property_computed=false;
   value.energy_computed=false;

   s="./property_script "+id+"\n";

   int r=system(s.c_str());
   if(r==0)
   {
      value.property_computed=true;
      value.energy_computed=true;
      {
         s=id+".result";
         ifstream gfile3 (s.c_str());
         gfile3 >> value.property ;
         gfile3.close();
      }
      {
         s=id+".energy";
         ifstream gfile3 (s.c_str());
         gfile3 >> value.energy;
         gfile3.close();
      }
      {
         s=id+".penalty";
         ifstream gfile3 (s.c_str());
         double p=0.0;
         long count=0;
         while(gfile3.good() && nconstraints>count) {
            gfile3 >> p;
            value.penalty[count]=p;
            count++;
         }
         gfile3.close();
      }
      {
         long i;
         s=id+".rconsts";
         long nconsts=A.count_constants();
         returnA=A;
         refvector<double> consts(nconsts);
         ifstream gfile3 (s.c_str());
         for(i=0;gfile3.good() && i<nconsts;i++)
            gfile3 >> consts[i];
         gfile3.close();
         if(i<nconsts)
            return value;

         s=id+".rvars";
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

template <class X>
void print(const X& D)
{

   cout << "Space size: " << D.get_space_size() << endl;
   cout << "# of bits:  " << D.get_bits() << endl;
}

template<class X>
ulong run(const X& D, ulong value=0, bool value_passed=false)
{
   if(!value_passed) {
      cout << "Enter a starting occupation (Enter for default):";
      cout.flush();
      cin.clear();
      cin >> value;
   }

   value=D.optimize(value);

   cout << "The optimized value is: " << D.get_value(value).property << " for configuration " << value << endl;
   return value;
}


//! Main execution routine. Handles command-line input and executes optimization.
int main(int argc, char *argv[])
{
   stringstream s;
   try {
      ifstream in;
      string filename;
      string method;
      string submethod;
      double T=0.0;
      double value=0.0;
      long max_steps=1;
      long tight_steps=1;
      bool pruned=false;
      bool precondition_flag=false;
      bool at_max_flag=false;
      bool at_current_flag=false;
      bool minimax_flag=false;
      bool checkinput=false;
      bool randomize = false;
      bool enumerateflag=false;
      bool value_passed=false;
      bool gbenreorderflag=false;

      cout << "Invocation:" << endl;
      for(int i=0;i<argc;i++)
         cout << argv[i] << endl;

      for (int i=1;i<argc;i++)
      {
         string command=argv[i];
         if(command=="-f" && i+1<argc)
            filename=argv[i+1];
         else if(command=="-m") {
            if(argc>i+1)
               method=argv[++i];
         }
         else if(command=="-sm" || command=="--sub-method" ) {
            if(argc>i+1)
               submethod=argv[++i];
         }
         else if(command=="-T")  {
            if(argc>i+1) {
               stringstream s;
               s << argv[++i];
               s >> T;
               cout << "Used temperature: " << T << endl;
            }
         }
         else if(command=="-p")  {
            pruned=true;
         }
         else if(command=="-MS")  {
            if(argc>i+1) {
               stringstream s;
               s << argv[++i];
               s >> max_steps;
               cout << "Maximum number of steps: " << max_steps << endl;
            }
         }
         else if(command=="-TS")  {
            if(argc>i+1) {
               stringstream s;
               s << argv[++i];
               s >> tight_steps;
               cout << "Number of tightenings of constraints: " << tight_steps << endl;
            }
         }
         else if(command=="--pre"){
            precondition_flag=true;
         }
         else if(command=="--atmax"){
            at_max_flag=true;
         }
         else if(command=="--atcurrent"){
            at_current_flag=true;
         }
         else if(command=="--minimax"){
            minimax_flag=true;
         }
         else if(command=="--checkinput") {
            checkinput=true;
         }
         else if(command=="--random-order") {
            randomize=true;
         }
         else if(command=="--gben-reorder" || command=="--gbenr") {
            gbenreorderflag=true;
         }
         else if(command=="--enumerate") {
            enumerateflag=true;
         }
         else if(command=="--start_compound" || command =="--sc") {
            value_passed=true;
            if(argc>i+1) {
               stringstream s;
               s << argv[++i];
               s >> value;
               cout << "Starting compound: " << value << endl;
            }
         }
      }

      if(filename!="")
         in.open(filename.c_str());
      else {
         string s2;
         cout << "Enter an input file name: ";
         cout.flush();
         getline(cin,s2);
         in.open(s2.c_str());
      }
      if(!in.good())
      {
         cerr << "Bad filename. exiting.\n";
         exit(1);
      }
      ChemGroup Complex;
      long nconstraints=0;
      {
         refvector<section> p=parse(in);
         long l = 0;
         while(l<p.size() && p[l].get_Name()!="ChemGroup")
            l++;
         if(l>=p.size())
            throw domain_error("No ChemGroup section");
         stringstream ss;
         ss.str("");
         ss<<p[l].get_buffer();
         ChemGroup Pattern(ss);
         for(l=0;l<p.size() && p[l].get_Name()!="nconstraints";l++);
         if(l<p.size()) {
            stringstream s;
            s.str("");
            s << p[l].get_buffer();
            cout << s.str() << endl;
            s >> nconstraints;
         }
         Complex=Pattern;
      }

      Complex.output();
      if(enumerateflag)
         Complex.enumerate();
      if(checkinput) {
         cout << "Input is fine\n";
         return 0;
      }

      if(randomize) Complex.randomize();

      if(method!="") {
         if(method=="SD" || 
               method=="sd" ||
               method=="steepest_descent")
         {
            cout << "Doing steepest_descent optimization" << endl;
            cerr << "Currently not debugged. Results may be unreliable." << endl;
            if(pruned)
            {
               typedef binary_steepest_descent<simple_prune<has_gradients_data<chem_opt> > > C_t;
               chem_opt t1;
               t1.compute_property(0);
               noprune<chem_opt> t2;
               t2.compute_property(0);
               has_gradients_data<noprune<chem_opt> > t3;
               t3.compute_property(0);
               binary_steepest_descent<simple_prune<has_gradients_data<chem_opt> > > t4;
               C_t C;
               (chem_opt&) C=Complex;
               C.set_number_of_constraints(nconstraints);
               cout << "top=binary_steepest_descent<has_gradients_data<simple_prune<chem_opt> > >\n";
               C.set_id("top");
               print<C_t>(C);
               run<C_t>(C,value,value_passed);
            }
            else
            {
               typedef binary_steepest_descent<has_gradients_data<noprune<chem_opt> > > C_t;
               C_t C;
               (chem_opt) C=Complex;
               C.set_number_of_constraints(nconstraints);
               cout << "top=binary_steepest_descent<has_gradients_data<noprune<chem_opt> > >\n";
               C.set_id("top");
               print<C_t>(C);
               run<C_t>(C,value,value_passed);
            }
         }
         else if(method=="GDMC" ||
               method=="gdmc") {
            cout << "Doing GDMC optimization" << endl;
            if(max_steps<2) {
               cout << " Please enter maximum number of computations: ";
               cin >> max_steps;
            }
            if(tight_steps<2) {
               cout << " Please enter maximum number of constraint tightening: ";
               cin >> tight_steps;
            }
            if(T<=0.0) {
               cout << " Please enter the temperature: ";
               cin >> T;
            }
            srandom(0);
            if(submethod=="BLS") {
               if(pruned)
               {
                  typedef binary_line_search<simple_prune<has_gradients_data<chem_opt> > > CC_t;
                  CC_t CC;
                  (chem_opt&) CC=Complex;
                  CC.set_number_of_constraints(nconstraints);
                  typedef binary_gdmc<binary_line_search<simple_prune<has_gradients_data<chem_opt> > > > C_t;
                  C_t C(CC);
                  C.T=T;
                  C.tight_steps=tight_steps;
                  C.max_steps=max_steps;
                  cout << "top=binary_gdmc<binary_line_search<simple_prune<has_gradients_data<chem_opt> > > >\n";
                  C.set_id("top");

                  print<CC_t>(CC);
                  run<C_t>(C,value,value_passed);
               }
               else
               {
                  typedef binary_line_search<noprune<has_gradients_data<chem_opt> > > CC_t;
                  CC_t CC;
                  (chem_opt&) CC=Complex;
                  CC.set_number_of_constraints(nconstraints);
                  typedef binary_gdmc<binary_line_search<noprune<has_gradients_data<chem_opt> > > > C_t;
                  C_t C(CC);
                  C.T=T;
                  C.tight_steps=tight_steps;
                  C.max_steps=max_steps;
                  cout << "top=binary_gdmc<binary_line_search<noprune<has_gradients_data<chem_opt> > > >\n";
                  C.set_id("top");

                  print<CC_t>(CC);
                  run<C_t>(C,value,value_passed);
               }
            }
            else if(submethod=="GBGLS") {

               //Compute the various bases
               long k=0;
               for(long i=0; i< Complex.Substituent_Groups_r.size(); i++)
                  k+=Complex.Substituent_Groups_r[i].allowed_Substituents_r.size();
               refvector<long> bases(k);
               k=0;
               for(long i=0; i< Complex.Substituent_Groups_r.size(); i++)
                  for(long j=0; j<Complex.Substituent_Groups_r[i].allowed_Substituents_r.size(); j++,k++)
                  {
                     bases[k]+=Complex.Substituent_Groups_r[i].allowed_Substituents_r[j].size();
                     if(bases[k]==0) bases[k]=1;
                  }

               if(pruned)
               {
                  typedef reorder_general_base<chem_opt> CC_t;
                  CC_t CC(bases);
                  (chem_opt&) CC=Complex;
                  CC.set_number_of_constraints(nconstraints);
                  CC.at_max=at_max_flag;
                  CC.at_current=at_current_flag;
                  CC.minimax=minimax_flag;
                  gen_base_grad_LS<CC_t,general_base_iterator<chem_opt> > CCC(CC);
                  CCC.precondition_flag=precondition_flag;
                  typedef gen_base_gdmc<gen_base_grad_LS<CC_t,general_base_iterator<chem_opt> > > C_t;
                  C_t C(CCC,bases);
                  C.T=T;
                  C.tight_steps=tight_steps;
                  C.max_steps=max_steps;
                  cout << "top=gen_base_gdmc<gen_base_grad_LS<reorder_general_base<chem_opt> > >\n";
                  C.set_id("top");

                  print<CC_t> (CC);
                  run<C_t>(C,value,value_passed);
               }
               else
               {
                  typedef noprune<chem_opt> CC_t;
                  CC_t CC;
                  (chem_opt&) CC=Complex;
                  CC.minimax = minimax_flag;
                  CC.set_number_of_constraints(nconstraints);
                  gen_base_grad_LS<CC_t,general_base_iterator<chem_opt> > CCC(CC);
                  typedef gen_base_gdmc<gen_base_grad_LS<CC_t,general_base_iterator<chem_opt> > > C_t;
                  C_t C(CCC,bases);
                  C.T=T;
                  C.tight_steps=tight_steps;
                  C.max_steps=max_steps;
                  cout << "top=gen_base_gdmc<gen_base_grad_LS<noprune<chem_opt> > >\n";
                  C.set_id("top");

                  print<CC_t>(CC);
                  run<C_t>(C,value,value_passed);
               }
            }
            else
            {
               cerr << "Please add a submethod\n";
               exit(0);
            }
         }
         else if(method=="GBEN" ||
               method=="gben") {
            cout << "Doing general base entropic optimization" << endl;
            if(max_steps<2) {
               cout << " Please enter maximum number of computations: ";
               cin >> max_steps;
            }
            if(submethod=="BLS") {
               cout << "WARNING: experimental\n";
               if(pruned)
               {
                  typedef binary_line_search<simple_prune<has_gradients_data<chem_opt> > > CC_t;
                  CC_t CC;
                  (chem_opt&) CC=Complex;
                  CC.set_number_of_constraints(nconstraints);
                  typedef binary_entropic<CC_t> C_t;
                  C_t C(CC);
                  C.max_steps=max_steps;
                  cout << "top=binary_entropic<binary_line_search<simple_prune<has_gradients_data<chem_opt> > > >\n";
                  C.set_id("top");

                  print<CC_t>(CC);
                  run<C_t>(C,value,value_passed);
               }
               else
               {
                  typedef binary_line_search<noprune<has_gradients_data<chem_opt> > > CC_t;
                  CC_t CC;
                  (chem_opt&) CC=Complex;
                  CC.set_number_of_constraints(nconstraints);
                  typedef binary_entropic<CC_t> C_t;
                  C_t C(CC);
                  C.max_steps=max_steps;
                  cout << "top=binary_entropic<binary_line_search<noprune<has_gradients_data<chem_opt> > > >\n";
                  C.set_id("top");

                  print<CC_t>(CC);
                  run<C_t>(C,value,value_passed);
               }
            }
            else if(submethod=="GBGLS") {

               //Compute the various bases
               long k=0;
               for(long i=0; i< Complex.Substituent_Groups_r.size(); i++)
                  k+=Complex.Substituent_Groups_r[i].allowed_Substituents_r.size();
               refvector<long> bases(k);
               k=0;
               for(long i=0; i< Complex.Substituent_Groups_r.size(); i++)
                  for(long j=0; j<Complex.Substituent_Groups_r[i].allowed_Substituents_r.size(); j++,k++)
                  {
                     bases[k]+=Complex.Substituent_Groups_r[i].allowed_Substituents_r[j].size();
                     if(bases[k]==0) bases[k]=1;
                  }

               if(pruned)
               {
                  typedef reorder_general_base<chem_opt> CC_t;
                  CC_t CC(bases);
                  (chem_opt&) CC=Complex;
                  CC.set_number_of_constraints(nconstraints);
                  CC.at_max=at_max_flag;
                  CC.at_current=at_current_flag;
                  CC.minimax=minimax_flag;
                  gen_base_grad_LS<CC_t,general_base_iterator<chem_opt> > CCC(CC);
                  CCC.precondition_flag=precondition_flag;
                  typedef gen_base_entropic<gen_base_grad_LS<CC_t,general_base_iterator<chem_opt> > > C_t;
                  C_t C(CCC,bases,gbenreorderflag);
                  C.max_steps=max_steps;
                  C.nruns=tight_steps;
                  cout << "top=gen_base_entropic<gen_base_grad_LS<reorder_general_base<chem_opt> > >\n";
                  C.set_id("top");
                  print<CC_t> (CC);
                  run<C_t>(C,value,value_passed);
               }
               else
               {
                  typedef noprune<chem_opt> CC_t;
                  CC_t CC;
                  (chem_opt&) CC=Complex;
                  CC.minimax = minimax_flag;
                  CC.set_number_of_constraints(nconstraints);
                  gen_base_grad_LS<CC_t,general_base_iterator<chem_opt> > CCC(CC);
                  typedef gen_base_entropic<gen_base_grad_LS<CC_t,general_base_iterator<chem_opt> > > C_t;
                  C_t C(CCC,bases,gbenreorderflag);
                  C.max_steps=max_steps;
                  C.nruns=tight_steps;
                  cout << "top=gen_base_entropic<gen_base_grad_LS<noprune<chem_opt> > >\n";
                  C.set_id("top");

                  print<CC_t>(CC);
                  run<C_t>(C,value,value_passed);
               }
            }
            else if(submethod=="GBLS") {

               //Compute the various bases
               long k=0;
               for(long i=0; i< Complex.Substituent_Groups_r.size(); i++)
                  k+=Complex.Substituent_Groups_r[i].allowed_Substituents_r.size();
               refvector<long> bases(k);
               k=0;
               for(long i=0; i< Complex.Substituent_Groups_r.size(); i++)
                  for(long j=0; j<Complex.Substituent_Groups_r[i].allowed_Substituents_r.size(); j++,k++)
                  {
                     bases[k]+=Complex.Substituent_Groups_r[i].allowed_Substituents_r[j].size();
                     if(bases[k]==0) bases[k]=1;
                  }

               if(pruned)
               {
                  typedef reorder_general_base<chem_opt> CC_t;
                  CC_t CC(bases);
                  CC.set_number_of_constraints(nconstraints);
                  (chem_opt&) CC=Complex;
                  CC.at_max=at_max_flag;
                  CC.at_current=at_current_flag;
                  CC.minimax=minimax_flag;
                  gen_base_LS<CC_t,general_base_iterator<chem_opt> > CCC(CC);
                  typedef gen_base_entropic<gen_base_LS<CC_t,general_base_iterator<chem_opt> > > C_t;
                  C_t C(CCC,bases);
                  C.max_steps=max_steps;
                  C.nruns=tight_steps;
                  cout << "top=gen_base_entropic<gen_base_LS<reorder_general_base<chem_opt> > >\n";
                  C.set_id("top");

                  print<CC_t> (CC);
                  run<C_t>(C,value,value_passed);
               }
               else
               {
                  typedef noprune<chem_opt> CC_t;
                  CC_t CC;
                  (chem_opt&) CC=Complex;
                  CC.minimax = minimax_flag;
                  CC.set_number_of_constraints(nconstraints);
                  gen_base_LS<CC_t, general_base_iterator<chem_opt> > CCC(CC);
                  typedef gen_base_entropic<gen_base_LS<CC_t,general_base_iterator<chem_opt> > > C_t;
                  C_t C(CCC,bases);
                  C.max_steps=max_steps;
                  C.nruns=tight_steps;
                  cout << "top=gen_base_entropic<gen_base_LS<noprune<chem_opt> > >\n";
                  C.set_id("top");

                  print<CC_t>(CC);
                  run<C_t>(C,value,value_passed);
               }
            }
            else
            {
               cerr << "Please add a submethod\n";
               exit(0);
            }
         }
         else if (method=="GBLS")
         {
            //Compute the various bases
            long k=0;
            for(long i=0; i< Complex.Substituent_Groups_r.size(); i++)
               k+=Complex.Substituent_Groups_r[i].allowed_Substituents_r.size();
            refvector<long> bases(k);
            k=0;
            for(long i=0; i< Complex.Substituent_Groups_r.size(); i++)
               for(long j=0; j<Complex.Substituent_Groups_r[i].allowed_Substituents_r.size(); j++,k++)
               {
                  bases[k]+=Complex.Substituent_Groups_r[i].allowed_Substituents_r[j].size();
                  if(bases[k]==0) bases[k]=1;
               }

            if(pruned)
            {
               typedef reorder_general_base<chem_opt> C_t;
               C_t CC(bases);
               CC.set_number_of_constraints(nconstraints);
               (chem_opt&) CC=Complex;
               CC.at_max=at_max_flag;
               CC.at_current=at_current_flag;
               CC.minimax=minimax_flag;
               gen_base_LS<C_t,general_base_iterator<chem_opt> > C(CC);
               cout << "top=gen_base_LS<reorder_general_base<chem_opt> >\n";
               C.set_id("top");
               print<C_t>(CC);
               run<gen_base_LS<C_t,general_base_iterator<chem_opt> > > (C,value,value_passed);
            }
            else
            {
               typedef noprune<chem_opt> C_t;
               C_t CC;
               (chem_opt&) CC=Complex;
               CC.minimax = minimax_flag;
               CC.set_number_of_constraints(nconstraints);
               gen_base_LS<C_t,general_base_iterator<chem_opt> > C(CC);
               cout << "top=gen_base_LS<noprune<chem_opt> >\n";
               C.set_id("top");
               print<C_t >(CC);
               run<gen_base_LS<C_t,general_base_iterator<chem_opt> > >(C,value,value_passed);
            }
         }
         else if (method=="GBGLS")
         {
            //Compute the various bases
            long k=0;
            for(long i=0; i< Complex.Substituent_Groups_r.size(); i++)
               k+=Complex.Substituent_Groups_r[i].allowed_Substituents_r.size();
            refvector<long> bases(k);
            k=0;
            for(long i=0; i< Complex.Substituent_Groups_r.size(); i++)
               for(long j=0; j<Complex.Substituent_Groups_r[i].allowed_Substituents_r.size(); j++,k++)
               {
                  bases[k]+=Complex.Substituent_Groups_r[i].allowed_Substituents_r[j].size();
                  if(bases[k]==0) bases[k]=1;
               }

            if(pruned)
            {
               typedef reorder_general_base<chem_opt> C_t;
               C_t CC(bases);
               (chem_opt&) CC=Complex;
               CC.set_number_of_constraints(nconstraints);
               CC.at_max=at_max_flag;
               CC.at_current=at_current_flag;
               CC.minimax=minimax_flag;
               gen_base_grad_LS<C_t,general_base_iterator<chem_opt> > C(CC);
               cout << "top=gen_base_grad_LS<reorder_general_base<chem_opt> >\n";
               C.set_id("top");
               C.precondition_flag=precondition_flag;
               print<C_t>(CC);
               run<gen_base_grad_LS<C_t,general_base_iterator<chem_opt> > >(C,value,value_passed);
            }
            else
            {
               typedef noprune<chem_opt> C_t;
               C_t CC;
               (chem_opt&) CC=Complex;
               CC.minimax = minimax_flag;
               CC.set_number_of_constraints(nconstraints);
               gen_base_grad_LS<C_t,general_base_iterator<chem_opt> > C(CC);
               cout << "top=gen_base_grad_LS<noprune<chem_opt> >\n";
               C.set_id("top");
               print<C_t >(CC);
               run<gen_base_grad_LS<C_t,general_base_iterator<chem_opt> > >(C,value,value_passed);
            }
         }
         else if(method=="BLS") {
            if(pruned)
            {
               typedef binary_line_search<simple_prune<chem_opt> > C_t;
               C_t C;
               (chem_opt&) C=Complex;
               C.set_number_of_constraints(nconstraints);
               cout << "top=binary_line_search<simple_prune<chem_opt> >\n";
               C.set_id("top");
               print<C_t>(C);
               run<C_t>(C,value,value_passed);
            }
            else
            {
               typedef binary_line_search<noprune<chem_opt> > C_t;
               C_t C;
               (chem_opt&) C=Complex;
               C.set_number_of_constraints(nconstraints);
               C.minimax = minimax_flag;
               cout << "top=binary_line_search<noprune<chem_opt> >\n";
               C.set_id("top");
               print<C_t>(C);
               run<C_t>(C,value,value_passed);
            }
         }
      }
      else
      {
         if(pruned)
         {
            typedef binary_line_search<simple_prune<chem_opt> > C_t;
            C_t C;
            (chem_opt&) C=Complex;
            C.set_number_of_constraints(nconstraints);
            cout << "top=binary_line_search<simple_prune<chem_opt> >\n";
            C.set_id("top");
            print<C_t>(C);
            run<C_t>(C,value,value_passed);
         }
         else
         {
            typedef binary_line_search<noprune<chem_opt> > C_t;
            C_t C;
            (chem_opt&) C=Complex;
            C.minimax = minimax_flag;
            C.set_number_of_constraints(nconstraints);
            cout << "top=binary_line_search<noprune<chem_opt> >\n";
            C.set_id("top");
            print<C_t>(C);
            run<C_t>(C,value,value_passed);
         }
      }
   } catch(exception& e) {
      cerr << e.what() << endl;
      cerr << "Highest level" << endl;
   }
}


