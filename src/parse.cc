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
#include "parse.hh"

//! Parse the input file for sections.
refvector<section> parse(istream& in)
/*!
 * Each section starts with a keyword(potentially empty) and an open parenthesis.
 * Each section is ended by a parenthesis.
 *
 */
{
   string serr="parse(istream& in):";
   string callederr=serr+"called";
   string interim;
   refvector<section> r;
   // remove white space and comments
   try {
      string s;
      bool atend=false;
      stringstream o;
      o.str("");
      while(!atend)
      {
         char c;
         in >> skipws >> c;
         // This is to allow comments
         if (c == '#') {
            c = in.get();
            while(c!='\n'&& c!='#' && !in.eof())
               c = in.get();
         }
         else {
            if(c!='\n' && c!='\r')
               o << c;
         }
         c=in.get();
         atend=!in.good();
         in.unget();
      }
      //! Now parse for sections.
      atend=false;
      char c;
      while (!atend) {
         stringstream p;
         p.str("");
         string interim="";
         int level=0;
         getline(o,interim,'(');
         level++;
         while(level>0)
         {
            c=o.get();
            if(!o.good()) {
               stringstream levelstring;
               levelstring << level;
               throw domain_error(serr+" file ended before section "+interim+" completed. Level = "+levelstring.str());
            }
            if(c==')') level--;
            if(c=='(') level++;
            p << c;
         }
         cout << "Found section #" << r.size() << ":" << interim << " = " << p.str() << endl;
         r.push_back(section(interim,p.str()));
         c = o.get();
         atend=!o.good() || o.eof();
         cout << atend << "\n";
         o.unget();
      }
   } catch(exception& e) {
      cerr << e.what() << endl;
      throw domain_error("ChemGroup::ChemGroup(istream& in)");
   }
   return r;
}


