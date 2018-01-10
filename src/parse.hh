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
#ifndef PARSE_HH_
#define PARSE_HH_
#include <BCR_CPP_LA/refcount.h>
using namespace linear_algebra;

//! Collects information on sections during parsing.
struct section {
private:
   string Name;
   string buffer;
public:
   section(const string& N, const string& b):
      Name(N),buffer(b) {};
   section():
      Name(""), buffer("") {};
   const string& get_Name() const {return Name;}
   const string& get_buffer() const {return buffer;}
   bool set_Name(const string& N) { Name=N; return true;}
   bool set_buffer(const string& b) {buffer=b; return true;}
};

//! Parse the input file for sections.
refvector<section> parse(istream& in);




#endif /* PARSE_HH_ */
