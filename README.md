# Introduction
This is the ARL Discrete Chemical Compound Space Optimization (ARL DCCSO) project.
ARL Discrete Chemical Compound Space Optimization (ARL DCCSO) is governed by the terms of the  Creative Commons Zero 1.0 Universal (CC0 1.0) Public Domain Dedication (the Agreement). You should have received a copy of the Agreement with a copy of this software. If not, see https://github.com/USArmyResearchLab/ARLDCCSO.

Your use or distribution of ARL Discrete Chemical Compound Space Optimization, in both source and binary form, in whole or in part, implies your agreement to abide by the terms set forth in the Agreement in full.

ARL Discrete Chemical Compound Space Optimization is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the  Creative Commons Zero 1.0 Universal (CC0 1.0) Public Domain Dedication for more details.

You may find the full license in the file LICENSE in this directory.

ARL DCCSO provides a framework to maximize properties of molecules using substitution patterns under constraints.

## Contributions

Due to legal issues, every contributor will need to have a signed Contributor License Agreement on file.
The ARL Contributor License Agreement (ARL Form 266) can be found [here](https://github.com/USArmyResearchLab/ARL-Open-Source-Guidance-and-Instructions/blob/develop/ARL%20Form%20-%20266.pdf). 
Each external contributor must execute and return a copy for each project that he or she intends to contribute to. 
Once ARL receives the executed form, it will remain in force permanently. 
Thus, external contributors need only execute the form once for each project that they plan on contributing to.
Contributions to ARL DCCSO will undergo review and will be released under the [Apache 2.0 license](http://apache.org/licenses/LICENSE-2.0).

# Building the Executable
## Prerequisites

Before compiling, BCR_CPP_LA header files are required. They may be downloaded at https://gitlab.com/crinders/BCR_CPP_LA

## Compilation

To compile, modify the makefiles in the `Build` subdirectory to reflect your directory structure and compiler needs and enter
`
make all
`
The executable will be called `DiscreteCCSOpt` unless you change it.

To compile the debug version of the code, do the same in the `Debug` subdirectory.
The debug binary will be called `DiscreteCCSOpt-d`

# Documentation
In order to generate the documentation use doxygen in the root directory:
`doxygen`
The documentation will then reside in the html/ and latex/ directories.


# Usage

Basic usage is:
`DiscreteCCSOpt`

The program will then prompt you for an input file. For specifics on the
input file format and other requirements, please refer to the documentation or 
the source files directly (which also include the documentation).
# ARL_Discrete_Chemical_Compound_Space_Optimization
