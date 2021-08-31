# SHS_simulation
SHS_simulation is a collection of different programs which can compute the Second Harmonic Intensity Scattered by a supramolecular structure.
Don't forget to cite this work:              
J. Chem. Inf. Model. 2020, 60, 12, pp 5912–5917

**********************************************************************************
#    PySHS V2 - An open source software about Second Harmonic Scattering 
#   developed at Institut Charles Gerhardt Montpellier - ENSCM
#   Dr Pierre-Marie GASSIN
#   Dr Gassin GASSIN
#    (june 2021)  
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
**********************************************************************************


This file describes the layout of the PySHS distribution directory along 
with instruction on how to run examples and set up your first calculation.

A C compiler has to be installed (on macOS or linux, gcc is normally already install, to check it, write in a command line: gcc --version  ).

Python 3.7 or later has to be installed for the script. See https://www.python.org for any question about python installation. 

After unpacking the pySHS distribution, the parent PySHS directory will contain
the following subdirectories:

Documentation
Examples 
Src
Work


It will also contain this file and license.txt. The license.txt file
describes the Gnu Public License under which pySHS is released. By
using pySHS, you agree to comply with the terms of the
license. Please read it carefully before you use the code. 

The following directories are included in the distribution:

Documentation
------------
Contains the documentation for the code, user_guide.pdf.
You can find most answers to your questions here. 

Examples
--------
Contains some examples of calculation and in particular the input and output files

Src
---
Contains the pySHS source code and some script used to generate the input files.

Work
----
It is the directory where the user can store the input and output file. 


To Compile the program:
-----------------------

To compile the SHS program enter in a command line:  
gcc -o Work/SHS Src/SHS.c -lm

This command will create the SHS executable file in the Work directory.

Do the same thing for the sphere_SHS program:

gcc -o Work/sphere_SHS Src/sphere_SHS.c -lm


Running a pySHS calculation
----------------------------
Depending of the pySHS program used, it requires one (sphere_SHS) or two input files (SHS) in order to run. 

 to launch, the SHS program, enter:
./Work/SHS <Keyword> <inputbeta> <input orientation> <outputfile>

where <Keyword> is equal to "polarplot_90" or "polarplot_180" or "angle_scattering" depending on the computation you want to perform.
inputbeta, input_orientation and outputfile are respectively the name of the input file containing the hyperpolarizability of the molecule, the  name of the input file containing the orientation and position of the molecules and the output file.


to launch, the spher_SHS program, enter:

./Work/sphere_SHS <Keyword> <inputbeta> <outputfile>

where <Keyword> is equal to "polarplot_90" or "polarplot_180" or "angle_scattering" depending on the computation you want to perform.
inputbeta, and outputfile are respectively the name of the input file containing the hyperpolarizability of the molecule,  and the output file. 

