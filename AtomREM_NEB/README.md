AtomREM and NEB
====================

[![License](https://img.shields.io/badge/License-GPLv3-orange.svg)] 

This is an explanation about the procedure of simulation using AtomREM and NEB methods to get a minimum energy path with DFT accuracy.


PROCEDURE
---

  **I) Fix stable state:**  First, we need to fix stable state and calculate cell parameters and atomic positions using DFT calculation.
For example, in /exmpl\_PW/ directory the script to calculate stable state for SH3 system under pressure 200 GPa for 16 atoms. 
To get the possibility for atoms to move according to a reaction, we defined the cell with the free structure in quantum espresso [1-3]

  **II) Preparation for DPMD potential:** To use AtomREM we may use classical potentials or Deep Potential molecular dynamics [4,5].  
To get the best accuracy and make the next step with a minimum inconsistency, we should use DPMD potential. 
The preparation of DPMD potential takes time and computational resources, but it is better because it also uses DFT calculation with 
approximation on deep learning neural network. 

  **III) AtomREM simulation:** To use AtomREM, kindly see instructions in the work [6]. After the simulation, we have possible reaction paths, which go to the initial and new stable states.
We should find a new one and copy files "atomic\_coord\_q.lammpstrj"  and "coordinates\_rank\_100XXX.lammpstrj" to a working directory. 


  **IV) Relaxation of the new structure:** Using DFT calculation make the relaxation of the new structure and find atomic positions together with new cell parameters.


  **V) Generate atomic positions for NEB calculation:** You can use "extract\_images.sh" script to generate images for NEB calculation.
Please, pay attantion, you should to define paths to two files: "atomic\_coord\_q.lammpstrj"  and "coordinates\_rank\_100XXX.lammpstrj" in the script.
After that, kindly use the script "make\_IN.sh". As a source file to make an input file for NEB calculation script uses the "H3S\_initial.in" file.  


  **VI) Preparation of NEB calculation:** The first one is using the cell parameters for a stable state. The second simulation is using new cell parameters for a new structure.
Here we use 18 images because we will use 18 nodes to accelerate the calculations. You can change this number as you need. 
For each calculation, we can use an acceleration procedure like: 

 - to calculate the NEB reaction path using short cutoff, for example, 20 Ry;

 - to repeat the NEB calculation from the last step, but using long cutoff, for example, 80 Ry (using "H3S\_initial\_2.in" file).

 To realize this procedure, you should change the input script for quantum espresso like:

 - generate new positions in input with "gen\_second\_neb.sh" script, which extracts the last images and add them to input script;

 - change the cutoff for larger value, as 80 Ry or more;

 - add in the quantum espresso input script: 

&ELECTRONS
  startingpot = 'file',
  startingwfc = 'file',

and start again from the same directory. Of course, the output directory should be the same as working to save all data to continue the simulation.



  **VII) Calculate the relaxation of the cell at each image of reaction path from NEB:** 
 The reason why we add relaxation step AFTER NEB is: 

 - NEB calculation uses relaxation for the atomic positions of images, but with the fix cell parameters.

 That is why after NEB we should make relaxation of the cell parameters for each image independently (with the fix atomic positions). 

To generate jobs for PW relaxation with the fix atomic positions you can use script "gen\_PW\_relax.sh", which generates input files 
for quantum espresso (from standard input "H3S.in") with the different atomic coordinates (extracted from NEB calculation). 
All calculations are saved to different folders, like IM\_1/, IM\_2/ , etc.

  **VIII) Calculate minimum energy path:**  To see energy minimum path we should to extract information about energies from PW relaxation
for each NEB calculation (with two different cell parameters). To extract energies, you can use script "extract\_minimum\_energies.sh".


NOTE: 
---

* the whole.sh script is an example to submit job for HPC system like ISSP supercomputer in University of Tokyo (based on SUSE Linux operation system), 
which is included all steps for NEB calculations and PW relaxation. 

* all scripts are in the folder /SCRIPTS/. 

* a calculation is needed the high computational cost, so this is the example for 18 nodes with 24 cores in each (432 processors in total).



References
---

1. QUANTUM ESPRESSO is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling 
at the nanoscale, which is based on density-functional theory, plane waves, and pseudopotentials. https://www.quantum-espresso.org
2. P. Giannozzi, S. Baroni, N. Bonini, et al. J.Phys.: Condens.Matter 21, 395502 (2009).  
3. P. Giannozzi, O. Andreussi, T. Brumme, et al. J.Phys.: Condens.Matter 29, 465901 (2017).
4. L. Zhang, J. Han, H. Wang, R. Car, Phys. Rev. Lett. 120 (2018) 143001. doi:10.1103/PhysRevLett.120.143001. https://link.aps.org/doi/10.1103/PhysRevLett.120.143001 
5. A deep learning package for many-body potential energy representation and molecular dynamics. https://github.com/deepmodeling/deepmd-kit
6. Yu.S. Nagornov, R. Akashi, arXiv:1907.13316. https://arxiv.org/abs/1907.13316


