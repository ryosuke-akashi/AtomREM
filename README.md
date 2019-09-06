[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/nagornovys/Cancer_cell_evolution/blob/master/LICENSE)

# AtomREM - Atom*istic* R*are* E*vent* M*anager*

---

**Nature of problem:** AtomREM has been developed to help to find a reaction path on the atomic level using weighted Langevin mechanics on the inverse potential landscape. The method is describe the low temperature transformation of complex systems without artificial forces or/and collective variables.

**Solution method:** Recently we designed a non-empirical scheme to search for the minimum-energy escape paths from the minima of the potential surface to unknown saddle points nearby with one dimensional application [J. Phys. Soc. Jpn. 87, 063801 (2018)]. The method is based on the Master equation and its solver algorithm is constructed to move the walkers up the surface through the potential valleys [Physica A, 528, 121481 (2019)]. The stochastic algorithm uses a birth/death stochastic processes for numerous of walker in combination with the Langevin equation. Each walker obey the statistics of average movement on the reaction path under biasing potential. Under this consideration the reaction path is derived as the average of the walker distribution from stable to saddle states.

---

**There is a _THEORY_ of method:**
                  J. Phys. Soc. Jpn. 87, 063801 (2018) - https://journals.jps.jp/doi/abs/10.7566/JPSJ.87.063801 

**There is an _ALGORITHM_ of method and examples of simulations:** 
                  Yu.S.Nagornov, R.Akashi Physica A, 528, 121481 (2019) - https://doi.org/10.1016/j.physa.2019.121481 

---

## LICENCE
GNU General Public License 3 (GPL)

## CITATION:
Publications making substantial use of AtomREM  or diagonalHESSIAN package for LAMMPS  should cite this software package and the following paper:
Yu.S.Nagornov, R.Akashi  "AtomREM: Non-empirical seeker of the minimum energy escape paths on many-dimensional potential landscapes without coarse graining"    https://arxiv.org/abs/1907.13316

## Examples of simulations 
Please, see the exampples of simulations:
https://www.youtube.com/watch?v=UGrnAlYLJY8&list=PLuI0M67EVLYlDUeZqBTsoju9sslH_Mdp2&index=1

## Instllation and detailed description of the program usage
Installation procedure and detailed description of usage of the program is in the manuscript by Yuri S. Nagornov, Ryosuke Akashi "AtomREM: Non-empirical seeker of the minimum energy escape paths on many-dimensional potential landscapes without coarse graining" at https://arxiv.org/abs/1907.13316

**NOTE**
The installation requires 600MB disk space, as examples/ directory includes several output examples for reference.

**Installation**
1) git clone https://github.com/ryosuke-akashi/AtomREM
2) Please read reference (https://arxiv.org/abs/1907.13316) for further installation steps.

---

## CONTENT OF PACKAGE:

 - lammps_pot - the AtomREM program code to find reaction pathways using the LAMMPS calculations of the atomic potentials, forces and laplacian (diagonal elements of Hessian matrix). The LAMMPS package is called as a library (or shared library). This version allows to use any potentials, which are available in LAMMPS. 

- analytic_pot - the version of the code with analytical calculation of atomic potentials, forces and laplacian with potential for Lennard-Jones function.

- examples - the examples of simulations of reaction paths for argon cluster with 7,13,38 atoms, for argon crystal and cyclobutane molecule. 

## Acknowledgment 

We thank **Taichi Kosugi** for advice on the coding. This research was supported by **MEXT as Exploratory Challenge on Post-K computer** _(Frontiers of Basic Science: Challenging the Limits)_. This research used computational resources of the K computer provided by the RIKEN Advanced Institute for Computational Science through the HPCI System Research project _(Project ID:hp160257, hp170244, hp180184, hp190176)_.

