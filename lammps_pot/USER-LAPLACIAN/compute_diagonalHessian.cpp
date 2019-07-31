/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This is a reductoin of code to calculate Hessian matrix, 
   but now this is ONLY diagonal elements of Hessian matrix.
   Output data "hessian" is global vector with number of elements = 3 * N_atoms
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "neighbor.h"
#include "compute_diagonalHessian.h"
#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "update.h"
#include "memory.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "comm.h"
#include "group.h"
#include "fix.h"
#include "min.h"
#include "finish.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDiagonalHessian::ComputeDiagonalHessian(LAMMPS *lmp, int narg, char **arg)
    : Compute(lmp, narg, arg) {
  if (domain->box_exist == 0)
    error->all(FLERR,"diagonalHessian compute before simulation box is defined");
  if (narg != 4)
    error->all(FLERR, "Illegal compute diagonalHessian command");
  
  //parse input arguments 
  natomgroup = group->count(igroup);
  epsilon = atof(arg[3]);
  iepsilon = 1 / epsilon;
  

  /* even though this is a massive 2d array, return the a vector instead.
   * we will explicitly manage the addressing in each dimension with a
   * preprocessor index macro. */
  vector_flag = 1;
  extvector = 0;

  /* these values will change if the system size changes. */
  int ndofs = natomgroup * 3;
  size_vector = ndofs ;    //   * ndofs;   /* later, please, change for diagonal elements ONLY, before:   * ndofs; */ 

  mylocalsize = 0;
  myglobalsize = 0;

  fglobal_ref = fglobal_new_pos = fglobal_new_neg = fglobal_copy = NULL;
  hessian = NULL;
}


/* ---------------------------------------------------------------------- */

ComputeDiagonalHessian::~ComputeDiagonalHessian() {
  free(fglobal_ref);
  free(fglobal_new_pos);
  free(fglobal_new_neg);
  free(fglobal_copy);
  free(hessian);
}

/* ---------------------------------------------------------------------- */

void ComputeDiagonalHessian::compute_vector(void) {
  invoked_vector = update->ntimestep;
  /* tags must be defined and consecutive. */
  if (atom->tag_enable == 0)
    error->all(FLERR,
               "Cannot use Hessian compute unless atoms have IDs");
  if (atom->tag_consecutive() == 0)
    error->all(FLERR,
               "Atom IDs must be consecutive for Hessian compute");

  /* get pointers to all the original data. */
  double **x = atom->x;
  double **f = atom->f;
  /* the global force and hessian arrays must be explicitly the correct size. */
  int needglobalsize = natomgroup;
  int ndofs = natomgroup * 3;
  bigint nhessianelements = ndofs ;    // * ndofs ;     /* later, please,  change for diagonal elements ONLY, before:   * ndofs; */ 
  if (needglobalsize != myglobalsize) {
    free (fglobal_ref); 
    free (fglobal_new_pos); 
    free (fglobal_new_neg); 
    free (fglobal_copy);
    free (hessian);
/*    fglobal_ref = (double *) malloc (ndofs * sizeof (double));   
    fglobal_new_pos = (double *) malloc (ndofs * sizeof (double));   
    fglobal_new_neg = (double *) malloc (ndofs * sizeof (double));   
    fglobal_copy = (double *) malloc (ndofs * sizeof (double));   
*/

    fglobal_ref = (double *) malloc ( 3 * sizeof (double));
    fglobal_new_pos = (double *) malloc ( 3 * sizeof (double));
    fglobal_new_neg = (double *) malloc ( 3 * sizeof (double));
    fglobal_copy = (double *) malloc ( 3 * sizeof (double));


    hessian = (double *) malloc (nhessianelements * sizeof (double));


    /* always be sure to set the output vector since the address keeps changing. */
    vector = hessian;

    myglobalsize = needglobalsize;
  }

  /* a lot of the hessian will be zero, so start there. */
  memset (hessian, 0, nhessianelements * sizeof (double));
  /* set the initial fglobal_new to zero */
/*  memset (&fglobal_new_pos[0], 0, myglobalsize * 3 * sizeof (double));
  memset (&fglobal_new_neg[0], 0, myglobalsize * 3 * sizeof (double));
*/
  memset (&fglobal_new_pos[0], 0, 3 * sizeof (double));
  memset (&fglobal_new_neg[0], 0, 3 * sizeof (double));



  /* set up a map if none exists so we can incrementally loop through all dofs
   * regardless of the location of the atom data. */
  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }

  /* no energy or virial updates. */
  int eflag = 0;
  int vflag = 0;

  /* allow pair and kspace compute to be turned off via modify flags. */
  if (force->pair && force->pair->compute_flag)
    pair_compute_flag = 1;
  else
    pair_compute_flag = 0;
  if (force->kspace && force->kspace->compute_flag)
    kspace_compute_flag = 1;
  else
    kspace_compute_flag = 0;

  /* do a standard force call to get the reference forces. */
  comm->forward_comm();
  force_clear();
  if (modify->n_pre_force) modify->pre_force(vflag);

  if (pair_compute_flag) force->pair->compute(eflag, vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
  }

  if (kspace_compute_flag) force->kspace->compute(eflag, vflag);
  if (force->newton) comm->reverse_comm();
  if (modify->n_post_force) modify->post_force(vflag);

  /* construct fglobal_ref by explicit scatter and reduce to preserve atom-id
   * ordering. */
  int m, reduce_m, id=0;
//  memset (&fglobal_copy[0], 0, myglobalsize * 3 * sizeof (double));
  memset (&fglobal_copy[0], 0, 3 * sizeof (double));

  for (int i = 1; i <= atom->natoms; i++) {
    m = atom->map(i);
    if (atom->mask[m] & groupbit) {
      //reduce_m = atom->tag[m] - 1;
      for (int j = 0; j < domain->dimension; j++){
//        fglobal_copy[idx2_c(id, j,myglobalsize )] = f[m][j];
        fglobal_copy[j] = f[m][j];  
      }
      id++;
    }
  }
//  MPI_Allreduce (fglobal_copy, fglobal_ref, ndofs, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce (fglobal_copy, fglobal_ref, 3, MPI_DOUBLE, MPI_SUM, world);

  /* set up a group of atom that are in the environment*/
  char **group_arg = new char*[4];
  group_arg[0] = (char *) "bg";
  group_arg[1] = (char *) "subtract";
  group_arg[2] = (char *) "all";
  group_arg[3] = (char *) group->names[igroup];
  group->assign(4,group_arg);
  delete [] group_arg; 
  
  int alligroup=group->find("all");
  int bgigroup=group->find("bg");
  
   
  /* do numerical hessian compute by forward differences. */
  int n, reduce_n=-1, index_a, index_b, global_atom_a=-1, global_atom_b=-1;
  double mass, difference, mass_weight, xstore;
  for (int i = 1; i <= atom->natoms; i++) {//for **
    m = atom->map(i);
    if (atom->mask[m]& groupbit) { // if **
      global_atom_a++;
     
      for (int j = 0; j < domain->dimension; j++) { //loop over dimensin j
       /* increment the dof by epsilon on the right task. */
         if (atom->mask[m]& groupbit) {
           xstore = x[m][j];
           x[m][j] += epsilon;
         }

        //  
	/* standard force call. */
        comm->forward_comm();
        force_clear();
        if (modify->n_pre_force) modify->pre_force(vflag);

        if (pair_compute_flag) force->pair->compute(eflag, vflag);

        if (atom->molecular) {
         	if (force->bond) force->bond->compute(eflag, vflag);
  	        if (force->angle) force->angle->compute(eflag, vflag);
  	        if (force->dihedral) force->dihedral->compute(eflag, vflag);
        	if (force->improper) force->improper->compute(eflag, vflag);
        }

        if (kspace_compute_flag) force->kspace->compute(eflag, vflag);

        if (force->newton) comm->reverse_comm();
        if (modify->n_post_force) modify->post_force(vflag);

 
       /* put the original position back. */
       if (atom->mask[m] & groupbit) x[m][j] = xstore;

       /* construct fglobal_new by explicit scatter and reduce to preserve
        * atom-id ordering. */
       memset (&fglobal_copy[0], 0, 3 * sizeof (double));

//       memset (&fglobal_copy[0], 0, myglobalsize * 3 * sizeof (double));
//       for (int k = i; k <= i; k++ )  {    //   1; k <= atom->natoms; k++) {

//         n = atom->map(k);
         if (atom->mask[m] & groupbit) {
           reduce_n++;
//           for (int l = j; l <= j ; l++)      //    0; l < domain->dimension; l++)
//             fglobal_copy[idx2_c(reduce_n, j, myglobalsize)] = f[m][j];
           fglobal_copy[j] = f[m][j];

         }
//       }

       reduce_n=-1;
//       MPI_Allreduce (fglobal_copy, fglobal_new_pos, ndofs, MPI_DOUBLE, MPI_SUM, world);
       MPI_Allreduce (fglobal_copy, fglobal_new_pos, 3, MPI_DOUBLE, MPI_SUM, world);	
      /* increment the dof by epsilon on the left task. */
      if (atom->mask[m] & groupbit) {
        xstore = x[m][j];
        x[m][j] -= epsilon;
      }

      /* standard force call. */
      //standard_force_call(pair_compute_flag,kspace_compute_flag,eflag,vflag);
      comm->forward_comm();
      force_clear();
      if (modify->n_pre_force) modify->pre_force(vflag);

      if (pair_compute_flag) force->pair->compute(eflag, vflag);

       if (atom->molecular) {
         if (force->bond) force->bond->compute(eflag, vflag);
         if (force->angle) force->angle->compute(eflag, vflag);
          if (force->dihedral) force->dihedral->compute(eflag, vflag);
         if (force->improper) force->improper->compute(eflag, vflag);
       }

      if (kspace_compute_flag) force->kspace->compute(eflag, vflag);

      if (force->newton) comm->reverse_comm();
      if (modify->n_post_force) modify->post_force(vflag);

      /* put the original position back. */
      if (atom->mask[m] & groupbit) x[m][j] = xstore;
      
      /* construct fglobal_new by explicit scatter and reduce to preserve
       * atom-id ordering. */
      memset (&fglobal_copy[0], 0, 3 * sizeof (double));
//      for (int k = i; k <= i; k++ ) {   //  1; k <= atom->natoms; k++) {
 
//        n = atom->map(k);
        if (atom->mask[m] & groupbit) {
          reduce_n++;
//          for (int l = j; l <= j ; l++)      //    0; l < domain->dimension; l++)
//            fglobal_copy[idx2_c(reduce_n, j, myglobalsize)] = f[m][j];
          fglobal_copy[j] = f[m][j];
        }
//      }
      reduce_n=-1;
//      MPI_Allreduce (fglobal_copy, fglobal_new_neg, ndofs, MPI_DOUBLE, MPI_SUM, world);
      MPI_Allreduce (fglobal_copy, fglobal_new_neg, 3, MPI_DOUBLE, MPI_SUM, world);
      /* compute the difference (not using symmetry so we can do an in-place
       * reduciton). */
//      index_a = j ;       //     + 3 * global_atom_a;
//      for (int k = i; k <= i; k++ ) {   //   1; k <= atom->natoms; k++) {
//        n = atom->map(k);
        if (atom->mask[m]& groupbit) {

//          global_atom_b++;
          /* once again, global arrays use 1-based indexing, so have to rebase
           * them to 0. */
//          for (int l = j; l <= j ; l++) {     //    0; l < domain->dimension; l++) {
//            index_b =  global_atom_a ;  //      l + 3 * global_atom_b;
            difference =
                fglobal_new_pos[j] - \
                fglobal_new_neg[j];
//                fglobal_new_pos[idx2_c(global_atom_b, j, myglobalsize)] - \
//                fglobal_new_neg[idx2_c(global_atom_b, j, myglobalsize)];
   //       hessian[idx2_c(index_a, index_b, ndofs)] =

            hessian[idx2_c(j , global_atom_a , 3)] =
               -0.5 * difference * iepsilon;
//          }
        }
//      }

    global_atom_b=-1;
    

  
    }//end for over dimension j
   
   }//end if **
  }// end for **

global_atom_a=-1;
  /* only reduce the hessian to the root task. */
  MPI_Reduce(MPI_IN_PLACE, hessian, nhessianelements, MPI_DOUBLE, MPI_SUM, 0, world);
  
  /* destroy the atom map. */
  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  }

}

void ComputeDiagonalHessian::force_clear() {
  size_t nbytes;
  int nlocal = atom->nlocal;

  nbytes = sizeof (double) * nlocal;
  if (force->newton) nbytes += sizeof (double) * atom->nghost;

  if (nbytes) memset (&atom->f[0][0], 0, 3 * nbytes);
}
/*void ComputeDiagonalHessian::standard_force_call(int pair_compute_flag, int kspace_compute_flag, int eflag, int vflag){

  comm->forward_comm();
  force_clear();
  if (modify->n_pre_force) modify->pre_force(vflag);

  if (pair_compute_flag) force->pair->compute(eflag, vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
  }

  if (kspace_compute_flag) force->kspace->compute(eflag, vflag);
  if (force->newton) comm->reverse_comm();
  if (modify->n_post_force) modify->post_force(vflag);


}*/
