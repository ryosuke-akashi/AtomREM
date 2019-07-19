/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* libwrapper = fortran wrappers for LAMMPS library functions.
   See README for compilation instructions */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "stdint.h"
#include "library.h"        /* this is a LAMMPS include file */

/* wrapper for creating a lammps instance from fortran.
   since fortran has no simple way to emit a C-compatible
   argument array, we don't support it. for simplicity,
   the address of the pointer to the lammps object is
   stored in a 64-bit integer on all platforms. */

void lammps_open_(MPI_Fint *comm, int64_t *ptr)
{
    void *obj;
    MPI_Comm ccomm;

    /* convert MPI communicator from fortran to c */
    ccomm = MPI_Comm_f2c(*comm);

//    lammps_open(0,NULL,ccomm,&obj);
//    To stop output from lammps and avoid huge out files
    char* opt[]={"", "-sc", "none"};
 
    lammps_open(3,opt,ccomm,&obj);

    *ptr = (int64_t) obj;
}

/* no-MPI version of the wrapper from above. */

void lammps_open_no_mpi_(int64_t *ptr)
{
    void *obj;

    lammps_open_no_mpi(0,NULL,&obj);
    *ptr = (int64_t) obj;
}

/* wrapper for shutting down a lammps instance from fortran. */

void lammps_close_(int64_t *ptr)
{
    void *obj;
    obj = (void *) *ptr;

    lammps_close(obj);
}

/* wrapper for passing an input file to lammps from fortran.
   since fortran strings are not zero terminated, we have
   to pass the length explicitly and make a copy that is. */

void lammps_file_(int64_t *ptr, char *fname, MPI_Fint *len)
{
    void *obj;
    char *cpy;

    obj = (void *) *ptr;

    cpy = (char *)calloc(*len + 1,sizeof(char));
    memcpy(cpy,fname,*len);

    lammps_file(obj,cpy);
    free(cpy);
}

/* wrapper for passing a line input to lammps from fortran.
   since fortran strings are not zero terminated, we have
   to pass the length explicitly and make a copy that is. */

void lammps_command_(int64_t *ptr, char *line, MPI_Fint *len)
{
    void *obj;
    char *cpy;

    obj = (void *) *ptr;

    cpy = (char *)calloc(*len + 1,sizeof(char));
    memcpy(cpy,line,*len);

    lammps_command(obj,cpy);
    free(cpy);
}

/* fortran wrapper to get the number of atoms from lammps.
   return values require an interface in fortran, so we 
   make the wrapper into a procedure. */

void lammps_get_natoms_(int64_t *ptr, MPI_Fint *natoms)
{
    void *obj;
    obj = (void *) *ptr;

    *natoms = lammps_get_natoms(obj);
}

void lammps_gather_atoms_(int64_t *ptr, char *name, MPI_Fint *len, double *x)
{
    void *obj;
    char *cpy;

    obj = (void *) *ptr;

    cpy = (char *)calloc(*len + 1,sizeof(char));
    memcpy(cpy, name, *len);
    lammps_gather_atoms(obj, cpy, 1, 3, x);
    free(cpy);
}

void lammps_scatter_atoms_(int64_t *ptr, char *name, MPI_Fint *len, double *x)
{
    void *obj;
    char *cpy;

    obj = (void *) *ptr;

    cpy = (char *)calloc(*len + 1,sizeof(char));
    memcpy(cpy, name, *len);
    lammps_scatter_atoms(obj, cpy, 1, 3, x);
    free(cpy);
}

void lammps_extract_compute_(int64_t *ptr, char *name, MPI_Fint *len1, double *pe )
{
    void *obj;
    char *cpy1;

    obj = (void *) *ptr;

    cpy1 = (char *)calloc(*len1 + 1,sizeof(char));
    memcpy(cpy1, name, *len1);

//    printf("obj= %s \n", obj);
    *pe = *(double *)lammps_extract_compute(obj, cpy1, 0, 0);
//    printf("pe= %f \n", *pe);
    free(cpy1);
}

void lammps_extract_compute_atom_vec_(int64_t *ptr, char *name, MPI_Fint *len1, double *pe, MPI_Fint *len2)
{
    void *obj;
    char *cpy1;
    double *vec;
    int i;

    obj = (void *) *ptr;

    cpy1 = (char *)calloc(*len1 + 1,sizeof(char));
    memcpy(cpy1, name, *len1);

    


    vec = (double *)lammps_extract_compute(obj, cpy1, 1, 1);
//    printf("%s\n",cpy1);
//    printf("obj= %s \n", obj);
//   printf("X()=%d\n",vec); 
    for(i = 0; i < *len2; i++){
     *(pe + i) = *(vec + i);
    }
//    printf("pe= %f \n", *pe);
    free(cpy1);
}


void lammps_extract_variable_atom_vec_(int64_t *ptr, char *name, MPI_Fint *len1, double *pe, MPI_Fint *len2)
{   
    void *obj;
    char *cpy1;
    double *vec;
    int i;
    
    obj = (void *) *ptr;
    
    cpy1 = (char *)calloc(*len1 + 1,sizeof(char));
    memcpy(cpy1, name, *len1);

    
    vec = (double *)lammps_extract_variable(obj, cpy1, "all");
    printf("%s\n",cpy1); 
    printf("obj= %s \n", obj);
   printf("X()=%d\n",vec); 
    for(i = 0; i < *len2; i++){
     *(pe + i) = *(vec + i);
    }

        free(cpy1);
        }



void lammps_extract_compute_hessian_(int64_t *ptr, char *name, MPI_Fint *len1, double *hsX, double *hsY, double *hsZ,  MPI_Fint *len2)
{
    void *obj;
    char *cpy1;
//  Hessian diagonal elements as a global vector length = 3 * N_atoms
    double *hes;

    int i;

//    arr = (double *)calloc(*len2 * 3, sizeof(double));
//    arr = (double **)calloc(*len2 * 3, sizeof(double));

    obj = (void *) *ptr;

    cpy1 = (char *)calloc(*len1 + 1,sizeof(char));
    memcpy(cpy1, name, *len1);

//    printf("obj= %s \n", obj);
//    printf("%s\n",cpy1);
    hes = (double *)lammps_extract_compute(obj, cpy1, 0, 1);

//    printf("obj= %s \n", obj);
   
    for(i = 0; i < *len2; i++){
//       printf("X(%d)=%d\n",i,hes);

        *(hsX + i) = *(hes + 3*i + 0 );
        *(hsY + i) = *(hes + 3*i + 1 );
        *(hsZ + i) = *(hes + 3*i + 2 );

    }
      free(cpy1);
}





/* wrapper to copy coordinates from lammps to fortran */

/* NOTE: this is now out-of-date, needs to be updated to lammps_gather_atoms()

void lammps_get_coords_(int64_t *ptr, double *coords)
{
    void *obj;
    obj = (void *) *ptr;

    lammps_get_coords(obj,coords);
}

*/

/* wrapper to copy coordinates from fortran to lammps */

/* NOTE: this is now out-of-date, needs to be updated to lammps_scatter_atoms()

void lammps_put_coords_(int64_t *ptr, double *coords)
{
    void *obj;
    obj = (void *) *ptr;

    lammps_put_coords(obj,coords);
}

*/
