/* alloc.c
 * 
 * Allocate, free and zero fields.
 */
#include "hydro.h"



double ***make_field(hydro_params p) {

   
  double *true_field = malloc((p.slicex+2)*(p.slicey+2)*(p.Lz)*sizeof(double));


  double ***field = (double ***)malloc((p.slicex+2)*sizeof(double **));
  int x, y;

  for(x=0;x<(p.slicex+2);x++) {
    field[x] = (double **)malloc((p.slicey+2)*sizeof(double *));
    for(y=0;y<(p.slicey+2);y++) {
      field[x][y] = &true_field[x*(p.slicey+2)*(p.Lz) + y*(p.Lz)];
    }
  }


  //  free(field[0][0]);
  field[0][0] = true_field;

  return field;
}




double ****make_vector(hydro_params p) {

   
  double *true_field = malloc(3*(p.slicex+2)*(p.slicey+2)
			      *(p.Lz)*sizeof(double));

  int x, y, i;

  double ****vector = (double ****) malloc(3*sizeof(double***));

  for(i=0;i<3;i++) {
    vector[i] = (double ***)malloc((p.slicex+2)*sizeof(double **));

    for(x=0;x<(p.slicex+2);x++) {
      vector[i][x] = (double **)malloc((p.slicey+2)*sizeof(double *));
      for(y=0;y<(p.slicey+2);y++) {
	vector[i][x][y]
	  = &true_field[i*(p.slicex+2)*(p.slicey+2)*(p.Lz) + x*(p.slicey+2)*(p.Lz) + y*(p.Lz)];
      }
    }
  }


  //  free(field[0][0]);
  vector[0][0][0] = true_field;
  
  
  return vector;
}




double ****make_tensor(hydro_params p) {


   
  double *true_field = malloc(TENSOR_CPTS*(p.slicex+2)*(p.slicey+2)*(p.Lz)*sizeof(double));

  int x, y, i;

  double ****tensor = (double ****) malloc(TENSOR_CPTS*sizeof(double***));

  for(i=0;i<TENSOR_CPTS;i++) {
    tensor[i] = (double ***)malloc((p.slicex+2)*sizeof(double **));

    for(x=0;x<(p.slicex+2);x++) {
      tensor[i][x] = (double **)malloc((p.slicey+2)*sizeof(double *));
      for(y=0;y<(p.slicey+2);y++) {
	tensor[i][x][y]
	  = &true_field[i*(p.slicex+2)*(p.slicey+2)*(p.Lz) + x*(p.slicey+2)*(p.Lz) + y*(p.Lz)];
      }
    }
  }


  //  free(field[0][0]);
  tensor[0][0][0] = true_field;
  
  
  return tensor;
}







// not convinced
void free_field(hydro_params p, double ***field) {


  int x, y, z;
  
  

  //  double *true_field = field[0][0];

  free(field[0][0]);

  //  fprintf(stderr,"Freeing true field at %p\n", true_field);

  for(x=0;x<(p.slicex+2);x++) {
      free(field[x]);

  }
  


  free(field);
}



void free_vector(hydro_params p, double ****vector) {


  int x, y, z, i;
  
  

  free(vector[0][0][0]);

  //  fprintf(stderr,"Freeing true field at %p\n", true_field);

  for(i=0;i<3;i++) {
    for(x=0;x<(p.slicex+2);x++) {
      free(vector[i][x]);
    }
    free(vector[i]);
  }
  


  free(vector);
}


void free_tensor(hydro_params p, double ****tensor) {


  int x, y, z, i;
  
  

  free(tensor[0][0][0]);

  //  fprintf(stderr,"Freeing true field at %p\n", true_field);

  for(i=0;i<TENSOR_CPTS;i++) {
    for(x=0;x<(p.slicex+2);x++) {
      free(tensor[i][x]);
    }
    free(tensor[i]);
  }
  


  free(tensor);
}





void alloc_fields(hydro_fields *f, hydro_params p) {

  // fields below also set up in initial condition functions
  f->phi = make_field(p);
  f->pifull = make_field(p);
  f->T = make_field(p);
  f->E = make_field(p);
  f->W = make_field(p);
  
  // pi gets initialised in backstep
  f->pi = make_field(p);

  // phiold initialised by evolve_field
  f->phiold = make_field(p);

  // equation of state (obtained by eq_of_state first time
  // and used in hydro)
  f->kappa = make_field(p);




  // pressure (obtained by eq_of_state first time
  // and used in hydro)
  f->p = make_field(p);

  f->V = make_vector(p); 

  /* (double ****)malloc(3*sizeof(double ***));

  f->V[0] = make_field(p);
  f->V[1] = make_field(p);
  f->V[2] = make_field(p);
  */


  /*  
  f->Z = (double ****)malloc(3*sizeof(double ***));

  f->Z[0] = make_field(p);
  f->Z[1] = make_field(p);
  f->Z[2] = make_field(p);
  */

  f->Z = make_vector(p);

  /*
  f->U = (double ****)malloc(3*sizeof(double ***));

  f->U[0] = make_field(p);
  f->U[1] = make_field(p);
  f->U[2] = make_field(p);
  */

  f->U = make_vector(p);

  /*
  f->F = (double ****)malloc(3*sizeof(double ***));

  f->F[0] = make_field(p);
  f->F[1] = make_field(p);
  f->F[2] = make_field(p);
  */

  f->F = make_vector(p);

  f->uij = make_tensor(p);
  f->udotij = make_tensor(p);

  // (calloc considered harmful)
}





void zero_fields(hydro_fields f, hydro_params p) {

  int x;
  int i;

  // We don't need to do this because create_gaussian_bubble
  // initialises everything, but it keeps valgrind quiet.
  for(x=0;x<p.fieldN;x++) {

    f.phi[0][0][x] = 0.0000;
    f.pifull[0][0][x] = 0.0000;
    f.T[0][0][x] = 0.0000;
    f.E[0][0][x] = 0.0000;
    f.Z[0][0][0][x] = 0.0000;
    f.Z[1][0][0][x] = 0.0000;
    f.Z[2][0][0][x] = 0.0000;
    f.U[0][0][0][x] = 0.0000;
    f.U[1][0][0][x] = 0.0000;
    f.U[2][0][0][x] = 0.0000;
    f.V[0][0][0][x] = 0.0000;
    f.V[1][0][0][x] = 0.0000;
    f.V[2][0][0][x] = 0.0000;
    f.W[0][0][x] = 0.0000;

    f.pi[0][0][x] = 0.0000;
    f.phiold[0][0][x] = 0.0000;
    f.kappa[0][0][x] = 0.0000;
    f.p[0][0][x] = 0.0000;

    for(i=0; i<TENSOR_CPTS; i++) {
      f.uij[i][0][0][x] = 0.0;
      f.udotij[i][0][0][x] = 0.0;
    }
  }

}

void free_fields(hydro_fields *f, hydro_params p) {

  free_field(p, f->phi);
  free_field(p, f->pifull);
  free_field(p, f->T);
  free_field(p, f->E);
  free_field(p, f->W);

  free_field(p, f->pi);

  free_field(p, f->phiold);

  free_field(p, f->kappa);
  free_field(p, f->p);

  /*
  free_field(p, f->V[0]);
  free_field(p, f->V[1]);
  free_field(p, f->V[2]);

  free(f->V);
  */
  free_vector(p, f->V);

  /*
  free_field(p, f->Z[0]);
  free_field(p, f->Z[1]);
  free_field(p, f->Z[2]);

  free(f->Z);
  */

  free_vector(p, f->Z);

  /*
  free_field(p, f->U[0]);
  free_field(p, f->U[1]);
  free_field(p, f->U[2]);

  free(f->U);
  */

  free_vector(p, f->U);

  /*
  free_field(p, f->F[0]);
  free_field(p, f->F[1]);
  free_field(p, f->F[2]);

  free(f->F);
  */

  free_vector(p, f->F);

  free_tensor(p, f->uij);

  free_tensor(p, f->udotij);

 
}

