/** @file fft_func.c
 *
 * Fourier transform functions for fields, vectors, tensors.
 *
 * Also includes initialisation and cleanup of fftw mpi requirements.
 *
 * When we perform fourier transforms we take fields that are defined
 * on the processor grid as pencils and then distribute them into
 * slabs. Once we perform the fourier transform fields exist in
 * momentum space as complex numbers defined on slabs.
 */
#include "hydro.h"


#ifdef FFT

/** Initialises the fftw fields and plans and mpi routines.
 *
 * Note this documentation should be extended.
 */
void fft_init(hydro_params p, fft_fields *fft_f){

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  
  //  int slab;

  MPI_Status status;

  int x, y, z;
  int i;

  
  printf0(p,"initialising fftw mpi\n");

  
  // First things first - initialise fftwf_mpi
  fftwf_mpi_init();


  printf0(p,"Creating map for nodes responsible for x coord in slab decomp\n");

  alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
					MPI_COMM_WORLD, &x_thickness, &x_start);

  
  // Now we initialise the map for which node is responsible for each
  // x coord in the slab decomposition.
  int *thicknesses = malloc(p.size*sizeof(int));
  int *starts = malloc(p.size*sizeof(int));
  fft_f->map = malloc(p.Lx*sizeof(int));
			    

  if(p.rank) {
    MPI_Send(&x_thickness, 1, MPI_INTEGER, 0, p.rank, MPI_COMM_WORLD);
    MPI_Send(&x_start, 1, MPI_INTEGER, 0, p.rank, MPI_COMM_WORLD);

  } else {
    thicknesses[0] = x_thickness;
    starts[0] = x_start;

    for(i=1; i<p.size; i++) {
      MPI_Recv(&thicknesses[i], 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, &status);
      MPI_Recv(&starts[i], 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, &status);

      //      fprintf(stderr,"rank %d, thickness %d, start %d\n", i, thicknesses[i], starts[i]);
    }


    for(x=0; x<p.Lx; x++) {
      fft_f->map[x] = -1;
      for(i=0; i<p.size; i++) {
	if(starts[i] <= x && (starts[i] + thicknesses[i]) > x) {
	  fft_f->map[x] = i;
	  //	  fprintf(stderr, "therefore node %d is responsible for x=%d\n", i, x);
	  break;
	}
      }
      if(fft_f->map[x] == -1) {
	fprintf(stderr,"cannot find a node responsible for x=%d!\n", x);
	die(-99);
      }

    }
  }
    
  printf0(p,"broadcasting our results\n");
  MPI_Bcast(fft_f->map, p.Lx, MPI_INTEGER, 0, MPI_COMM_WORLD);


  // Allocate in and out fields and plan for fourier transforms.
  

  fft_f->in = fftwf_alloc_complex(alloc_local);
  fft_f->out = fftwf_alloc_complex(alloc_local);

  fft_f->plan = fftwf_mpi_plan_dft_3d(p.Lx, p.Ly, p.Lz,
				     fft_f->in, fft_f->out, MPI_COMM_WORLD,
				     FFTW_FORWARD, FFTW_ESTIMATE);

  if (p.uetcstart>=0){
    fft_f->initial_Tij = (fftwf_complex **)malloc(6*sizeof(fftwf_complex *));
  
    for(i=0;i<TENSOR_CPTS;i++) {
      fft_f->initial_Tij[i] = fftwf_alloc_complex(alloc_local);
    }
  }

  
}

/** Cleanup the fft_fields fftwf fields, plan and mpi routines.
 */
void fft_finalise(hydro_params p, fft_fields *fft_f){

  int i;

  fftwf_destroy_plan(fft_f->plan);
  
  fftwf_free(fft_f->in);
  fftwf_free(fft_f->out);

  if (p.uetcstart>=0){

    for(i=0;i<TENSOR_CPTS;i++)
      fftwf_free(fft_f->initial_Tij[i]);

    free(fft_f->initial_Tij);
  }

  free(fft_f->map);
  fftwf_mpi_cleanup();
}


/** Populate `out` field in fft_fields with unnormalised fourier
 * transform of supplied real_field.
 *
 * Function which takes in real_field which is a float defined on the
 * lattice in real space (pencil decomposition), and populates the
 * `in` field in fft_fields with real_field now in slab
 * decomposition. Then executes the fftw_plan to populate the `out`
 * field in fft_fields with the `in` field transformed into momentum
 * space. Note that no normalisation is applied during this procedure.
 *
 */
void fft_scalar(hydro_params p, fft_fields fft_f, float ***real_field) {

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  
  //  int slab;

  MPI_Status status;

  int x, y, z;
  int i;

  //printf0(p, "FFT: alloc_local\n");
  
  alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
					MPI_COMM_WORLD, &x_thickness, &x_start);


  float *trim_field = (float *)malloc(p.slicex*p.slicey*p.Lz*sizeof(float));

  //printf0(p, "FFT: populating trim field\n");
  
  float checken = 0.0;
  for(x=1; x<=p.slicex; x++) {
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {

	trim_field[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z] = real_field[x][y][z];
      }
    }
  }

  float *slice = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));

  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;

  int ry;

  for(x=0; x<p.Lx; x++) {

    // Check whether we are the target node for this slab
    if(fft_f.map[x] == p.rank) {

      //fprintf(stderr, "%d: my slice, x=%d\n", p.rank, x);
      for(ry = 0; ry < ny; ry++) {

	if((x/p.slicex == p.myposx) && (ry == p.myposy)) {
	  memcpy(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 &trim_field[(x-p.shiftx)*p.slicey*p.Lz],
		 p.slicey*p.Lz*sizeof(float));
	  continue;
	}

	MPI_Recv(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 p.slicey*p.Lz,
		 MPI_FLOAT,
		 ry*nx + x/p.slicex,
		 x*ny + ry,
		 MPI_COMM_WORLD,
		 &status);

	
	//fprintf(stderr, "%d: got pencil with x=%d\n", p.rank,
	//	x);
	

      }

    } else if((x >= p.shiftx) && (x < (p.shiftx + p.slicex))) {

      
    //fprintf(stderr, "%d: my pencil, x=%d to dest %d/%d\n",
    //	      p.rank, x, x, fft_f.map[x]);
      

      MPI_Send(&trim_field[(x-p.shiftx)*p.slicey*p.Lz],
	       p.slicey*p.Lz,
	       MPI_FLOAT,
	       fft_f.map[x],
	       x*ny + p.myposy,
	       MPI_COMM_WORLD);
    }

    // then wait until that slice of x is sorted
    
  }
  // Alloc slab array and do communication to get it into place

  //printf0(p, "FFT: Done slicing\n");



  for(x=0;x<x_thickness;x++) {
    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {
	
	fft_f.in[x*p.Ly*p.Lz + y*p.Lz + z][0] = slice[x*p.Ly*p.Lz + y*p.Lz + z];
	fft_f.in[x*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;
      }
    }
  }



  
  fftwf_mpi_execute_dft(fft_f.plan, fft_f.in, fft_f.out);

  //printf0(p, "Basic FFT done\n");

  
  // Perform cleanup on declared fields.
  
  free(trim_field);
  
  free(slice);

}

/** Populate outcpts with normalised fourier transform of
 * `vector_field`. Uses `fft_scalar` on each component of the vector.
 *
 * Note that this will take `vector_field` which is decomposed in
 * pencils and populate outcpts with the fourier transform decomposed
 * in slabs.
 *
 */
void fft_vector(hydro_params p, fft_fields fft_f, float ****vector_field,
		fftwf_complex **outcpts) {

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  int x, y, z;
  int i;

  alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
					MPI_COMM_WORLD, &x_thickness, &x_start);


  
  float fft_norm = (1.0/(((float)p.Lx)*((float)p.Ly)*((float)p.Lz)));

  
  for(i = 0; i < 3; i++) {
    fft_scalar(p, fft_f, vector_field[i]);

    for(x=0;x<x_thickness;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  outcpts[i][x*p.Ly*p.Lz + y*p.Lz + z][0] 
	    = fft_norm*fft_f.out[x*p.Ly*p.Lz + y*p.Lz + z][0];
	  outcpts[i][x*p.Ly*p.Lz + y*p.Lz + z][1]
	    = fft_norm*fft_f.out[x*p.Ly*p.Lz + y*p.Lz + z][1];

	  //	  total += outcpts[i][x*p.Ly*p.Lz + y*p.Lz + z][0];
	}
      }
    }

    //printf0(p, "FFT vector: Completed fft for cpt %d\n", i);

  }
}

/** Populate outcpts with normalised fourier transform of
 * `tensor_field`. Uses `fft_scalar` on each component of the tensor.
 *
 * Takes an additional normalisation `norm`, so we can maintain
 * consistency for the gravitational wave power spectrum which is
 * normalised by an additional \f$1/(\sqrt{32 \pi})\f$.
 *
 * Note that this will take `tensor_field` which is decomposed in
 * pencils and populate outcpts with the fourier transform decomposed
 * in slabs.
 */
void fft_tensor(hydro_params p, fft_fields fft_f, float ****tensor_field,
		fftwf_complex **outcpts, float norm) {

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  
  int x, y, z;
  int i;

  alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
					MPI_COMM_WORLD, &x_thickness, &x_start);



  
  float fft_norm = (1.0/(((float)p.Lx)*((float)p.Ly)*((float)p.Lz)));
  float tot_norm = fft_norm*norm;
  for(i = 0; i < TENSOR_CPTS; i++) {
    fft_scalar(p, fft_f, tensor_field[i]);

    for(x=0;x<x_thickness;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  outcpts[i][x*p.Ly*p.Lz + y*p.Lz + z][0] 
	    = tot_norm*fft_f.out[x*p.Ly*p.Lz + y*p.Lz + z][0];
	  outcpts[i][x*p.Ly*p.Lz + y*p.Lz + z][1]
	    = tot_norm*fft_f.out[x*p.Ly*p.Lz + y*p.Lz + z][1];

	  //	  total += outcpts[i][x*p.Ly*p.Lz + y*p.Lz + z][0];
	}
      }
    }

    //printf0(p, "FFT tensor: Completed fft for cpt %d\n", i);

  }
}


/** Populate outcpts with momentum space fourier transform of
 * temperature current \f$ J^{i}=T U^{i} \f$.
 *
 */
void fft_J(hydro_fields f, hydro_params p, fft_fields fft_f,
	   fftwf_complex **outcpts){
#ifndef SCALAR
  int x, y, z, i;
  float ****J = make_vector(p);
  float Tnode;

  // Construct temperature current.
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	Tnode = 0.125*(f.T[x][y][z] + f.T[x-1][y][z]
		       + f.T[x][y-1][z] + f.T[x][y][(z-1+p.Lz)%p.Lz]
		       + f.T[x-1][y-1][z] + f.T[x-1][y][(z-1+p.Lz)%p.Lz]
		       + f.T[x][y-1][(z-1+p.Lz)%p.Lz]
		       + f.T[x-1][y-1][(z-1+p.Lz)%p.Lz]
		       );
	for(i =0; i < 3; i++){
	  J[i][x][y][z] = f.U[i][x][y][z]*Tnode;
	}
      }
    }
  }

  halo_field(J[0], p);
  halo_field(J[1], p);
  halo_field(J[2], p);

  fft_vector(p, fft_f, J, outcpts);
  
  free_vector(p, J);
    
    
#endif //!SCALAR
}


/** Populate outcpts with momentum space fourier transform of
 * \f$ X^{i}=\sqrt{w}U^{i} \f$.
 *
 */
 void fft_X(hydro_fields f, hydro_params p, fft_fields fft_f,
	    fftwf_complex **outcpts){
#ifndef SCALAR
  int x, y, z, i;
  float ****X = make_vector(p);
  float sqrt_enth_node;

  // Construct X vector.
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	sqrt_enth_node = pow((f.kappa[x][y][z]*f.E[x][y][z]
			      /f.W[x][y][z]
			      + f.kappa[x-1][y][z]*f.E[x-1][y][z]
			      /f.W[x-1][y][z]
			      + f.kappa[x][y-1][z]*f.E[x][y-1][z]
			      /f.W[x][y-1][z]
			      + f.kappa[x][y][(z-1+p.Lz)%p.Lz]
			      *f.E[x][y][(z-1+p.Lz)%p.Lz]
			      /f.W[x][y][(z-1+p.Lz)%p.Lz]
			      + f.kappa[x-1][y-1][z]*f.E[x-1][y-1][z]
			      /f.W[x-1][y-1][z]
			      + f.kappa[x-1][y][(z-1+p.Lz)%p.Lz]
			      *f.E[x-1][y][(z-1+p.Lz)%p.Lz]
			      /f.W[x-1][y][(z-1+p.Lz)%p.Lz]
			      + f.kappa[x][y-1][(z-1+p.Lz)%p.Lz]
			      *f.E[x][y-1][(z-1+p.Lz)%p.Lz]
			      /f.W[x][y-1][(z-1+p.Lz)%p.Lz]
			      + f.kappa[x-1][y-1][(z-1+p.Lz)%p.Lz]
			      *f.E[x-1][y-1][(z-1+p.Lz)%p.Lz]
			      /f.W[x-1][y-1][(z-1+p.Lz)%p.Lz]
			      )*0.125, 0.5);
	for(i =0; i < 3; i++){
	  X[i][x][y][z] = sqrt_enth_node*f.U[i][x][y][z];
	}
      }
    }
  }

  halo_field(X[0], p);
  halo_field(X[1], p);
  halo_field(X[2], p);
  
  fft_vector(p, fft_f, X, outcpts);
  
  free_vector(p, X);
    
    
#endif //!SCALAR
}

 
/** Create internal energy field e = E/W and then perform fft of this field.
 *
 * Not a fan of this implementation, must be something neater/other 
 * place for this to go.
 */
void fft_e(hydro_fields f, hydro_params p, fft_fields fft_f){
#ifndef SCALAR

  int x, y, z;
  float ***e = make_field(p);
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	e[x][y][z] = f.E[x][y][z]/f.W[x][y][z];
      }
    }
  }
  halo_field(e, p);

  fft_scalar(p, fft_f, e);

  free_field(p, e);
#endif //!SCALAR
}

  
#endif // FFT
