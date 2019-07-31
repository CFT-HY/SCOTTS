/** @file silage_slice.c
 *
 * Module to write slices of silo data for visualisations.
 *
 * Slicing is done through a specified x coordinate, typically at
 * `x=p.siloslicecoord` where `p.siloslicecoord` is specified in the
 * input file.
 *
 * write_silo_slice_step() makes a silo file containing a slice
 * of the scalar field phi and the kinetic energy density. It calls make_kinetic()
 * to populate an array with the kinetic energy density and make_slice() to create
 * a slice out of some field, which is then populated on the rank 0 processor.
 *
 */
#include "hydro.h"

#ifdef SILO

#ifndef SCALAR
/** Populate an array with kinetic energy density.
 *
 */
void make_kinetic(hydro_fields f, hydro_params p, float ***temp) {
  float vol = p.dx*p.dx*p.dx;

  int x, y, z;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	temp[x][y][z] = f.kappa[x][y][z]*(f.E[x][y][z]/f.W[x][y][z])
          *(f.W[x][y][z]*f.W[x][y][z]-1.0)*vol;

      }
    }
  }

  halo_field(temp, p);
}

/** Make ordinary vorticity i.e vort = curl V
 */
void make_vort(hydro_fields f, hydro_params p, float ****temp){

  int x, y, z;
  float ****Vcell= make_vector(p);
  // Construct T velocity (centered at cell)

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	Vcell[0][x][y][z] = 0.5*(f.V[0][x][y][z] + f.V[0][x+1][y][z]);
	Vcell[1][x][y][z] = 0.5*(f.V[1][x][y][z] + f.V[1][x][y+1][z]);
	Vcell[2][x][y][z] = 0.5*(f.V[2][x][y][z] + f.V[2][x][y][(z+1)%p.Lz]);

      }
    }
  }

  halo_field(Vcell[0], p);
  halo_field(Vcell[1], p);
  halo_field(Vcell[2], p);

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
  	temp[0][x][y][z] = (Vcell[2][x][y+1][z] - Vcell[2][x][y-1][z]
			    - Vcell[1][x][y][(z+1)%p.Lz]
			    + Vcell[1][x][y][(z-1+p.Lz)%p.Lz])/(2*p.dx);
	temp[1][x][y][z] = (Vcell[0][x][y][(z+1)%p.Lz]
			    - Vcell[0][x][y][(z-1+p.Lz)%p.Lz]
			    - Vcell[2][x+1][y][z] + Vcell[2][x-1][y][z])/(2*p.dx);
	temp[2][x][y][z] = (Vcell[1][x+1][y][z] - Vcell[1][x-1][y][z]
			    - Vcell[0][x][y+1][z] + Vcell[0][x][y-1][z])/(2*p.dx);

      }
    }
  }

  free_vector(p, Vcell);


  halo_field(temp[0], p);
  halo_field(temp[1], p);
  halo_field(temp[2], p);


}

/** Make divergence of velocity field.
 */
void make_Vdiv(hydro_fields f, hydro_params p, float ***temp){
  int x, y, z;
    for(x = 1; x <= p.slicex; x++) {
      for(y = 1; y <= p.slicey; y++) {
	for(z = 0; z < p.Lz; z++) {
	  temp[x][y][z] = (f.V[0][x+1][y][z] - f.V[0][x][y][z]
			   + f.V[1][x][y+1][z] - f.V[1][x][y][z]
			   + f.V[2][x][y][z+1] - f.V[2][x][y][z])/p.dx;
	}
      }
    }

    halo_field(temp, p);
}

/** Make T vorticity i.e (Tvort = curl T U)
 */
void make_Tvort(hydro_fields f, hydro_params p, float ****temp){
  int x, y, z;

  float ****J = make_vector(p);

  // Construct temperature current (centered at cell)

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	J[0][x][y][z] = 0.5*(f.V[0][x][y][z] + f.V[0][x+1][y][z]
				)*f.T[x][y][z]*f.W[x][y][z];
	J[1][x][y][z] = 0.5*(f.V[1][x][y][z] + f.V[1][x][y+1][z]
				)*f.T[x][y][z]*f.W[x][y][z];
	J[2][x][y][z] = 0.5*(f.V[2][x][y][z] + f.V[2][x][y][(z+1)%p.Lz]
				)*f.T[x][y][z]*f.W[x][y][z];
      }
    }
  }

  halo_field(J[0], p);
  halo_field(J[1], p);
  halo_field(J[2], p);


  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	temp[0][x][y][z] = (J[2][x][y+1][z] - J[2][x][y-1][z]
			    - J[1][x][y][(z+1)%p.Lz]
			    + J[1][x][y][(z-1+p.Lz)%p.Lz])/(2*p.dx);
	temp[1][x][y][z] = (J[0][x][y][(z+1)%p.Lz]
			    - J[0][x][y][(z-1+p.Lz)%p.Lz]
			    - J[2][x+1][y][z] + J[2][x-1][y][z])/(2*p.dx);
	temp[2][x][y][z] = (J[1][x+1][y][z] - J[1][x-1][y][z]
			    - J[0][x][y+1][z] + J[0][x][y-1][z])/(2*p.dx);
      }
    }
  }
  free_vector(p, J);

  halo_field(temp[0], p);
  halo_field(temp[1], p);
  halo_field(temp[2], p);
}

/** Make divergence of temperature current \f$ J^{i} = TU \f$.
 */
void make_Jdiv(hydro_fields f, hydro_params p, float ***temp){
  int x, y, z;
  float ****J = make_vector(p);

  // Construct temperature current (centered at cell)

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	J[0][x][y][z] = 0.5*(f.V[0][x][y][z] + f.V[0][x+1][y][z]
				)*f.T[x][y][z]*f.W[x][y][z];
	J[1][x][y][z] = 0.5*(f.V[1][x][y][z] + f.V[1][x][y+1][z]
				)*f.T[x][y][z]*f.W[x][y][z];
	J[2][x][y][z] = 0.5*(f.V[2][x][y][z] + f.V[2][x][y][(z+1)%p.Lz]
				)*f.T[x][y][z]*f.W[x][y][z];
      }
    }
  }

  halo_field(J[0], p);
  halo_field(J[1], p);
  halo_field(J[2], p);

  for(x = 1; x <= p.slicex; x++) {
      for(y = 1; y <= p.slicey; y++) {
	for(z = 0; z < p.Lz; z++) {
	  temp[x][y][z] = (J[0][x+1][y][z] - J[0][x][y][z]
			   + J[1][x][y+1][z] - J[1][x][y][z]
			   + J[2][x][y][z+1] - J[2][x][y][z])/p.dx;
	}
      }
    }

    halo_field(temp, p);
    free_vector(p, J);
}

/** Populate an array with velocity magnitude.
 *
 */
void make_vel(hydro_fields f, hydro_params p, float ***temp) {
  float vol = p.dx*p.dx*p.dx;

  int x, y, z;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	temp[x][y][z] = sqrt(0.25*((f.V[0][x][y][z] + f.V[0][x+1][y][z])
				   *(f.V[0][x][y][z] + f.V[0][x+1][y][z])
				   + (f.V[1][x][y][z] + f.V[1][x][y+1][z])
				   *(f.V[1][x][y][z] + f.V[1][x][y+1][z])
				   + (f.V[2][x][y][z] + f.V[2][x][y][z+1])
				   *(f.V[2][x][y][z] + f.V[2][x][y][z+1])));

      }
    }
  }

  halo_field(temp, p);
}

/** Populate an array with Z magnitude.
 *
 */
void make_Z(hydro_fields f, hydro_params p, float ***temp) {
  float vol = p.dx*p.dx*p.dx;

  int x, y, z;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

   temp[x][y][z] = sqrt(f.Z[0][x][y][z]*f.Z[0][x][y][z]
                + f.Z[1][x][y][z]*f.Z[1][x][y][z]
                + f.Z[2][x][y][z]*f.Z[2][x][y][z]);

      }
    }
  }

  halo_field(temp, p);
}


#endif // !SCALAR

/** Populate an array with a scalar quantity in a slice through `temp`.
 *
 * This function takes a 3 dimensional array `temp`
 * which contains a quanitity defined over the simulation box
 * and populates the 1 dimensional array `slice` with a slice through
 * `x=xcoord`.
 *
 * Only processors containing the coordinate `x=xcoord`
 * will participate in this.
 *
 * First the 1 dimensional array `trim` is populated with the data on
 * the local slice. The size of trim is `p.slicey*p.Lz`, and it increases
 * first in `z` and then in `y`.
 *
 * Then the array slice is populated on the rank 0 processor
 * with `trim` from the other processors.
 * The size of `slice` is `p.Lz*p.Ly` and it increases first in `z`
 * and then in `y`.
 *
 *
 */

void make_slice(hydro_fields f, hydro_params p, int xcoord, float *slice, float ***temp)
{


  int x, y, z;

  x=xcoord;

  float *trim = (float *)malloc(p.slicey*p.Lz*sizeof(float));

  if(x/p.slicex == p.myposx){
    int localx=x%p.slicex+1;
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {
	trim[(y-1)*p.Lz + z] = temp[localx][y][z];
      }
    }
  }


  int ry;
  int ny = p.Ly/p.slicey;
  int nx = p.Lx/p.slicex;
  MPI_Status status;

  {

    // Check whether we are the target node for this slab
    if(!p.rank) {

      for(ry = 0; ry < ny; ry++) {

	if((x/p.slicex == p.myposx) && (ry == p.myposy)) {

	  memcpy(&slice[ry*p.slicey*p.Lz],
		 &trim[0],
		 p.slicey*p.Lz*sizeof(float));
	  continue;
	}


	MPI_Recv(&slice[ry*p.slicey*p.Lz],
		 p.slicey*p.Lz,
		 MPI_FLOAT,
		 ry*nx + x/p.slicex,
		 ry,
		 MPI_COMM_WORLD,
		 &status);
      }


    } else if((x >= p.shiftx) && (x < (p.shiftx + p.slicex))) {


      MPI_Send(&trim[0],
               p.slicey*p.Lz,
               MPI_FLOAT,
	       0,
               p.myposy,
	       MPI_COMM_WORLD);
    }

  }

  free(trim);

}


/** Write a silo file containing slices through the simulation box.
 *
 * This function saves in a silo file slices through
 * the coordinate `x=xcoord`.
 *
 * The `slice` arrays are of length `p.Lz*p.Ly`, and increase first
 * in `z` and then in `y`.
 *
 * Contributors:
 * - 2010-2017 David Weir
 * - 2018-     Daniel Cutting
 */
void write_silo_slice_step(hydro_fields f, hydro_params p, int step, int xcoord)
{


  int x, y, z, i;


  DBfile *dbfile = NULL;
  int *meshsize = NULL;
  float **mesh = NULL;
  DBoptlist *dboptlist = NULL;

  //// set up silage ////

  if(!p.rank) {

    /* Create a unique filename for the new Silo file */
    char filename[600];
    sprintf(filename, "%s/slice-%06d.silo", p.silodir, step);
    fprintf(stderr,"Writing SLICE to %s\n",filename);

    dbfile = DBCreate(filename, DB_CLOBBER, DB_LOCAL,
		      "time step", CPMODE);


    dboptlist = DBMakeOptlist(1);

    int col_major = DB_COLMAJOR;

    DBAddOption(dboptlist, DBOPT_MAJORORDER, &col_major);


    meshsize = (int *)malloc(2*sizeof(int));

    int sizey, sizez;

    sizey = p.Ly;
    sizez = p.Lz;

    meshsize[0] = sizey;
    meshsize[1] = sizez;

    mesh = (float **)malloc(2*sizeof(float *));

    mesh[0] = (float *)malloc(sizey*sizeof(float));
    mesh[1] = (float *)malloc(sizez*sizeof(float));

    for(x=0; x<sizey; x++) {
      mesh[0][x] = p.dx*((float)(x));
    }
    for(x=0; x<sizez; x++) {
      mesh[1][x] = p.dx*((float)(x));
    }

    DBPutQuadmesh(dbfile, "quadmesh", NULL, mesh, meshsize, 2,
		  DB_FLOAT, DB_COLLINEAR, NULL);

  }


  /////// prepare slice ///////

  float ***temp = make_field(p);
  float *slice = (float *)malloc(p.Ly*p.Lz*sizeof(float));

#ifndef SCALAR

  make_kinetic(f, p, temp);

  make_slice(f, p, xcoord, slice, temp);

  if(!p.rank){
    fprintf(stderr,"Slice coord x=%d\n",xcoord);
  }

  //// write kinetic slice ///


  if(!p.rank) {
    DBPutQuadvar1(dbfile, "kinetic", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }



  //// prepare fluid velocity ////

  make_vel(f, p, temp);

  make_slice(f, p, xcoord, slice, temp);


  //// write fluid velocity slice ///


  if(!p.rank) {
    DBPutQuadvar1(dbfile, "V", "quadmesh", slice, meshsize, 2,
         NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }


  //// prepare Z magnitude slice ////


  make_Z(f, p, temp);

  make_slice(f, p, xcoord, slice, temp);

  // write Z magnitude slice.

  if(!p.rank) {
    DBPutQuadvar1(dbfile, "Z", "quadmesh", slice, meshsize, 2,
         NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }


  //// prepare Vdiv slice

  make_Vdiv(f, p, temp);

  make_slice(f, p, xcoord, slice, temp);

  // write Vdiv slice.

  if(!p.rank) {
    DBPutQuadvar1(dbfile, "Vdiv", "quadmesh", slice, meshsize, 2,
         NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }

  //// prepare Jdiv slice

  make_Jdiv(f, p, temp);

  make_slice(f, p, xcoord, slice, temp);

  // write Vdiv slice.

  if(!p.rank) {
    DBPutQuadvar1(dbfile, "Jdiv", "quadmesh", slice, meshsize, 2,
         NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }

  free_field(p, temp);

  //// Write all vel components slice:

  // Vx slice
  make_slice(f, p, xcoord, slice, f.V[0]);
  if(!p.rank) {
    DBPutQuadvar1(dbfile, "Vx", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }

  // Vy slice
  make_slice(f, p, xcoord, slice, f.V[1]);
  if(!p.rank) {
    DBPutQuadvar1(dbfile, "Vy", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }

  // Vz slice
  make_slice(f, p, xcoord, slice, f.V[2]);
  if(!p.rank) {
    DBPutQuadvar1(dbfile, "Vz", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }

  //// Write all Z components slice:

  // Zx slice
  make_slice(f, p, xcoord, slice, f.Z[0]);
  if(!p.rank) {
    DBPutQuadvar1(dbfile, "Zx", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }

  // Zy slice
  make_slice(f, p, xcoord, slice, f.Z[1]);
  if(!p.rank) {
    DBPutQuadvar1(dbfile, "Zy", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }

  // Zz slice
  make_slice(f, p, xcoord, slice, f.Z[2]);
  if(!p.rank) {
    DBPutQuadvar1(dbfile, "Zz", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }


  float ****temp_vec = make_vector(p);

  //// Write all vorticity (vort) components slice:


  make_vort(f, p, temp_vec);

  // vortx slice
  make_slice(f, p, xcoord, slice, temp_vec[0]);
  if(!p.rank) {
    DBPutQuadvar1(dbfile, "vortx", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }

  // vorty slice
  make_slice(f, p, xcoord, slice, temp_vec[1]);
  if(!p.rank) {
    DBPutQuadvar1(dbfile, "vorty", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }

  // vortz slice
  make_slice(f, p, xcoord, slice, temp_vec[2]);
  if(!p.rank) {
    DBPutQuadvar1(dbfile, "vortz", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }


  //// Write all T vorticity (Tvort) components slice:


  make_Tvort(f, p, temp_vec);

  // Tvortx slice
  make_slice(f, p, xcoord, slice, temp_vec[0]);
  if(!p.rank) {
    DBPutQuadvar1(dbfile, "Tvortx", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }

  // Tvorty slice
  make_slice(f, p, xcoord, slice, temp_vec[1]);
  if(!p.rank) {
    DBPutQuadvar1(dbfile, "Tvorty", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }

  // Tvortz slice
  make_slice(f, p, xcoord, slice, temp_vec[2]);
  if(!p.rank) {
    DBPutQuadvar1(dbfile, "Tvortz", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }


  free_vector(p, temp_vec);





  //// make T slice ///

  make_slice(f, p, xcoord, slice, f.T);

  //// write T slice ///


  if(!p.rank) {
    DBPutQuadvar1(dbfile, "T", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }


  //// make E slice ///

  make_slice(f, p, xcoord, slice, f.E);

  //// write E slice ///


  if(!p.rank) {
    DBPutQuadvar1(dbfile, "E", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }



#endif // !SCALAR

  //// make phi slice ///

  make_slice(f, p, xcoord, slice, f.phi);

  //// write phi slice ///


  if(!p.rank) {
    DBPutQuadvar1(dbfile, "phi", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }


  free(slice);



  if(!p.rank) {
    DBClose(dbfile);

    free(mesh[0]);
    free(mesh[1]);
    free(mesh);
    free(meshsize);
  }



}



#endif // SILO
