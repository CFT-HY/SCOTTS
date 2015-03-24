/* silage.c
 *
 * Write data for the purposes of visualisation. Note that checkpointing
 * (see checkpoint.c) also uses silo as a wrapper around hdf5 storage.
 */
#include "hydro.h"

#ifdef SILO


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

void make_slice(hydro_fields f, hydro_params p, float *slice, float ***temp)
{
  int x, y, z;

  float *trim = (float *)malloc(p.slicex*p.slicey*p.Lz*sizeof(float));
  
  for(x=1; x<=p.slicex; x++) {
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {
	trim[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z] = temp[x][y][z];
      }
    }
  }
  
  x = 0;

  int ry;
  int ny = p.Ly/p.slicey;
  int nx = p.Lx/p.slicex;
  MPI_Status status;

  {

    // Check whether we are the target node for this slab
    if(!p.rank) {

      for(ry = 0; ry < ny; ry++) {

	if((x/p.slicex == p.myposx) && (ry == p.myposy)) {

	  memcpy(&slice[(x-0)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 &trim[(x-p.shiftx)*p.slicey*p.Lz],
		 p.slicey*p.Lz*sizeof(float));
	  continue;
	}


	MPI_Recv(&slice[(x-0)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 p.slicey*p.Lz,
		 MPI_FLOAT,
		 ry*nx + x/p.slicex,
		 x*ny + ry,
		 MPI_COMM_WORLD,
		 &status);
      }


    } else if((x >= p.shiftx) && (x < (p.shiftx + p.slicex))) {


      MPI_Send(&trim[(x-p.shiftx)*p.slicey*p.Lz],
               p.slicey*p.Lz,
               MPI_FLOAT,
	       x/1,
               x*ny + p.myposy,
	       MPI_COMM_WORLD);
    }

  }

  free(trim);

  }

/* void write_silo_step(hydro_fields f, hydro_params p, int step)
 *
 * write fields to a silo file for use with Visit
 */
void write_silo_slice_step(hydro_fields f, hydro_params p, int step)
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

  make_kinetic(f, p, temp);

  make_slice(f, p, slice, temp);

  //// write kinetic slice ///


  if(!p.rank) {
    DBPutQuadvar1(dbfile, "kinetic", "quadmesh", slice, meshsize, 2,
		  NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
  }


  free_field(p, temp);

  //// make phi slice ///

  make_slice(f, p, slice, f.phi);

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

