/* silage.c
 *
 * Write data for the purposes of visualisation. Note that checkpointing
 * (see checkpoint.c) also uses silo as a wrapper around hdf5 storage.
 */
#include "hydro.h"

#ifdef SILO

/* void silo_init(hydro_params p)
 *
 * set up silo data repository directory
 */
void silo_init(hydro_params p)
{

  //  if(p.silostep > 0) {
  //    fprintf(stderr,"p.silostep is %d\n", p.silostep);
  
  //  char silodir[100];
  //  sprintf(silodir,"silage-%d",(int)getpid());
  if(mkdir(p.silodir,07777) && !p.rank)
    perror(p.silodir);
  
  if(!p.rank)
    fprintf(stderr,"created directory %s, where the silo will go\n",
	    p.silodir);

  //  DBSetCompression("METHOD=GZIP");

  //  }
}


/* void write_silo_step(hydro_fields f, hydro_params p, int step)
 *
 * write fields to a silo file for use with Visit
 */
void write_silo_step(hydro_fields f, hydro_params p, int step)
{

  int x, y, z, i;

  DBfile *dbfile = NULL;

  //  char silodir[100];
  //  sprintf(silodir,"silage-%d",(int)getpid());

  /* Create a unique filename for the new Silo file */
  char filename[600];
  sprintf(filename, "%s/output-%d-%06d.silo", p.silodir, p.rank, step);  
  fprintf(stderr,"Writing step to %s\n",filename);  
  dbfile = DBCreate(filename, DB_CLOBBER, DB_LOCAL,
                    "time step", CPMODE);

  DBoptlist *dboptlist = NULL;

  dboptlist = DBMakeOptlist(1);

  int col_major = DB_COLMAJOR;

  DBAddOption(dboptlist, DBOPT_MAJORORDER, &col_major);

  int *meshsize = (int *)malloc(3*sizeof(int));

  int sizex, sizey, sizez;

  sizex = p.slicex+2;
  sizey = p.slicey+2;
  sizez = p.Lz;




  meshsize[0] = sizex;
  meshsize[1] = sizey;
  meshsize[2] = sizez;

  double **mesh = (double **)malloc(3*sizeof(double *));

  
  mesh[0] = (double *)malloc(sizex*sizeof(double));
  mesh[1] = (double *)malloc(sizey*sizeof(double));
  mesh[2] = (double *)malloc(sizez*sizeof(double));


  for(x=0; x<sizex; x++) {
    mesh[0][x] = p.dx*((double)(x + p.shiftx - 1));
  }
  for(x=0; x<sizey; x++) {
    mesh[1][x] = p.dx*((double)(x + p.shifty - 1));
  }
  for(x=0; x<sizez; x++) {
      mesh[2][x] = p.dx*((double)x);
  }

  


  DBPutQuadmesh(dbfile, "quadmesh", NULL, mesh, meshsize, 3,
		DB_DOUBLE, DB_COLLINEAR, NULL);

  DBPutQuadvar1(dbfile, "phi", "quadmesh", f.phi[0][0], meshsize, 3,
		NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);

#ifndef SCALAR
  DBPutQuadvar1(dbfile, "E", "quadmesh", f.E[0][0], meshsize, 3,
		NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);
#endif
  
  double ***temp = make_field(p);

  double vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	
	temp[x][y][z] = 0.0;
        
	temp[x][y][z] += 0.5*((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)
	    *((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)*vol;
        
	temp[x][y][z] += 0.5*((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)
	    *((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)*vol;
        
	temp[x][y][z] += 0.5*((f.phi[x][y][(z+1)%p.Lz] - f.phi[x][y][z])/p.dx)
	    *((f.phi[x][y][(z+1)%p.Lz] - f.phi[x][y][z])/p.dx)*vol;
        
      }
    }
  }

  halo_field(temp, p);
  
  DBPutQuadvar1(dbfile, "gradphi", "quadmesh", temp[0][0], meshsize, 3,
		NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);


#ifndef SCALAR

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	
	temp[x][y][z] = f.kappa[x][y][z]*(f.E[x][y][z]/f.W[x][y][z])
          *(f.W[x][y][z]*f.W[x][y][z]-1.0)*vol;
        
      }
    }
  }
  
  halo_field(temp, p);


  DBPutQuadvar1(dbfile, "kinetic", "quadmesh", temp[0][0], meshsize, 3,
		NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);
#endif


  free_field(p, temp);


  // Don't store vector quantities because they are huuuuge on 1024^3...  

  /*
  char *ux_name = "Ux";
  char *uy_name = "Uy";
  char *uz_name = "Uz";

  char *u_names[3];
  u_names[0] = ux_name;
  u_names[1] = uy_name;
  u_names[2] = uz_name;

  DBPutQuadvar(dbfile, "U", "quadmesh", 3, u_names, f.U, meshsize, 3,
	       NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);

  */

  /*
  fprintf(stderr, "Writing V\n");

  char *vx_name = "Vx";
  char *vy_name = "Vy";
  char *vz_name = "Vz";

  char *v_names[3];
  v_names[0] = vx_name;
  v_names[1] = vy_name;
  v_names[2] = vz_name;


  double **Vtemp = (double **)malloc(3*sizeof(double *));
  Vtemp[0] = f.V[0][0][0];
  Vtemp[1] = f.V[1][0][0];
  Vtemp[2] = f.V[2][0][0];


  DBPutQuadvar(dbfile, "V", "quadmesh", 3, v_names, Vtemp, meshsize, 3,
  	       NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);

  */



  /*
  char *zx_name = "Zx";
  char *zy_name = "Zy";
  char *zz_name = "Zz";

  char *z_names[3];
  z_names[0] = zx_name;
  z_names[1] = zy_name;
  z_names[2] = zz_name;

  DBPutQuadvar(dbfile, "Z", "quadmesh", 3, z_names, f.Z, meshsize, 3,
	       NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);

  */


  DBClose(dbfile);

  // Only needed if we are recording the vector quantities.
  //  free(Vtemp);

  for(i=0;i<3;i++) {
    free(mesh[i]);
  }
  free(meshsize);
  free(mesh);








  /* Write the root, which specifies where the other files
   * may be found and what they contain.
   */
  if(!p.rank) {
    sprintf(filename, "%s/root-%06d.silo", p.silodir, step);  
    fprintf(stderr,"Writing root to %s\n",filename);  
    dbfile = DBCreate(filename, DB_CLOBBER, DB_LOCAL,
		      "time step", CPMODE);
    

    char **names = (char **)malloc(p.size*sizeof(char *));
    int *types = (int *)malloc(p.size*sizeof(int));

    for(i=0; i<p.size;i++) {
      names[i] = (char *)malloc(600*sizeof(char));
    }


    for(i=0;i<p.size;i++) {
      types[i] = DB_QUAD_RECT;
      sprintf(names[i], "output-%d-%06d.silo:quadmesh", i, step);
    }

    DBPutMultimesh(dbfile, "quadmesh", p.size, names, types, NULL);

#ifndef SCALAR
    for(i=0;i<p.size;i++) {
      types[i] = DB_QUADVAR;
      sprintf(names[i], "output-%d-%06d.silo:E", i, step);
    }

    DBPutMultivar(dbfile, "E", p.size, names, types, NULL);
#endif

    for(i=0;i<p.size;i++) {
      types[i] = DB_QUADVAR;
      sprintf(names[i], "output-%d-%06d.silo:phi", i, step);
    }

    DBPutMultivar(dbfile, "phi", p.size, names, types, NULL);


    for(i=0;i<p.size;i++) {
      types[i] = DB_QUADVAR;
      sprintf(names[i], "output-%d-%06d.silo:gradphi", i, step);
    }

    DBPutMultivar(dbfile, "gradphi", p.size, names, types, NULL);


#ifndef SCALAR
    for(i=0;i<p.size;i++) {
      types[i] = DB_QUADVAR;
      sprintf(names[i], "output-%d-%06d.silo:kinetic", i, step);
    }

    DBPutMultivar(dbfile, "kinetic", p.size, names, types, NULL);
#endif


    for(i=0;i<p.size;i++) {
      free(names[i]);
    }

    DBClose(dbfile);

    free(names);
    free(types);
	
  }

}  
  


#endif // SILO

