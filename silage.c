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

  int row_major = DB_ROWMAJOR;

  DBAddOption(dboptlist, DBOPT_MAJORORDER, &row_major);

  int hi_offset[3];
  int lo_offset[3];
  
  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;
  int this_posx = p.rank % nx;
  int this_posy = (p.rank - this_posx)/nx;
  // z-dir
  lo_offset[0] = 0;
  hi_offset[0] = 0;
  if(this_posx == 0)
    lo_offset[2] = 1;
  else
    lo_offset[2] = 0;
  hi_offset[2] = 1;
  if(this_posy == 0)
    lo_offset[1] = 1;
  else
    lo_offset[1] = 0;
  hi_offset[1] = 1;
  DBAddOption(dboptlist, DBOPT_HI_OFFSET, (void *)hi_offset);
  DBAddOption(dboptlist, DBOPT_LO_OFFSET, (void *)lo_offset);

  int *meshsize = (int *)malloc(3*sizeof(int));

  int sizex, sizey, sizez;
  
  sizex = p.slicex+2;
  sizey = p.slicey+2;
  sizez = p.Lz;




  meshsize[0] = sizez;
  meshsize[1] = sizey;
  meshsize[2] = sizex;

  float **mesh = (float **)malloc(3*sizeof(float *));

  
  mesh[0] = (float *)malloc(sizez*sizeof(float));
  mesh[1] = (float *)malloc(sizey*sizeof(float));
  mesh[2] = (float *)malloc(sizex*sizeof(float));


  for(x=0; x<sizex; x++) {
    mesh[2][x] = p.dx*((float)(x + p.shiftx - 1));
  }
  for(x=0; x<sizey; x++) {
    mesh[1][x] = p.dx*((float)(x + p.shifty - 1));
  }
  for(x=0; x<sizez; x++) {
      mesh[0][x] = p.dx*((float)x);
  }

  


  DBPutQuadmesh(dbfile, "quadmesh", NULL, mesh, meshsize, 3,
		DB_FLOAT, DB_COLLINEAR, NULL);

  DBPutQuadvar1(dbfile, "phi", "quadmesh", f.phi[0][0], meshsize, 3,
		NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);

#ifndef SCALAR
  DBPutQuadvar1(dbfile, "E", "quadmesh", f.E[0][0], meshsize, 3,
		NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);

  
  DBPutQuadvar1(dbfile, "T", "quadmesh", f.T[0][0], meshsize, 3,
		NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);
#endif
  
  float ***temp = make_field(p);

  float vol = p.dx*p.dx*p.dx;

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
		NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);


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
		NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);


#endif


  free_field(p, temp);


  // Don't store vector quantities because they are huuuuge on 1024^3...  

  fprintf(stderr, "Writing vectors\n");
  
  char *ux_name = "Ux";
  char *uy_name = "Uy";
  char *uz_name = "Uz";

  char *u_names[3];
  u_names[2] = ux_name;
  u_names[1] = uy_name;
  u_names[0] = uz_name;

  float **vec_temp = (float **)malloc(3*sizeof(float *));
  
  vec_temp[0] = f.U[2][0][0];
  vec_temp[1] = f.U[1][0][0];
  vec_temp[2] = f.U[0][0][0];

  DBPutQuadvar(dbfile, "U", "quadmesh", 3, u_names, vec_temp, meshsize, 3,
	       NULL, 0, DB_FLOAT, DB_NODECENT, NULL);

  

  
  

  char *vx_name = "Vx";
  char *vy_name = "Vy";
  char *vz_name = "Vz";

  char *v_names[3];
  v_names[2] = vx_name;
  v_names[1] = vy_name;
  v_names[0] = vz_name;


  vec_temp[0] = f.V[2][0][0];
  vec_temp[1] = f.V[1][0][0];
  vec_temp[2] = f.V[0][0][0];


  DBPutQuadvar(dbfile, "V", "quadmesh", 3, v_names, vec_temp, meshsize, 3,
  	       NULL, 0, DB_FLOAT, DB_NODECENT, dboptlist);


  
  char *zx_name = "Zx";
  char *zy_name = "Zy";
  char *zz_name = "Zz";

  char *z_names[3];
  z_names[2] = zx_name;
  z_names[1] = zy_name;
  z_names[0] = zz_name;

  vec_temp[2] = f.Z[2][0][0];
  vec_temp[1] = f.Z[1][0][0];
  vec_temp[0] = f.Z[0][0][0];


  DBPutQuadvar(dbfile, "Z", "quadmesh", 3, z_names, vec_temp, meshsize, 3,
	       NULL, 0, DB_FLOAT, DB_NODECENT, NULL);

  

  DBClose(dbfile);

  // Only needed if we are recording the vector quantities.
  free(vec_temp);
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
    
    DBoptlist *mesh_optlist = NULL;
    DBoptlist *trphisq_optlist = NULL;

    mesh_optlist = DBMakeOptlist(2);
    trphisq_optlist = DBMakeOptlist(3);

    const int two = 2;
    const int six = 6;

    double *spatial_extents;


    int node_posx, node_posy, node_shiftx, node_shifty;

    spatial_extents = (double *) malloc(6*p.size*sizeof(double));

    for(i=0; i<p.size; i++) {


      node_posx = i % nx;
      node_posy = (i - node_posx)/nx;

      node_shiftx = (p.slicex)*(node_posx);
      node_shifty = (p.slicey)*(node_posy);




      spatial_extents[i*6 + 0] = 0.0;

      if(node_posy == 0)
	spatial_extents[i*6 + 1] = 0; // y min
      else
	spatial_extents[i*6 + 1] = node_shifty - 1; // y min

      if(node_posx == 0)
	      spatial_extents[i*6 + 2] = 0; // x min
      else
	spatial_extents[i*6 + 2] = node_shiftx - 1; // x min

      spatial_extents[i*6 + 3] = p.Lz - 1;

      spatial_extents[i*6 + 4] = node_shifty + p.slicey - 1; // y max
      spatial_extents[i*6 + 5] = node_shiftx + p.slicex - 1; // x max



    }

    DBAddOption(mesh_optlist, DBOPT_EXTENTS_SIZE, (void *)&six);

    DBAddOption(mesh_optlist, DBOPT_EXTENTS, (void *)spatial_extents);


    for(i=0;i<p.size;i++) {
      types[i] = DB_QUAD_RECT;
      sprintf(names[i], "output-%d-%06d.silo:quadmesh", i, step);
    }

    DBPutMultimesh(dbfile, "quadmesh", p.size, names, types, mesh_optlist);

#ifndef SCALAR
    for(i=0;i<p.size;i++) {
      types[i] = DB_QUADVAR;
      sprintf(names[i], "output-%d-%06d.silo:E", i, step);
    }

    DBPutMultivar(dbfile, "E", p.size, names, types, NULL);

    for(i=0;i<p.size;i++) {
      types[i] = DB_QUADVAR;
      sprintf(names[i], "output-%d-%06d.silo:T", i, step);
    }

    DBPutMultivar(dbfile, "T", p.size, names, types, NULL);
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

    for(i=0;i<p.size;i++) {
      types[i] = DB_QUADVAR;
      sprintf(names[i], "output-%d-%06d.silo:V", i, step);
    }

    DBPutMultivar(dbfile, "V", p.size, names, types, NULL);
   
    for(i=0;i<p.size;i++) {
      types[i] = DB_QUADVAR;
      sprintf(names[i], "output-%d-%06d.silo:U", i, step);
    }

    DBPutMultivar(dbfile, "U", p.size, names, types, NULL);

    for(i=0;i<p.size;i++) {
      types[i] = DB_QUADVAR;
      sprintf(names[i], "output-%d-%06d.silo:Z", i, step);
    }

    DBPutMultivar(dbfile, "Z", p.size, names, types, NULL);
    
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

