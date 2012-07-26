#include "hydro.h"

#ifdef SILO

/* silo_init
 *
 * set up silo data repository directory
 */
void silo_init(hydro_params p)
{

  //  if(p.silostep > 0) {
  //    fprintf(stderr,"p.silostep is %d\n", p.silostep);
  
  char silodir[100];
  sprintf(silodir,"silage-%d",(int)getpid());
  if(mkdir(silodir,07777))
    perror(silodir);
  
  fprintf(stderr,"created directory %s, where the silo will go\n",
	  silodir);
  //  }
}


/* write_silo_step
 *
 * write fields to a silo file for use with Visit
 */
void write_silo_step(hydro_fields f, hydro_params p, int step)
{

  int x, i;

  DBfile *dbfile = NULL;

  char silodir[100];
  sprintf(silodir,"silage-%d",(int)getpid());

  /* Create a unique filename for the new Silo file */
  char filename[100];
  sprintf(filename, "%s/output%06d.silo", silodir, step);  
  fprintf(stderr,"Writing step to %s\n",filename);  
  dbfile = DBCreate(filename, DB_CLOBBER, DB_LOCAL,
                    "time step", DB_PDB);
  

  int *meshsize = (int *)malloc(3*sizeof(int));


  meshsize[0] = p.Lx;
  meshsize[1] = p.Ly;
  meshsize[2] = p.Lz;

  double **mesh = (double **)malloc(p.N*sizeof(double *));

  
  mesh[0] = (double *)malloc(p.Lx*sizeof(double));
  mesh[1] = (double *)malloc(p.Ly*sizeof(double));
  mesh[2] = (double *)malloc(p.Lz*sizeof(double));

  for(x=0; x<p.Lx; x++) {
      mesh[0][x] = p.dx*((double)x);
  }
  for(x=0; x<p.Ly; x++) {
      mesh[1][x] = p.dx*((double)x);
  }
  for(x=0; x<p.Lz; x++) {
      mesh[2][x] = p.dx*((double)x);
  }

  


  DBPutQuadmesh(dbfile, "quadmesh", NULL, mesh, meshsize, 3,
		DB_DOUBLE, DB_COLLINEAR, NULL);

  DBPutQuadvar1(dbfile, "E", "quadmesh", f.E, meshsize, 3,
		NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);

  DBPutQuadvar1(dbfile, "phi", "quadmesh", f.phi, meshsize, 3,
		NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);


  char *ux_name = "Ux";
  char *uy_name = "Uy";
  char *uz_name = "Uz";

  char *u_names[3];
  u_names[0] = ux_name;
  u_names[1] = uy_name;
  u_names[2] = uz_name;

  DBPutQuadvar(dbfile, "U", "quadmesh", 3, u_names, f.U, meshsize, 3,
	       NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);

  char *vx_name = "Vx";
  char *vy_name = "Vy";
  char *vz_name = "Vz";

  char *v_names[3];
  v_names[0] = vx_name;
  v_names[1] = vy_name;
  v_names[2] = vz_name;

  DBPutQuadvar(dbfile, "V", "quadmesh", 3, v_names, f.V, meshsize, 3,
	       NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);


  char *zx_name = "Zx";
  char *zy_name = "Zy";
  char *zz_name = "Zz";

  char *z_names[3];
  z_names[0] = zx_name;
  z_names[1] = zy_name;
  z_names[2] = zz_name;

  DBPutQuadvar(dbfile, "Z", "quadmesh", 3, z_names, f.Z, meshsize, 3,
	       NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);


  DBClose(dbfile);


  for(i=0;i<3;i++) {
    free(mesh[i]);
  }
  free(meshsize);
  free(mesh);

}  
  


#endif // SILO

