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


  for(i=0;i<3;i++)
    meshsize[i] = p.L;

  double **mesh = (double **)malloc(p.N*sizeof(double *));

  for(i=0;i<3;i++)
    mesh[i] = (double *)malloc(p.N*sizeof(double));

  for(x=0; x<p.L; x++) {
    for(i=0; i<3; i++) {
      mesh[i][x] = p.dx*((double)x);
    }
  }

  


  DBPutQuadmesh(dbfile, "quadmesh", NULL, mesh, meshsize, 3,
		DB_DOUBLE, DB_COLLINEAR, NULL);

  DBPutQuadvar1(dbfile, "E", "quadmesh", f.E, meshsize, 3,
		NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);

  DBPutQuadvar1(dbfile, "Zx", "quadmesh", f.Zx, meshsize, 3,
		NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);

  DBPutQuadvar1(dbfile, "Zy", "quadmesh", f.Zy, meshsize, 3,
		NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);

  DBPutQuadvar1(dbfile, "Zz", "quadmesh", f.Zz, meshsize, 3,
		NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);


  DBClose(dbfile);


  for(i=0;i<3;i++) {
    free(mesh[i]);
  }

  free(meshsize);
  free(mesh);

}  
  


#endif // SILO

