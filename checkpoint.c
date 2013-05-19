/* checkpoint.c
 * 
 * Checkpoint and restore configuration.
 */
#include "hydro.h"


/* int load_checkpoint(hydro_fields f, hydro_params p)
 *
 * Restore checkpointed data to hydro_fields struct. Returns
 * the iteration step to which the state was restored.
 */
int load_checkpoint(hydro_fields f, hydro_params p)
{

  char filename[600];

  sprintf(filename, "%s/checkpoint-%06d", p.checkpointdir, p.rank);

  DBfile *dbfile = DBOpen(filename, CPMODE, DB_READ);

  DBquadvar *qvptr;
  
  int sizex, sizey, sizez, step;

  DBReadVar(dbfile, "sizex", &sizex);
  DBReadVar(dbfile, "sizey", &sizey);
  DBReadVar(dbfile, "sizez", &sizez);
  DBReadVar(dbfile, "step", &step);


  qvptr = DBGetQuadvar(dbfile, "phi");
  memcpy(f.phi[0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  free(qvptr);

  qvptr = DBGetQuadvar(dbfile, "pifull");
  memcpy(f.pifull[0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  free(qvptr);

  qvptr = DBGetQuadvar(dbfile, "pi");
  memcpy(f.pi[0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  free(qvptr);

  qvptr = DBGetQuadvar(dbfile, "phiold");
  memcpy(f.phiold[0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  free(qvptr);



#ifndef SCALAR
  qvptr = DBGetQuadvar(dbfile, "T");
  memcpy(f.T[0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  free(qvptr);

  qvptr = DBGetQuadvar(dbfile, "E");
  memcpy(f.E[0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  free(qvptr);

  qvptr = DBGetQuadvar(dbfile, "W");
  memcpy(f.W[0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  free(qvptr);

  qvptr = DBGetQuadvar(dbfile, "kappa");
  memcpy(f.kappa[0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  free(qvptr);

  qvptr = DBGetQuadvar(dbfile, "p");
  memcpy(f.p[0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  free(qvptr);



  qvptr = DBGetQuadvar(dbfile, "V");
  memcpy(f.V[0][0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  memcpy(f.V[1][0][0], qvptr->vals[1], sizex*sizey*sizez*sizeof(double));
  memcpy(f.V[2][0][0], qvptr->vals[2], sizex*sizey*sizez*sizeof(double));
  free(qvptr);

  qvptr = DBGetQuadvar(dbfile, "Z");
  memcpy(f.Z[0][0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  memcpy(f.Z[1][0][0], qvptr->vals[1], sizex*sizey*sizez*sizeof(double));
  memcpy(f.Z[2][0][0], qvptr->vals[2], sizex*sizey*sizez*sizeof(double));
  free(qvptr);

  qvptr = DBGetQuadvar(dbfile, "U");
  memcpy(f.U[0][0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  memcpy(f.U[1][0][0], qvptr->vals[1], sizex*sizey*sizez*sizeof(double));
  memcpy(f.U[2][0][0], qvptr->vals[2], sizex*sizey*sizez*sizeof(double));
  free(qvptr);

  qvptr = DBGetQuadvar(dbfile, "F");
  memcpy(f.F[0][0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  memcpy(f.F[1][0][0], qvptr->vals[1], sizex*sizey*sizez*sizeof(double));
  memcpy(f.F[2][0][0], qvptr->vals[2], sizex*sizey*sizez*sizeof(double));
  free(qvptr);
#endif // SCALAR



  qvptr = DBGetQuadvar(dbfile, "uij");
  memcpy(f.uij[0][0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  memcpy(f.uij[1][0][0], qvptr->vals[1], sizex*sizey*sizez*sizeof(double));
  memcpy(f.uij[2][0][0], qvptr->vals[2], sizex*sizey*sizez*sizeof(double));
  memcpy(f.uij[3][0][0], qvptr->vals[3], sizex*sizey*sizez*sizeof(double));
  memcpy(f.uij[4][0][0], qvptr->vals[4], sizex*sizey*sizez*sizeof(double));
  memcpy(f.uij[5][0][0], qvptr->vals[5], sizex*sizey*sizez*sizeof(double));
  free(qvptr);

  qvptr = DBGetQuadvar(dbfile, "udotij");
  memcpy(f.udotij[0][0][0], qvptr->vals[0], sizex*sizey*sizez*sizeof(double));
  memcpy(f.udotij[1][0][0], qvptr->vals[1], sizex*sizey*sizez*sizeof(double));
  memcpy(f.udotij[2][0][0], qvptr->vals[2], sizex*sizey*sizez*sizeof(double));
  memcpy(f.udotij[3][0][0], qvptr->vals[3], sizex*sizey*sizez*sizeof(double));
  memcpy(f.udotij[4][0][0], qvptr->vals[4], sizex*sizey*sizez*sizeof(double));
  memcpy(f.udotij[5][0][0], qvptr->vals[5], sizex*sizey*sizez*sizeof(double));
  free(qvptr);

  // No checkpointing support for UETCs
 
  DBClose(dbfile);

  printf0(p, "Restarted with checkpoint to step %d\n", step);

  return step;
}


/* int usable_checkpoint(hydro_fields f, hydro_params p)
 *
 * Do we have a usable checkpoint directory? Return 1 if so, 0 otherwise.
 */
int usable_checkpoint(hydro_fields f, hydro_params p)
{
  char filename[600];

  int happy = 1;

  printf0(p, "Checking for valid checkpoint\n");


  sprintf(filename, "%s/checkpoint-integrity", p.checkpointdir);

  if(access(filename,R_OK) != 0) {
    printf0(p, "Unable to access \"%s\"\n", filename);
    happy = 0;
  } else {

    sprintf(filename, "%s/checkpoint-%06d", p.checkpointdir, p.rank);


    DBfile *dbfile = DBOpen(filename, CPMODE, DB_READ);


    if(dbfile == NULL) {
      printf0(p, "File pointer is null\n");
      happy = 0;
    } else {

      int sizex, sizey, sizez, step;


      DBReadVar(dbfile, "sizex", &sizex);
      DBReadVar(dbfile, "sizey", &sizey);
      DBReadVar(dbfile, "sizez", &sizez);
      DBReadVar(dbfile, "step", &step);

      if(! ((sizex == p.slicex+2)  && (sizey == p.slicey+2) 
	    && (sizez == p.Lz)) ) {
	printf0(p, "Sizes do not match (%d vs %d), (%d vs %d), (%d vs %d)\n", 
		sizex, p.slicex+2, sizey, p.slicey+2, sizez, p.Lz);
	happy = 0;
      } else {
	printf0(p, "Looks like we can restart to step %d\n", step);
      }
      
    }

    DBClose(dbfile);
  }


  // Every rank was responsible for checking its own file, are they
  // all happy?
  return reduce_and(happy, p);

}


/* void checkpoint(hydro_fields f, hydro_params p, int step)
 *
 * Checkpoints the state of the system (except the initial Tij needed for
 * UETC calculations), for given timestep step.
 */
void checkpoint(hydro_fields f, hydro_params p, int step)
{

  DBSetCompression("METHOD=GZIP");

  char filename[600];

  MPI_Barrier(MPI_COMM_WORLD);

  double start = clock();


  if(!p.rank) {
    sprintf(filename, "%s/checkpoint-integrity", p.checkpointdir);
    printf0(p, "Removing checkpoint integrity file\n");
    if(unlink(filename) == -1) {
      printf0(p, "Problem removing checkpoint integrity file \"%s\"!\n",
	      filename);
    }
  }

  int x, i;


  /* Create a unique filename for the new Silo file */

  sprintf(filename, "%s/checkpoint-%06d", p.checkpointdir, p.rank);  

  printf0(p,"Checkpointing in %s\n", p.checkpointdir);


  DBfile *dbfile = DBCreate(filename, DB_CLOBBER, DB_LOCAL,
			    "checkpoint", CPMODE);


  DBoptlist *dboptlist = NULL;

  dboptlist = DBMakeOptlist(1);

  int col_major = DB_COLMAJOR;

  DBAddOption(dboptlist, DBOPT_MAJORORDER, &col_major);


  int *meshsize = (int *)malloc(3*sizeof(int));

  // Write metadata (magic then step)

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


  int scalardims = 1;

  DBWrite(dbfile, "sizex", &sizex, &scalardims, 1, DB_INT);
  DBWrite(dbfile, "sizey", &sizey, &scalardims, 1, DB_INT);
  DBWrite(dbfile, "sizez", &sizez, &scalardims, 1, DB_INT);
  DBWrite(dbfile, "step", &step, &scalardims, 1, DB_INT);


  DBPutQuadvar1(dbfile, "phi", "quadmesh", f.phi[0][0], meshsize, 3,
                NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);

  DBPutQuadvar1(dbfile, "pifull", "quadmesh", f.pifull[0][0], meshsize, 3,
                NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);

  DBPutQuadvar1(dbfile, "pi", "quadmesh", f.pi[0][0], meshsize, 3,
                NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);

  DBPutQuadvar1(dbfile, "phiold", "quadmesh", f.phiold[0][0], meshsize, 3,
                NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);


#ifndef SCALAR

  DBPutQuadvar1(dbfile, "T", "quadmesh", f.T[0][0], meshsize, 3,
                NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);

  DBPutQuadvar1(dbfile, "E", "quadmesh", f.E[0][0], meshsize, 3,
                NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);

  DBPutQuadvar1(dbfile, "W", "quadmesh", f.W[0][0], meshsize, 3,
                NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);

  DBPutQuadvar1(dbfile, "kappa", "quadmesh", f.kappa[0][0], meshsize, 3,
                NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);

  DBPutQuadvar1(dbfile, "p", "quadmesh", f.p[0][0], meshsize, 3,
                NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);





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




  char *zx_name = "Zx";
  char *zy_name = "Zy";
  char *zz_name = "Zz";

  char *z_names[3];

  z_names[0] = zx_name;
  z_names[1] = zy_name;
  z_names[2] = zz_name;

  double **Ztemp = (double **)malloc(3*sizeof(double *));
  Ztemp[0] = f.Z[0][0][0];
  Ztemp[1] = f.Z[1][0][0];
  Ztemp[2] = f.Z[2][0][0];
  
  DBPutQuadvar(dbfile, "Z", "quadmesh", 3, z_names, Ztemp, meshsize, 3, 
	       NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);       



  char *ux_name = "Ux";
  char *uy_name = "Uy";
  char *uz_name = "Uz";

  char *U_names[3];

  U_names[0] = ux_name;
  U_names[1] = uy_name;
  U_names[2] = uz_name;

  double **Utemp = (double **)malloc(3*sizeof(double *));
  Utemp[0] = f.U[0][0][0];
  Utemp[1] = f.U[1][0][0];
  Utemp[2] = f.U[2][0][0];
  
  DBPutQuadvar(dbfile, "U", "quadmesh", 3, U_names, Utemp, meshsize, 3, 
	       NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);       



  char *fx_name = "Fx";
  char *fy_name = "Fy";
  char *fz_name = "Fz";

  char *f_names[3];

  f_names[0] = fx_name;
  f_names[1] = fy_name;
  f_names[2] = fz_name;

  double **Ftemp = (double **)malloc(3*sizeof(double *));
  Ftemp[0] = f.F[0][0][0];
  Ftemp[1] = f.F[1][0][0];
  Ftemp[2] = f.F[2][0][0];
  
  DBPutQuadvar(dbfile, "F", "quadmesh", 3, f_names, Ftemp, meshsize, 3, 
	       NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);       

#endif // SCALAR


  char *u11_name = "U11";
  char *u21_name = "U21";
  char *u31_name = "U31";
  char *u22_name = "U22";
  char *u32_name = "U32";
  char *u33_name = "U33";

  char *u_names[6];

  u_names[0] = u11_name;
  u_names[1] = u21_name;
  u_names[2] = u31_name;
  u_names[3] = u22_name;
  u_names[4] = u32_name;
  u_names[5] = u33_name;

  double **utemp = (double **)malloc(6*sizeof(double *));
  utemp[0] = f.uij[0][0][0];
  utemp[1] = f.uij[1][0][0];
  utemp[2] = f.uij[2][0][0];
  utemp[3] = f.uij[3][0][0];
  utemp[4] = f.uij[4][0][0];
  utemp[5] = f.uij[5][0][0];
  
  DBPutQuadvar(dbfile, "uij", "quadmesh", 6, u_names, utemp, meshsize, 3, 
	       NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);       





  char *udot11_name = "Udot11";
  char *udot21_name = "Udot21";
  char *udot31_name = "Udot31";
  char *udot22_name = "Udot22";
  char *udot32_name = "Udot32";
  char *udot33_name = "Udot33";

  char *udot_names[6];

  udot_names[0] = udot11_name;
  udot_names[1] = udot21_name;
  udot_names[2] = udot31_name;
  udot_names[3] = udot22_name;
  udot_names[4] = udot32_name;
  udot_names[5] = udot33_name;

  double **udottemp = (double **)malloc(6*sizeof(double *));
  udottemp[0] = f.udotij[0][0][0];
  udottemp[1] = f.udotij[1][0][0];
  udottemp[2] = f.udotij[2][0][0];
  udottemp[3] = f.udotij[3][0][0];
  udottemp[4] = f.udotij[4][0][0];
  udottemp[5] = f.udotij[5][0][0];
  
  DBPutQuadvar(dbfile, "udotij", "quadmesh", 6,
	       udot_names, udottemp, meshsize, 3, 
	       NULL, 0, DB_DOUBLE, DB_NODECENT, dboptlist);       


  DBClose(dbfile);

  MPI_Barrier(MPI_COMM_WORLD);

  // write ghost overlap details?

  double end = clock();

  printf0(p, "checkpointing took %lf\n",
            ((double) (end - start)) / CLOCKS_PER_SEC);


  if(!p.rank) {
    sprintf(filename, "%s/checkpoint-integrity", p.checkpointdir);
    printf0(p, "Checkpointing completed on all nodes\n");
    if(fclose(fopen(filename,"w"))) {
      printf0(p, "Problem creating integrity file \"%s\"!\n", filename);
    }
  }
}  
  

