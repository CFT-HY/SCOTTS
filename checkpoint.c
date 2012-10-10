#include "hydro.h"

int load_checkpoint(hydro_fields f, hydro_params p)
{
  char filename[600];

  sprintf(filename, "%s/checkpoint-%06d", p.checkpointdir, p.rank);

  FILE *fp = fopen(filename, "r");

  char *frontmatter = "HYDRO";  
  fseek(fp, strlen(frontmatter)*sizeof(char), SEEK_SET);
  
  int sizex, sizey, sizez, step;
  
  if(fread(&sizex, sizeof(int), 1, fp) != 1) {
    printf0(p, "Error reading sizex\n");
    die(100);
  }


  if(fread(&sizey, sizeof(int), 1, fp) != 1) {
    printf0(p, "Error reading sizey\n");
    die(100);
  }

  if(fread(&sizez, sizeof(int), 1, fp) != 1) {
    printf0(p, "Error reading sizez\n");
    die(100);
  }


  if(fread(&step, sizeof(int), 1, fp) != 1) {
    printf0(p, "Error reading step\n");
    die(100);
  }


  if(fread(f.phi[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    printf0(p, "Error reading phi\n");
    die(100);
  }

  if(fread(f.pi[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    printf0(p, "Error reading pi\n");
    die(100);
  }

  if(fread(f.pifull[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    printf0(p, "Error reading pifull\n");
    die(100);
  }

  if(fread(f.T[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    printf0(p, "Error reading T\n");
    die(100);
  }

  if(fread(f.kappa[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    printf0(p, "Error reading kappa\n");
    die(100);
  }

  if(fread(f.phiold[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    printf0(p, "Error reading phiold\n");
    die(100);
  }

  if(fread(f.p[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    printf0(p, "Error reading p\n");
    die(100);
  }

  if(fread(f.E[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    printf0(p, "Error reading E\n");
    die(100);
  }

  if(fread(f.W[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    printf0(p, "Error reading W\n");
    die(100);
  }

  if(fread(f.V[0][0][0], sizeof(double), sizex*sizey*sizez*3, fp) != sizex*sizey*sizez*3) {
    printf0(p, "Error reading V\n");
    die(100);
  }

  if(fread(f.U[0][0][0], sizeof(double), sizex*sizey*sizez*3, fp) != sizex*sizey*sizez*3) {
    printf0(p, "Error reading U\n");
    die(100);
  }

  if(fread(f.F[0][0][0], sizeof(double), sizex*sizey*sizez*3, fp) != sizex*sizey*sizez*3) {
    printf0(p, "Error reading F\n");
    die(100);
  }

  if(fread(f.Z[0][0][0], sizeof(double), sizex*sizey*sizez*3, fp) != sizex*sizey*sizez*3) {
    printf0(p, "Error reading Z\n");
    die(100);
  }

  char endmatter[6];
  if(fread(endmatter, sizeof(char), 5, fp) != 5) {
    printf0(p, "Did not read endmatter\n");
    die(100);
  } else {
    if(strncmp("ORDYH",endmatter,5) != 0) {
      printf0(p, "Endmatter incorrect\n");
      die(100);
    }
  }
  if(fgetc(fp) != EOF) {
    printf0(p, "Not reached eof after endmatter\n");
    die(100);    
  }


  printf0(p, "Done loading checkpoint, at step %d\n", step);

  return step;
}



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

    FILE *fp = fopen(filename, "r");

    if(fp == NULL) {
      printf0(p, "File pointer is null\n");
      happy = 0;
    } else {

      char *frontmatter = "HYDRO";

      fseek(fp, strlen(frontmatter)*sizeof(char), SEEK_SET);

      int sizex, sizey, sizez, step;
 
      if(fread(&sizex, sizeof(int), 1, fp) != 1)
	happy = 0;

      if(fread(&sizey, sizeof(int), 1, fp) != 1)
	happy = 0;

      if(fread(&sizez, sizeof(int), 1, fp) != 1)
	happy = 0;

      if(fread(&step, sizeof(int), 1, fp) != 1)
	happy = 0;

      if(!happy) {
	printf0(p, "Short item count in header stuff\n");
      } else {

	if(! ((sizex == p.slicex+2)  && (sizey == p.slicey+2) && (sizez == p.Lz)) ) {
	  printf0(p, "Sizes do not match (%d vs %d), (%d vs %d), (%d vs %d)\n", sizex, p.slicex+2, sizey, p.slicey+2, sizez, p.Lz);
	  happy = 0;
	} else {
	  printf0(p, "Looks like we can restart to step %d\n", step);

	  
	  fseek(fp, 21*sizex*sizey*sizez*sizeof(double), SEEK_CUR);
	  char endmatter[6];
	  if(fread(endmatter, sizeof(char), 5, fp) != 5) {
	    printf0(p, "Did not read endmatter\n");
	    happy = 0;
	  } else {
	    if(strncmp("ORDYH",endmatter,5) != 0) {
	      printf0(p, "Endmatter incorrect\n");
	      happy = 0;
	    }
	  }
	  
	  if(fgetc(fp) != EOF) {
	    printf0(p, "Not reached eof after endmatter\n");
	    happy = 0;
	  }
	  
	}
	
      }
    }
  }

  //  fprintf(stderr,"Happiness %d\n", happy);

  return reduce_and(happy, p);

}

void checkpoint(hydro_fields f, hydro_params p, int step)
{

  char filename[600];

  MPI_Barrier(MPI_COMM_WORLD);

  double start = clock();



  if(!p.rank) {
    sprintf(filename, "%s/checkpoint-integrity", p.checkpointdir);
    printf0(p, "Removing checkpoint integrity file\n");
    if(unlink(filename) == -1) {
      printf0(p, "Problem removing checkpoint integrity file \"%s\"!\n", filename);
    }
  }

  int x, i;


  //  char silodir[100];
  //  sprintf(silodir,"silage-%d",(int)getpid());

  /* Create a unique filename for the new Silo file */

  sprintf(filename, "%s/checkpoint-%06d", p.checkpointdir, p.rank);  

  printf0(p,"Checkpointing in %s\n",p.checkpointdir);


  FILE *fp = fopen(filename, "w");


  // Write metadata (magic then step)

  int sizex, sizey, sizez;

  sizex = p.slicex+2;
  sizey = p.slicey+2;
  sizez = p.Lz;

  char *frontmatter = "HYDRO";

  fwrite(frontmatter,sizeof(char),strlen(frontmatter),fp);

  fwrite(&sizex, sizeof(int), 1, fp);
  fwrite(&sizey, sizeof(int), 1, fp);
  fwrite(&sizez, sizeof(int), 1, fp);
  fwrite(&step, sizeof(int), 1, fp);

  fwrite(f.phi[0][0], sizeof(double), sizex*sizey*sizez, fp);
  fwrite(f.pi[0][0], sizeof(double), sizex*sizey*sizez, fp);
  fwrite(f.pifull[0][0], sizeof(double), sizex*sizey*sizez, fp);
  fwrite(f.T[0][0], sizeof(double), sizex*sizey*sizez, fp);
  fwrite(f.kappa[0][0], sizeof(double), sizex*sizey*sizez, fp);
  fwrite(f.phiold[0][0], sizeof(double), sizex*sizey*sizez, fp);
  fwrite(f.p[0][0], sizeof(double), sizex*sizey*sizez, fp);
  fwrite(f.E[0][0], sizeof(double), sizex*sizey*sizez, fp);
  fwrite(f.W[0][0], sizeof(double), sizex*sizey*sizez, fp);

  fwrite(f.V[0][0][0], sizeof(double), sizex*sizey*sizez*3, fp);
  fwrite(f.U[0][0][0], sizeof(double), sizex*sizey*sizez*3, fp);
  fwrite(f.F[0][0][0], sizeof(double), sizex*sizey*sizez*3, fp);
  fwrite(f.Z[0][0][0], sizeof(double), sizex*sizey*sizez*3, fp);

#ifdef GW
#error not implemented yet
#endif // GW


  // Write end metadata (magic)

  char *endmatter = "ORDYH";

  fwrite(endmatter,sizeof(char),strlen(endmatter),fp);

  fclose(fp);

  MPI_Barrier(MPI_COMM_WORLD);

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
  

