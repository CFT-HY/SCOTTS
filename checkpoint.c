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
    fprintf(stderr, "Error reading sizex\n");
    exit(100);
  }


  if(fread(&sizey, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "Error reading sizey\n");
    exit(100);
  }

  if(fread(&sizez, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "Error reading sizez\n");
    exit(100);
  }


  if(fread(&step, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "Error reading step\n");
    exit(100);
  }


  if(fread(f.phi[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    fprintf(stderr, "Error reading phi\n");
    exit(100);
  }

  if(fread(f.pi[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    fprintf(stderr, "Error reading pi\n");
    exit(100);
  }

  if(fread(f.pifull[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    fprintf(stderr, "Error reading pifull\n");
    exit(100);
  }

  if(fread(f.T[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    fprintf(stderr, "Error reading T\n");
    exit(100);
  }

  if(fread(f.kappa[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    fprintf(stderr, "Error reading kappa\n");
    exit(100);
  }

  if(fread(f.phiold[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    fprintf(stderr, "Error reading phiold\n");
    exit(100);
  }

  if(fread(f.p[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    fprintf(stderr, "Error reading p\n");
    exit(100);
  }

  if(fread(f.E[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    fprintf(stderr, "Error reading E\n");
    exit(100);
  }

  if(fread(f.W[0][0], sizeof(double), sizex*sizey*sizez, fp) != sizex*sizey*sizez) {
    fprintf(stderr, "Error reading W\n");
    exit(100);
  }

  if(fread(f.V[0][0][0], sizeof(double), sizex*sizey*sizez*3, fp) != sizex*sizey*sizez*3) {
    fprintf(stderr, "Error reading V\n");
    exit(100);
  }

  if(fread(f.U[0][0][0], sizeof(double), sizex*sizey*sizez*3, fp) != sizex*sizey*sizez*3) {
    fprintf(stderr, "Error reading U\n");
    exit(100);
  }

  if(fread(f.F[0][0][0], sizeof(double), sizex*sizey*sizez*3, fp) != sizex*sizey*sizez*3) {
    fprintf(stderr, "Error reading F\n");
    exit(100);
  }

  if(fread(f.Z[0][0][0], sizeof(double), sizex*sizey*sizez*3, fp) != sizex*sizey*sizez*3) {
    fprintf(stderr, "Error reading Z\n");
    exit(100);
  }

  char endmatter[6];
  if(fread(endmatter, sizeof(char), 5, fp) != 5) {
    fprintf(stderr, "Did not read endmatter\n");
    exit(100);
  } else {
    if(strncmp("ORDYH",endmatter,5) != 0) {
      fprintf(stderr, "Endmatter incorrect\n");
      exit(100);
    }
  }
  if(fgetc(fp) != EOF) {
    fprintf(stderr, "Not reached eof after endmatter\n");
    exit(100);    
  }


  if(!p.rank)
    fprintf(stderr, "Done loading checkpoint, at step %d\n", step);

  return step;
}



int usable_checkpoint(hydro_fields f, hydro_params p)
{
  char filename[600];

  int happy = 1;

  if(!p.rank) {
    fprintf(stderr, "Checking for valid checkpoint\n");
  }

  sprintf(filename, "%s/checkpoint-integrity", p.checkpointdir);

  if(access(filename,R_OK) != 0) {
    fprintf(stderr, "Unable to access \"%s\"\n", filename);
    happy = 0;
  } else {

    sprintf(filename, "%s/checkpoint-%06d", p.checkpointdir, p.rank);

    FILE *fp = fopen(filename, "r");

    if(fp == NULL) {
      fprintf(stderr, "File pointer is null\n");
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
	fprintf(stderr, "Short item count in header stuff\n");
      } else {

	if(! ((sizex == p.slicex+2)  && (sizey == p.slicey+2) && (sizez == p.Lz)) ) {
	  fprintf(stderr, "Sizes do not match (%d vs %d), (%d vs %d), (%d vs %d)\n", sizex, p.slicex+2, sizey, p.slicey+2, sizez, p.Lz);
	  happy = 0;
	} else {
	  if(!p.rank) {
	    fprintf(stderr, "Looks like we can restart to step %d\n", step);
	  }
	  
	  fseek(fp, 21*sizex*sizey*sizez*sizeof(double), SEEK_CUR);
	  char endmatter[6];
	  if(fread(endmatter, sizeof(char), 5, fp) != 5) {
	    fprintf(stderr, "Did not read endmatter\n");
	    happy = 0;
	  } else {
	    if(strncmp("ORDYH",endmatter,5) != 0) {
	      fprintf(stderr, "Endmatter incorrect\n");
	      happy = 0;
	    }
	  }
	  
	  if(fgetc(fp) != EOF) {
	    fprintf(stderr, "Not reached eof after endmatter\n");
	    happy = 0;
	  }
	  
	}
	
      }
    }
  }

  fprintf(stderr,"Happiness %d\n", happy);

  return reduce_and(happy, p);

}

void checkpoint(hydro_fields f, hydro_params p, int step)
{

  char filename[600];

  MPI_Barrier(MPI_COMM_WORLD);

  double start = clock();



  if(!p.rank) {
    sprintf(filename, "%s/checkpoint-integrity", p.checkpointdir);
    fprintf(stderr, "Removing checkpoint integrity file\n");
    if(unlink(filename) == -1) {
      fprintf(stderr, "Problem removing checkpoint integrity file \"%s\"!\n", filename);
    }
  }

  int x, i;


  //  char silodir[100];
  //  sprintf(silodir,"silage-%d",(int)getpid());

  /* Create a unique filename for the new Silo file */

  sprintf(filename, "%s/checkpoint-%06d", p.checkpointdir, p.rank);  

  if(!p.rank)
    fprintf(stderr,"Checkpointing in %s\n",p.checkpointdir);


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

  if(!p.rank)
    fprintf(stderr,"checkpointing took %lf\n",
            ((double) (end - start)) / CLOCKS_PER_SEC);


  if(!p.rank) {
    sprintf(filename, "%s/checkpoint-integrity", p.checkpointdir);
    fprintf(stderr, "Checkpointing completed on all nodes\n");
    if(fclose(fopen(filename,"w"))) {
      fprintf(stderr, "Problem creating integrity file \"%s\"!\n", filename);
    }
  }
}  
  

