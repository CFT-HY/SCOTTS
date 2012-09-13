/* parameters.c
 *
 * Load parameters from infile, in form "key value".
 */
#include "hydro.h"


/* hydro_params get_parameters()
 *
 * Mostly boilerplate code for parsing lines of the form "key value" 
 * from a file.
 *
 * I previously did this call-by-reference, but the parameters struct
 * (defined in hydro.h) makes the calls a lot neater, even if what is
 * being read in is a little bit more opaque here.
 *
 * For clarity, I try to keep the order in which variables are initialised
 * and parsed the same as for the struct. This helps keep this messy
 * repetitive code as tidy as possible.
 *
 * Note, empty lines and lines where the first non-space character is
 * '#' are not parsed.
 *
 * A significant outstanding issue is the limited line length due to our
 * naive use of fgets, but this is unlikely to cause trouble
 * with the sort of parameters we are using.
 */
void get_parameters(char *infile, hydro_params *p)
{


  int set_dx = 0;
  int set_dt = 0;
 
  int set_Lx = 0;
  int set_Ly = 0;
  int set_Lz = 0;

  int set_steps = 0;

  int set_Cav = 0;
  int set_C = 0;

  int set_Lheat = 0;
  int set_sigma = 0;
  int set_lcorr = 0;

  int set_interval = 0;
  int set_silointerval = 0;

  int set_initial = 0;

  int set_silodir = 0;

  int set_nucleation = 0;

  int set_bubbles = 0;

  int set_beta = 0;

  int set_scale = 0;

  char key[100];
  char value[100];
  char option[100];
  char total[300];

  int ret;
  char *retptr;

  FILE *fp;

  int i;

  fp = fopen(infile,"r");

  while(!feof(fp)) {

    // gets is dodgy, fgets less so
    retptr = fgets(total,198,fp);

    // probably EOF
    if(retptr == NULL)
      continue;

    // Not handled gracefully
    if(index(total,'\n') == NULL) {
      fprintf(stderr,"Line too long!\n");
      continue;
    }
    
    ret = sscanf(total,"%99s%99s%99s",key,value,option);


    // Not of the form "<key> <value>"
    // -- line is blank
    if((ret == EOF))
      continue;

    // Not of the form "<key> <value>"
    // -- other reason
    if (ret < 2) {
      // Get rid of the newline
      // -- doesn't leak memory, we're on the stack
      *(index(total, '\n')) = '\0';

      fprintf(stderr,"Unable to parse input line: \"%s\"\n", total);
      continue;
    }
    
    // Comment
    if(key[0] == '#') {
      continue;
    }


    // Check for parameters (in order)
    if(!strcasecmp(key,"dx")) {
      p->dx = strtod(value,NULL);
      set_dx = 1;
    }
    else if(!strcasecmp(key,"dt")) {
      p->dt = strtod(value,NULL);
      set_dt = 1;
    }    
    else if(!strcasecmp(key,"Lx")) {
      p->Lx = strtol(value,NULL,10);
      set_Lx = 1;
    } 
    else if(!strcasecmp(key,"Ly")) {
      p->Ly = strtol(value,NULL,10);
      set_Ly = 1;
    } 
    else if(!strcasecmp(key,"Lz")) {
      p->Lz = strtol(value,NULL,10);
      set_Lz = 1;
    } 
    else if(!strcasecmp(key,"steps")) {
      p->steps = strtol(value,NULL,10);
      set_steps = 1;
    }
    else if(!strcasecmp(key,"Cav")) {
      p->Cav = strtod(value,NULL);
      set_Cav = 1;
    }
    else if(!strcasecmp(key,"C")) {
      p->C = strtod(value,NULL);
      set_C = 1;
    }
    else if(!strcasecmp(key,"Lheat")) {
      p->Lheat = strtod(value,NULL);
      set_Lheat = 1;
    }
    else if(!strcasecmp(key,"sigma")) {
      p->sigma = strtod(value,NULL);
      set_sigma = 1;
    }
    else if(!strcasecmp(key,"lcorr")) {
      p->lcorr = strtod(value,NULL);
      set_lcorr = 1;
    }
    else if(!strcasecmp(key,"interval")) {
      p->interval = strtol(value,NULL,10);
      set_interval = 1;
    }
    else if(!strcasecmp(key,"silointerval")) {
      p->silointerval = strtol(value,NULL,10);
      set_silointerval = 1;
    }
    else if(!strcasecmp(key,"bubbles")) {
      p->bubbles = strtol(value,NULL,10);
      set_bubbles = 1;
    }
    else if(!strcasecmp(key,"beta")) {
      p->beta = strtod(value,NULL);
      set_beta = 1;
    }
    else if(!strcasecmp(key,"scale")) {
      p->scale = strtod(value,NULL);
      set_scale = 1;
    }
    else if(!strcasecmp(key,"initial")) {
      if(!strcasecmp(value, "shocktube")) {
	p->initial = INIT_SHOCK_TUBE;
      } else if(!strcasecmp(value, "bubble")) {
	p->initial = INIT_BUBBLE;
      } else {
	fprintf(stderr, "warning, unrecognised value for initial"
		" (%s); defaulting to bubble.\n", value);
	p->initial = INIT_BUBBLE;
      }
      set_initial = 1;
    }
    else if(!strcasecmp(key,"nucleation")) {
      if(!strcasecmp(value, "exp")) {
	p->nucleation = NUC_EXP;
      } else if(!strncasecmp(value, "list", 4)) {
	p->nucleation = NUC_LIST;
	char *curr = option;
	//	fprintf(stderr,"Curr is %s\n", curr);
	char *next = curr;

	p->n_nucsteps = 1;
	for(i=0;i<strlen(curr);i++) {
	  if( curr[i] == ',')
	    p->n_nucsteps++;
	}
	//	fprintf(stderr,"nucsteps is %d\n", p->n_nucsteps);

	i = 0;

	p->nucsteps = (int *)malloc((p->n_nucsteps)*sizeof(int));

	while(strlen(curr) && strlen(next)) {
	  p->nucsteps[i] = strtol(curr,&next,10);
	  //	  fprintf(stderr,"bubble at step %d, next is %s\n", p->nucsteps[i],next);
	  curr = next+sizeof(char);
	  //	  fprintf(stderr,"strlen next is %d and curr is %d\n", strlen(next), strlen(curr));
	  i++;
	}
	

      } else {
	fprintf(stderr, "warning, unrecognised value for nucleation"
		" (%s); defaulting to exp.\n", value);
	p->nucleation = NUC_EXP;
      }
      set_nucleation = 1;
    }
    else if(!strcasecmp(key,"silodir")) {
     
      if(strlen(value) > 500)
	fprintf(stderr,
		"Warning: silodir name \"%s\" may be too long!\n",
		value);

      strncpy(p->silodir, value, 500);
     
      set_silodir = 1;
    }


  }
  
  // Check parameters were set (in order)
  if(!set_dx) {
    fprintf(stderr, "Did not set parameter \'dx\'\n");
    exit(100);
  } else if(!set_dt) {
    fprintf(stderr, "Did not set parameter \'dt\'\n");
    exit(100);
  } else if(!set_Lx) {
    fprintf(stderr, "Did not set parameter \'Lx\'\n");
    exit(100);
  } else if(!set_Ly) {
    fprintf(stderr, "Did not set parameter \'Ly\'\n");
    exit(100);
  } else if(!set_Lz) {
    fprintf(stderr, "Did not set parameter \'Lz\'\n");
    exit(100);
  } else if(!set_steps) {
    fprintf(stderr, "Did not set parameter \'steps\'\n");
    exit(100);
  } else if(!set_Cav) {
    fprintf(stderr, "Did not set parameter \'Cav\'\n");
    exit(100);
  } else if(!set_C) {
    fprintf(stderr, "Did not set parameter \'C\'\n");
    exit(100);
  } else if(!set_Lheat) {
    fprintf(stderr, "Did not set parameter \'Lheat\'\n");
    exit(100);
  } else if(!set_sigma) {
    fprintf(stderr, "Did not set parameter \'sigma\'\n");
    exit(100);
  } else if(!set_lcorr) {
    fprintf(stderr, "Did not set parameter \'lcorr\'\n");
    exit(100);
  } else if(!set_interval) {
    fprintf(stderr, "Did not set parameter \'interval\'\n");
    exit(100);
  } else if(!set_silointerval) {
    fprintf(stderr, "Did not set parameter \'silointerval\'\n");
    exit(100);
  } else if(!set_bubbles) {
    fprintf(stderr, "Did not set parameter \'bubbles\'\n");
    exit(100);
  } else if(!set_beta) {
    fprintf(stderr, "Did not set parameter \'beta\'\n");
    exit(100);
  } else if(!set_scale) {
    fprintf(stderr, "Did not set parameter \'scale\'\n");
    exit(100);
  } else if(!set_initial) {
    fprintf(stderr, "Did not set parameter \'initial\'\n");
    exit(100);
  } else if(!set_nucleation) {
    fprintf(stderr, "Did not set parameter \'nucleation\'\n");
    exit(100);
  } else if(!set_silodir) {
    fprintf(stderr, "Did not set parameter \'silodir\'\n");
    exit(100);
  }



  if(!p->rank) {
    // Report what we found
    fprintf(stderr,"-- Read parameters from %s:\n"
	    "-- dx %g, dt %g, steps %d\n"
	    "-- Lx %d, Ly %d, Lz %d\n"
	    "-- Cav %g, C %g,\n"
	    "-- Lheat %g, sigma %g, lcorr %g\n"
	    "-- interval %d, silointerval %d, bubbles %d, beta %g, scale %g\n"
	    "-- silodir \"%s\"\n",
	    infile,
	    p->dx, p->dt, p->steps,
	    p->Lx, p->Ly, p->Lz,
	    p->Cav, p->C,
	    p->Lheat, p->sigma, p->lcorr,
	    p->interval, p->silointerval, p->bubbles, p->beta, p->scale,
	    p->silodir);
    
    if(p->initial == INIT_SHOCK_TUBE) {
      fprintf(stderr, "-- shock tube\n");
    } else if(p->initial == INIT_BUBBLE) {
      fprintf(stderr, "-- bubble\n");
    } else {
      fprintf(stderr, "-- warning, somehow have unknown initial conds.\n");
    }

    if(p->nucleation == NUC_EXP) {
      fprintf(stderr, "<Random nucleation>\n");
    } else if(p->nucleation == NUC_LIST) {
      fprintf(stderr, "<List nucleation\n");
      fprintf(stderr, "At steps: ");
      for(i=0;i<p->n_nucsteps;i++)
	fprintf(stderr, "%d, ", p->nucsteps[i]);

      fprintf(stderr,"bubbles will be nucleated>\n");
    } else {
      fprintf(stderr, "<Warning, somehow have unknown nucleation process.\n");
    }
    
  }

  fclose(fp);

}

