/* parameters.c
 *
 * Load parameters from stdin, in form "key value".
 */
#include "hydro.h"


/* hydro_params get_parameters()
 *
 * Mostly boilerplate code for parsing lines of the form "key value" 
 * from a file (currently stdin).
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
hydro_params get_parameters()
{

  hydro_params parameters;

  int set_dx = 0;
  int set_dt = 0;
 
  int set_N = 0;
  int set_steps = 0;

  int set_Cav = 0;
  int set_C = 0;

  int set_Lheat = 0;
  int set_sigma = 0;
  int set_lcorr = 0;

  int set_interval = 0;

  int set_initial = 0;

  char key[100];
  char value[100];
  char total[200];

  int ret;
  char *retptr;

  while(!feof(stdin)) {

    // gets is dodgy, fgets less so
    retptr = fgets(total,198,stdin);

    // probably EOF
    if(retptr == NULL)
      continue;

    // Not handled gracefully
    if(index(total,'\n') == NULL) {
      fprintf(stderr,"Line too long!\n");
      continue;
    }
    
    ret = sscanf(total,"%99s%99s",key,value);


    // Not of the form "<key> <value>"
    // -- line is blank
    if((ret == EOF))
      continue;

    // Not of the form "<key> <value>"
    // -- other reason
    if (ret != 2) {
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
      parameters.dx = strtod(value,NULL);
      set_dx = 1;
    }
    else if(!strcasecmp(key,"dt")) {
      parameters.dt = strtod(value,NULL);
      set_dt = 1;
    }    
    else if(!strcasecmp(key,"N")) {
      parameters.N = strtol(value,NULL,10);
      set_N = 1;
    } 
    else if(!strcasecmp(key,"steps")) {
      parameters.steps = strtol(value,NULL,10);
      set_steps = 1;
    }
    else if(!strcasecmp(key,"Cav")) {
      parameters.Cav = strtod(value,NULL);
      set_Cav = 1;
    }
    else if(!strcasecmp(key,"C")) {
      parameters.C = strtod(value,NULL);
      set_C = 1;
    }
    else if(!strcasecmp(key,"Lheat")) {
      parameters.Lheat = strtod(value,NULL);
      set_Lheat = 1;
    }
    else if(!strcasecmp(key,"sigma")) {
      parameters.sigma = strtod(value,NULL);
      set_sigma = 1;
    }
    else if(!strcasecmp(key,"lcorr")) {
      parameters.lcorr = strtod(value,NULL);
      set_lcorr = 1;
    }
    else if(!strcasecmp(key,"interval")) {
      parameters.interval = strtod(value,NULL);
      set_interval = 1;
    }
    else if(!strcasecmp(key,"initial")) {
      if(!strcasecmp(value, "shocktube")) {
	parameters.initial = INIT_SHOCK_TUBE;
      } else if(!strcasecmp(value, "bubble")) {
	parameters.initial = INIT_BUBBLE;
      } else {
	fprintf(stderr, "warning, unrecognised value for initial" \
		" (%s); defaulting to bubble.\n", value);
	parameters.initial = INIT_BUBBLE;
      }
      set_initial = 1;
    }


  }
  
  // Check parameters were set (in order)
  if(!set_dx) {
    fprintf(stderr, "Did not set parameter \'dx\'\n");
    exit(100);
  } else if(!set_dt) {
    fprintf(stderr, "Did not set parameter \'dt\'\n");
    exit(100);
  } else if(!set_N) {
    fprintf(stderr, "Did not set parameter \'N\'\n");
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
  } else if(!set_initial) {
    fprintf(stderr, "Did not set parameter \'initial\'\n");
    exit(100);
  }


  // Report what we found
  fprintf(stderr,"-- Read parameters from stdin:\n" \
	  "-- dx %g, dt %g, N %d, steps %d\n" \
	  "-- Cav %g, C %g,\n" \
	  "-- Lheat %g, sigma %g, lcorr %g\n" \
	  "-- interval %d\n",
	  parameters.dx, parameters.dt, parameters.N, parameters.steps, \
	  parameters.Cav, parameters.C, \
	  parameters.Lheat, parameters.sigma, parameters.lcorr, \
	  parameters.interval);

  if(parameters.initial == INIT_SHOCK_TUBE) {
    fprintf(stderr, "-- shock tube\n");
  } else if(parameters.initial == INIT_BUBBLE) {
    fprintf(stderr, "-- bubble\n");
  } else {
    fprintf(stderr, "-- warning, somehow have unknown initial conds.\n");
  }

  return parameters;

}

