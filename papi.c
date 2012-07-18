#include "hydro.h"

#ifdef PAPI

// Naughty global variable!
int EventSet = PAPI_NULL;

void papi_init()
{

  int retval;


  fprintf(stderr, "Initialising PAPI counters...\n");

  retval = PAPI_library_init(PAPI_VER_CURRENT);

  if (retval != PAPI_VER_CURRENT && retval > 0) {
    fprintf(stderr,"PAPI library version mismatch!\n");
    exit(1);
  }

  if (retval < 0) {
    fprintf(stderr, "Initialization error!\n");
    exit(1);
  }

  /* Create the Event Set */
  if (PAPI_create_eventset(&EventSet) != PAPI_OK)
    exit(1);

  /* Add Total Instructions Executed to our Event Set */
  if (PAPI_add_event(EventSet, PAPI_TOT_INS) != PAPI_OK) {
    fprintf(stderr, "Error adding PAPI_TOT_INS to event set!\n");
    exit(1);
  }

  if (PAPI_add_event(EventSet, PAPI_L1_TCM) != PAPI_OK) {
    fprintf(stderr, "Error adding PAPI_L1_TCM to event set!\n");
    exit(1);
  }

  if (PAPI_add_event(EventSet, PAPI_L2_TCM) != PAPI_OK) {
    fprintf(stderr, "Error adding PAPI_L2_TCM to event set!\n");
    exit(1);
  }

  if (PAPI_add_event(EventSet, PAPI_L3_TCM) != PAPI_OK) {
    fprintf(stderr, "Error adding PAPI_L3_TCM to event set!\n");
    exit(1);
  }

  /* Start counting events in the Event Set */
  if (PAPI_start(EventSet) != PAPI_OK) {
    fprintf(stderr, "Error calling PAPI_start with chosen event set!\n");
    exit(1);
  }


}

void papi_finalise()
{

  long_long values[4];

  fprintf(stderr, "Finalising PAPI counters...\n");


  /* Stop the counting of events in the Event Set */
  if (PAPI_stop(EventSet, values) != PAPI_OK) {
    fprintf(stderr, "Did not safely stop PAPI!\n");
    exit(1);
  }

  fprintf(stderr, "Total instructions: %lld\n", values[0]);
  fprintf(stderr, "L1 misses: %lld\n", values[1]);
  fprintf(stderr, "L2 misses: %lld\n", values[2]);
  fprintf(stderr, "L3 misses: %lld\n", values[3]);
}

#endif // PAPI
