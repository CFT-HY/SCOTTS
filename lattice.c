#include "hydro.h"

/** Write a single field to a file.
 *
 * Writes the lattice field `field` to the already open file specified
 * by the file pointer `fp`. Each lattice site on the current node is
 * associated with `size` numbers - for example, `size = DIRS*SU2` for
 * an \f$\mathrm{SU}(2)\f$ gauge field.
 *
 * Transfers data 'finger by finger' (or pencil by pencil) from each
 * node to the root node, and then writes it to a single file.  The
 * resulting file can be read and used by any number of cores - the
 * layout of the resulting file does not depend on the
 * parallelisation.
 *
 * A file written in this way can be read with read_field().
 */
void write_float_field(hydro_fields f, hydro_params p, FILE *fp, float *field, int size) {

  float *finger;

  int x, y, z;

  MPI_Status status;

  finger = malloc(p.Lz * size * sizeof(float));

  for (x = 0; x < p.Lx; x++) {
    for (y = 0; y < p.Ly; y++) {

      // do we have this
      if (((x >= p.shiftx) && x < (p.shiftx + p.slicex)) &&
          ((y >= p.shifty) && y < (p.shifty + p.slicey))) {
        // yes, so are we the root node?
        if (!p.rank) {
          memcpy(finger,
                 &field[((x - p.shiftx + 1) * (p.slicey + 2) * p.Lz +
                         (y - p.shifty + 1) * p.Lz + 0) *
                        size],
                 p.Lz * size * sizeof(float));
        } else {
          MPI_Send(&field[((x - p.shiftx + 1) * (p.slicey + 2) * p.Lz +
                           (y - p.shifty + 1) * p.Lz + 0) *
                          size],
                   p.Lz * size, MPI_FLOAT, 0, x * p.Ly + y, MPI_COMM_WORLD);
        }
      } else {
        if (!p.rank) {
          MPI_Recv(finger, p.Lz * size, MPI_FLOAT,
                   (y / p.slicey) * (p.Lx / p.slicex) + x / p.slicex,
                   x * p.Ly + y, MPI_COMM_WORLD, &status);
        } else {
          // do nothing, not my job to send
        }
      }

      if (!p.rank)
        fwrite(finger, sizeof(float), p.Lz * size, fp);
    }
  }

  free(finger);
}

void write_int_field(hydro_fields f, hydro_params p, FILE *fp, int *field, int size) {

  int *finger;

  int x, y, z;

  MPI_Status status;

  finger = malloc(p.Lz * size * sizeof(int));

  for (x = 0; x < p.Lx; x++) {
    for (y = 0; y < p.Ly; y++) {

      // do we have this
      if (((x >= p.shiftx) && x < (p.shiftx + p.slicex)) &&
          ((y >= p.shifty) && y < (p.shifty + p.slicey))) {
        // yes, so are we the root node?
        if (!p.rank) {
          memcpy(finger,
                 &field[((x - p.shiftx + 1) * (p.slicey + 2) * p.Lz +
                         (y - p.shifty + 1) * p.Lz + 0) *
                        size],
                 p.Lz * size * sizeof(int));
        } else {
          MPI_Send(&field[((x - p.shiftx + 1) * (p.slicey + 2) * p.Lz +
                           (y - p.shifty + 1) * p.Lz + 0) *
                          size],
                   p.Lz * size, MPI_INTEGER, 0, x * p.Ly + y, MPI_COMM_WORLD);
        }
      } else {
        if (!p.rank) {
          MPI_Recv(finger, p.Lz * size, MPI_INTEGER,
                   (y / p.slicey) * (p.Lx / p.slicex) + x / p.slicex,
                   x * p.Ly + y, MPI_COMM_WORLD, &status);
        } else {
          // do nothing, not my job to send
        }
      }

      if (!p.rank)
        fwrite(finger, sizeof(int), p.Lz * size, fp);
    }
  }

  free(finger);
}


/** Dump the sweeping field to a file.
 */
void dump_sweep_field(hydro_fields f, hydro_params p) {
  FILE *fp;

  if (!p.rank) {
    fp = fopen("sweep_dump.dat", "w");
  } else {
    fp = NULL;
  }

  write_int_field(f, p, fp, &f.sweep[0][0][0], 1);

  if (!p.rank) {
    fclose(fp);
    fprintf(stderr, "Sweep dump field\n");
  }
}