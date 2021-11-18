
/** @file transport.c
 *
 * This file contains the functions that perform advection.
 *
 * We have two versions, depending on whether the field being advected
 * lives on lattice sites or between them.
 *
 * Currently we use simple donor cell advection, with compiler flags
 * enabling 2nd order flux reconstruction with a corresponding flux limiter
 * for advection.
 * We may need to look into fancier
 * advection methods to avoid getting strong shocks.
 */
#include "hydro.h"

// Not sure if I need to set FLX to zero if not defined.
// When testing I could actually just put e.g VANLEER in place of FL1 in
// expression below but I worry that might be compiler dependent.
// Maybe this should go in hydro.h
#ifdef VANLEER
#define FL1 1
#endif
#ifdef SUPERBEE
#define FL2 1
#endif
#ifdef MINMOD
#define FL3 1
#endif
#ifdef OSPRE
#define FL4 1
#endif
#ifdef VANALBADA
#define FL5 1
#endif
#ifdef MONOCENT
#define FL6 1
#endif

#if (FL1 + FL2 + FL3 + FL4 + FL5 + FL6 == 1)
#define TRANSPORT
#elif (FL1 + FL2 + FL3 + FL4 + FL5 + FL6 > 1)
#error "Only one flux limiter compilation option allowed."
#endif

/** Donor cell advection for `f.E` in direction `dir`.
 *
 * Does nothing ifdef `SCALAR`.
 *
 * Simple donor cell i.e the value of `f.E` at the boundary is taken
 * to be the same as the value in the body of the cell. It may be more
 * consistent to swap to piecewise-linear advection.
 */
void donor_E_dir(hydro_fields f, hydro_params p, int dir)
{

#ifndef SCALAR
  int x, y, z;

  //  halo_field(f.V[dir],p);
  //  halo_field(f.E,p);

  int dx = 0;
  int dy = 0;
  int dz = 0;

  if (dir == 0)
    dx = 1;
  if (dir == 1)
    dy = 1;
  if (dir == 2)
    dz = 1;

  for (x = 1; x <= p.slicex; x++) {
    for (y = 1; y <= p.slicey; y++) {
      for (z = 0; z < p.Lz; z++) {
        if (f.V[dir][x][y][z] >= 0.0)
          f.F[dir][x][y][z]
              = f.V[dir][x][y][z] * f.E[x - dx][y - dy][(z - dz + p.Lz) % p.Lz];
        else
          f.F[dir][x][y][z] = f.V[dir][x][y][z] * f.E[x][y][z];
      }
    }
  }

  halo_field(f.F[dir], p);

  for (x = 1; x <= p.slicex; x++) {
    for (y = 1; y <= p.slicey; y++) {
      for (z = 0; z < p.Lz; z++) {
        f.E[x][y][z]
            = f.E[x][y][z]
              + p.dt
                    * (f.F[dir][x][y][z]
                       - f.F[dir][x + dx][y + dy][(z + dz + p.Lz) % p.Lz])
                    / p.dx;
      }
    }
  }

  halo_field(f.E, p);
#endif
}

/** Donor cell advection for `f.Z` in direction `dir`.
 *
 * Does nothing ifdef `SCALAR`.
 *
 *
 * Simple donor cell i.e the value of `f.Z` at the boundary is taken
 * to be the same as the value in the body of the cell. It may be more
 * consistent to swap to piecewise-linear advection.
 *
 * Needs to calculate velocity in the body of the cell as that is
 * where the interfaces for the momentum flux live.
 *
 * Advects each component of Z in the direction `dir`.
 */
void donor_Z_dir(hydro_fields f, hydro_params p, int dir)
{

#ifndef SCALAR

  int x, y, z, i;
  // Used to define the other three cells for face averaging in direction
  // general manner.
  int xcell1, xcell2, xcell3;
  int ycell1, ycell2, ycell3;
  int zcell1, zcell2, zcell3;

  int dx = 0;
  int dy = 0;
  int dz = 0;

  if (dir == 0) {
    dx = 1;
  }
  if (dir == 1) {
    dy = 1;
  }
  if (dir == 2) {
    dz = 1;
  }

  float Vface;

  for (x = 1; x <= p.slicex; x++) {
    xcell1 = x - (dy + dz);
    xcell2 = x;
    xcell3 = x - (dy + dz);
    for (y = 1; y <= p.slicey; y++) {
      ycell1 = y - (dx);
      ycell2 = y - (dz);
      ycell3 = y - (dx + dz);
      for (z = 0; z < p.Lz; z++) {
        zcell1 = z;
        zcell2 = z - (dx + dy);
        zcell3 = z - (dx + dy);

        // need to find V on node
        // centered face instead of zone centered.
        Vface = 0.125
                * ((f.V[dir][x][y][z]
                    + f.V[dir][x - dx][y - dy][(z - dz + p.Lz) % p.Lz])
                   + (f.V[dir][xcell1][ycell1][(zcell1 + p.Lz) % p.Lz]
                      + f.V[dir][xcell1 - dx][ycell1 - dy]
                           [(zcell1 - dz + p.Lz) % p.Lz])
                   + (f.V[dir][xcell2][ycell2][(zcell2 + p.Lz) % p.Lz]
                      + f.V[dir][xcell2 - dx][ycell2 - dy]
                           [(zcell2 - dz + p.Lz) % p.Lz])
                   + (f.V[dir][xcell3][ycell3][(zcell3 + p.Lz) % p.Lz]
                      + f.V[dir][xcell3 - dx][ycell3 - dy]
                           [(zcell3 - dz + p.Lz) % p.Lz]));

        for (i = 0; i < 3; i++) {
          if (Vface >= 0.0) {
            f.F[i][x][y][z]
                = Vface * f.Z[i][x - dx][y - dy][(z - dz + p.Lz) % p.Lz];
          } else {
            f.F[i][x][y][z] = Vface * f.Z[i][x][y][z];
          }
        }
      }
    }
  }

  halo_field(f.F[0], p);
  halo_field(f.F[1], p);
  halo_field(f.F[2], p);

  for (x = 1; x <= p.slicex; x++) {
    for (y = 1; y <= p.slicey; y++) {
      for (z = 0; z < p.Lz; z++) {
        for (i = 0; i < 3; i++) {
          f.Z[i][x][y][z] = (f.Z[i][x][y][z]
                             - p.dt
                                   * (f.F[i][x + dx][y + dy][(z + dz) % p.Lz]
                                      - f.F[i][x][y][z])
                                   / p.dx);
        }
      }
    }
  }
  halo_field(f.Z[0], p);
  halo_field(f.Z[1], p);
  halo_field(f.Z[2], p);

#endif
}

/** Flux limiter, multiple options that can be specified using compiler flags.
 * An attempt to go beyond the Van Leer limiter from Anninos and Fragile.
 *
 */
static inline float flux_limiter(float r)
{
#if defined VANLEER
  if (r > 0) {
    return 2 * r / (1 + r);
  } else {
    return 0;
  }
#elif defined SUPERBEE
  return fmaxf(fmaxf(0, fminf(2 * r, 1)), fminf(r, 2));
#elif defined MINMOD
  return fmaxf(0, fminf(1, r));
#elif defined OSPRE
  return 1.5 * (r * (r + 1)) / (r * (r + 1) + 1);
#elif defined VANALBADA
  return r * (r + 1) / (r * r + 1);
#elif defined MONOCENT
  return fmaxf(0, fminf(0.5 * (r + 1), fminf(2 * r, 2)));
#else
  fprintf(stderr, "In flux limiter but no flux scheme defined.\n"
                  "Should never reach here, dying\n");
  die(101);
  return -1;
#endif
}

/** Advection for `f.E` in direction `dir` using monotonic interpolation for the
 * flux at the boundary.
 *
 * Does nothing ifdef `SCALAR`.
 *
 * Advects using a similar algorithm to that described in
 * Anninos and Fragile (2002), but now with additional flux limiter options.
 */
void transport_E_dir(hydro_fields f, hydro_params p, int dir)
{
#ifndef SCALAR
  int x, y, z;
  float r;
  float phi;
  float*** delta = make_field(p);
  // Guarding against division by zero when constructing the flux limiter
  // argument r. Either because Delta is exactly zero or denormalised.
  float epsilon = 1e-20;
  //  halo_field(f.V[dir],p);
  //  halo_field(f.E,p);
  int dx = 0;
  int dy = 0;
  int dz = 0;

  if (dir == 0)
    dx = 1;
  if (dir == 1)
    dy = 1;
  if (dir == 2)
    dz = 1;

  // Construct backwards derivative of E
  for (x = 1; x <= p.slicex; x++) {
    for (y = 1; y <= p.slicey; y++) {
      for (z = 0; z < p.Lz; z++) {
        delta[x][y][z]
            = (f.E[x][y][z] - f.E[x - dx][y - dy][(z - dz + p.Lz) % p.Lz])
              / p.dx;
      }
    }
  }

  halo_field(delta, p);

  // This loop constructs the terms in the brackets of (3-12) using (3-13)
  // from Anninos and Fragile (2002).
  // Our F_i corresponds to \tilde{D}_i*V_i. We allow for alternative flux
  // limiters to just Van Leer.

  // This should be equivalent to the notes in section 4.4.3 of the advection
  // notes included in doc folder.
  for (x = 1; x <= p.slicex; x++) {
    for (y = 1; y <= p.slicey; y++) {
      for (z = 0; z < p.Lz; z++) {
        if (f.V[dir][x][y][z] > 0) {
          r = delta[x - dx][y - dy][(z - dz + p.Lz) % p.Lz]
              / (delta[x][y][z] + epsilon);
          phi = flux_limiter(r);
          f.F[dir][x][y][z] = (f.V[dir][x][y][z]
                               * (f.E[x - dx][y - dy][(z - dz + p.Lz) % p.Lz]
                                  + 0.5 * (p.dx - f.V[dir][x][y][z] * p.dt) * phi
                                        * delta[x][y][z]));
        } else {
          r = delta[x + dx][y + dy][(z + dz) % p.Lz]
              / (delta[x][y][z] + epsilon);
          phi = flux_limiter(r);
          f.F[dir][x][y][z] = (f.V[dir][x][y][z]
                               * (f.E[x][y][z]
                                  - 0.5 * (p.dx + f.V[dir][x][y][z] * p.dt) * phi
                                        * delta[x][y][z]));
        }
      }
    }
  }

  halo_field(f.F[dir], p);

  for (x = 1; x <= p.slicex; x++) {
    for (y = 1; y <= p.slicey; y++) {
      for (z = 0; z < p.Lz; z++) {
        f.E[x][y][z]
            = f.E[x][y][z]
              + p.dt
                    * (f.F[dir][x][y][z]
                       - f.F[dir][x + dx][y + dy][(z + dz + p.Lz) % p.Lz])
                    / p.dx;
      }
    }
  }

  halo_field(f.E, p);

  free_field(p, delta);
#endif
}

/** Advection for `f.Z` in direction `dir` using monotonic interpolation for the
 * flux at the boundary.
 *
 * Does nothing ifdef `SCALAR`.
 *
 * Advects using a similar  algorithm to that described in
 * Anninos and Fragile (2002), but now with additional flux limiters.
 * This advects each component of Z in the direction `dir`.
 */
void transport_Z_dir(hydro_fields f, hydro_params p, int dir)
{

#ifndef SCALAR
  int x, y, z, i;
  // Used to define the other three cells for face averaging in direction
  // general manner.
  int xcell1, xcell2, xcell3;
  int ycell1, ycell2, ycell3;
  int zcell1, zcell2, zcell3;

  float r;
  float phi;
  // Guarding against division by zero when constructing the flux limiter
  // argument r. Either because Delta is exactly zero or denormalised.
  float epsilon = 1e-20;
  int dx = 0;
  int dy = 0;
  int dz = 0;

  if (dir == 0) {
    dx = 1;
  }
  if (dir == 1) {
    dy = 1;
  }
  if (dir == 2) {
    dz = 1;
  }

  float Vface;
  float**** delta = make_vector(p);

  // Construct backwards derivative.
  for (x = 1; x <= p.slicex; x++) {
    for (y = 1; y <= p.slicey; y++) {
      for (z = 0; z < p.Lz; z++) {
        for (i = 0; i < 3; i++) {
          delta[i][x][y][z] = (f.Z[i][x][y][z]
                               - f.Z[i][x - dx][y - dy][(z - dz + p.Lz) % p.Lz])
                              / p.dx;
        }
      }
    }
  }

  halo_field(delta[0], p);
  halo_field(delta[1], p);
  halo_field(delta[2], p);
  for (x = 1; x <= p.slicex; x++) {
    xcell1 = x - (dy + dz);
    xcell2 = x;
    xcell3 = x - (dy + dz);
    for (y = 1; y <= p.slicey; y++) {
      ycell1 = y - (dx);
      ycell2 = y - (dz);
      ycell3 = y - (dx + dz);
      for (z = 0; z < p.Lz; z++) {
        zcell1 = z;
        zcell2 = z - (dx + dy);
        zcell3 = z - (dx + dy);

        // Need to find V on node
        // centered face instead of zone centered.
        Vface = 0.125
                * ((f.V[dir][x][y][z]
                    + f.V[dir][x - dx][y - dy][(z - dz + p.Lz) % p.Lz])
                   + (f.V[dir][xcell1][ycell1][(zcell1 + p.Lz) % p.Lz]
                      + f.V[dir][xcell1 - dx][ycell1 - dy]
                           [(zcell1 - dz + p.Lz) % p.Lz])
                   + (f.V[dir][xcell2][ycell2][(zcell2 + p.Lz) % p.Lz]
                      + f.V[dir][xcell2 - dx][ycell2 - dy]
                           [(zcell2 - dz + p.Lz) % p.Lz])
                   + (f.V[dir][xcell3][ycell3][(zcell3 + p.Lz) % p.Lz]
                      + f.V[dir][xcell3 - dx][ycell3 - dy]
                           [(zcell3 - dz + p.Lz) % p.Lz]));

        if (Vface >= 0.0) {
          r = delta[i][x - dx][y - dy][(z - dz + p.Lz) % p.Lz]
              / (delta[i][x][y][z] + epsilon);
          phi = flux_limiter(r);
          f.F[i][x][y][z]
              = (Vface
                 * (f.Z[i][x - dx][y - dy][(z - dz + p.Lz) % p.Lz]
                    + 0.5 * (p.dx - Vface * p.dt) * phi * delta[i][x][y][z]));
        } else {
          r = delta[i][x + dx][y + dy][(z + dz) % p.Lz]
              / (delta[i][x][y][z] + epsilon);
          phi = flux_limiter(r);
          f.F[i][x][y][z]
              = (Vface
                 * (f.Z[i][x][y][z]
                    - 0.5 * (p.dx + Vface * p.dt) * phi * delta[i][x][y][z]));
        }
      }
    }
  }

  halo_field(f.F[0], p);
  halo_field(f.F[1], p);
  halo_field(f.F[2], p);

  // Eq (3-12)
  for (x = 1; x <= p.slicex; x++) {
    for (y = 1; y <= p.slicey; y++) {
      for (z = 0; z < p.Lz; z++) {
        for (i = 0; i < 3; i++) {
          f.Z[i][x][y][z] = (f.Z[i][x][y][z]
                             - p.dt
                                   * (f.F[i][x + dx][y + dy][(z + dz) % p.Lz]
                                      - f.F[i][x][y][z])
                                   / p.dx);
        }
      }
    }
  }
  halo_field(f.Z[0], p);
  halo_field(f.Z[1], p);
  halo_field(f.Z[2], p);

  free_vector(p, delta);

#endif
}

/** In W&M, transport for momentum is done by using the inertial density
 * transfer to reconstruct the momentum transfer.
 *
 * We use second order flux reconstruction as in the other transport functions.
 * Must call after advecting E in all directions so that f.F is populated.
 * After calling this function f.F will be populated with flux of kappa*E in
 * direction dir.
 *
 * I think if doing half step updates that velocities need to be updated after
 * calling this for all directions.
 */
void transport_Z_WM_dir(hydro_fields f, hydro_params p, int dir)
{
#ifndef SCALAR

  float phi;
  float r;
  // Guarding against division by zero when constructing the flux limiter
  // argument r. Either because Delta is exactly zero or denormalised.
  float epsilon = 1e-20;

  int x, y, z, i;
  int dx = 0;
  int dy = 0;
  int dz = 0;

  // Used to define the other three cells for face averaging in direction
  // general manner.
  int xcell1, xcell2, xcell3;
  int ycell1, ycell2, ycell3;
  int zcell1, zcell2, zcell3;

  // U^2 at x,y,z, and at x,y,z shifted one site in the negative direction dir.
  float Usq, Usqminus;

  float**** delta = make_vector(p);
  // Transfer of Z constructed from F_node and interpolation of U.
  float**** F_node_Z = make_vector(p);
  // Flux of kappa*E in direction dir on node centered face
  float F_node;
  // Velocity in direction dir on node centered face
  float Vface;

  // Assume F is populated with flux of E. Want to turn this into flux of
  // \kappa E. Might as well reuse field rather than requiring even more
  // temporary arrays.
  if (dir == 0) {
    dx = 1;
  }
  if (dir == 1) {
    dy = 1;
  }
  if (dir == 2) {
    dz = 1;
  }
  for (x = 1; x <= p.slicex; x++) {
    for (y = 1; y <= p.slicey; y++) {
      for (z = 0; z < p.Lz; z++) {
        f.F[dir][x][y][z]
            = 0.5 * f.F[dir][x][y][z]
              * (f.kappa[x][y][z]
                 + f.kappa[x - dx][y - dy][(z - dz + p.Lz) % p.Lz]);

        // Will also need backward difference of U later:
        // Construct backwards derivative of U[i].
        for (i = 0; i < 3; i++) {
          delta[i][x][y][z] = (f.U[i][x][y][z]
                               - f.U[i][x - dx][y - dy][(z - dz + p.Lz) % p.Lz])
                              / p.dx;
        }
      }
    }
  }

  halo_field(f.F[dir], p);
  halo_field(delta[0], p);
  halo_field(delta[1], p);
  halo_field(delta[2], p);

  for (x = 1; x <= p.slicex; x++) {
    xcell1 = x - (dy + dz);
    xcell2 = x;
    xcell3 = x - (dy + dz);
    for (y = 1; y <= p.slicey; y++) {
      ycell1 = y - (dx);
      ycell2 = y - (dz);
      ycell3 = y - (dx + dz);
      for (z = 0; z < p.Lz; z++) {
        zcell1 = z;
        zcell2 = z - (dx + dy);
        zcell3 = z - (dx + dy);

        // Now we need to construct node face centered version of f.F in
        // direction dir (not including
        // U[i] yet).
        F_node = 0.125
                 * ((f.F[dir][x][y][z]
                     + f.F[dir][x - dx][y - dy][(z - dz + p.Lz) % p.Lz])
                    + (f.F[dir][xcell1][ycell1][(zcell1 + p.Lz) % p.Lz]
                       + f.F[dir][xcell1 - dx][ycell1 - dy]
                            [(zcell1 - dz + p.Lz) % p.Lz])
                    + (f.F[dir][xcell2][ycell2][(zcell2 + p.Lz) % p.Lz]
                       + f.F[dir][xcell2 - dx][ycell2 - dy]
                            [(zcell2 - dz + p.Lz) % p.Lz])
                    + (f.F[dir][xcell3][ycell3][(zcell3 + p.Lz) % p.Lz]
                       + f.F[dir][xcell3 - dx][ycell3 - dy]
                            [(zcell3 - dz + p.Lz) % p.Lz]));

        // Also need to find V on node
        // centered face instead of zone centered.
        // Pick one of the methods below, W&M advocate for method 2.
        // Method 1: average over three velocities
        /*
              Vface = 0.125
                               * ((f.V[dir][x][y][z]
                                   + f.V[dir][x - dx][y - dy][(z - dz + p.Lz) %
           p.Lz])
                                  + (f.V[dir][xcell1][ycell1][(zcell1 + p.Lz) %
           p.Lz]
                                     + f.V[dir][xcell1 - dx][ycell1 - dy]
                                          [(zcell1 - dz + p.Lz) % p.Lz])
                                  + (f.V[dir][xcell2][ycell2][(zcell2 + p.Lz) %
           p.Lz]
                                     + f.V[dir][xcell2 - dx][ycell2 - dy]
                                          [(zcell2 - dz + p.Lz) % p.Lz])
                                  + (f.V[dir][xcell3][ycell3][(zcell3 + p.Lz) %
           p.Lz]
                                     + f.V[dir][xcell3 - dx][ycell3 - dy]
                                          [(zcell3 - dz + p.Lz) % p.Lz]));
        */
        // Method 2: Construct from 4 velocities.

        Usq = (f.U[0][x][y][z] * f.U[0][x][y][z]
               + f.U[1][x][y][z] * f.U[1][x][y][z]
               + f.U[2][x][y][z] * f.U[2][x][y][z]);
        Usqminus = (f.U[0][x - dx][y - dy][(z - dz + p.Lz) % p.Lz]
                        * f.U[0][x - dx][y - dy][(z - dz + p.Lz) % p.Lz]
                    + f.U[1][x - dx][y - dy][(z - dz + p.Lz) % p.Lz]
                          * f.U[1][x - dx][y - dy][(z - dz + p.Lz) % p.Lz]
                    + f.U[2][x - dx][y - dy][(z - dz + p.Lz) % p.Lz]
                          * f.U[2][x - dx][y - dy][(z - dz + p.Lz) % p.Lz]);
        Vface = 0.5
                * (f.U[dir][x][y][z] / sqrt(1 - Usq)
                   + f.U[dir][x - dx][y - dy][(z - dz + p.Lz) % p.Lz]
                         / sqrt(1 - Usqminus));

        // Construct full flux for Z from F_node and interpolated U[i].
        // F_node_Z[i] will be(Eq 2.130) in W&M shifted one site in the negative
        // direction dir and multiplied by dx^2/dt.

        for (i = 0; i < 3; i++) {

          if (Vface >= 0.0) {
            r = delta[i][x - dx][y - dy][(z - dz + p.Lz) % p.Lz]
                / (delta[i][x][y][z] + epsilon);
            phi = flux_limiter(r);
            F_node_Z[i][x][y][z]
                = F_node
                  * (f.U[i][x - dx][y - dy][(z - dz + p.Lz) % p.Lz]
                     + 0.5 * (p.dx - Vface * p.dt) * phi * delta[i][x][y][z]);
          } else {
            r = delta[i][x + dx][y + dy][(z + dz) % p.Lz]
                / (delta[i][x][y][z] + epsilon);
            phi = flux_limiter(r);
            F_node_Z[i][x][y][z]
                = F_node
                  * (f.U[i][x][y][z]
                     - 0.5 * (p.dx + Vface * p.dt) * phi * delta[i][x][y][z]);
          }
        }
      }
    }
  }

  halo_field(F_node_Z[0], p);
  halo_field(F_node_Z[1], p);
  halo_field(F_node_Z[2], p);
  for (x = 1; x <= p.slicex; x++) {
    for (y = 1; y <= p.slicey; y++) {
      for (z = 0; z < p.Lz; z++) {
        for (i = 0; i < 3; i++) {
          f.Z[i][x][y][z]
              = (f.Z[i][x][y][z]
                 - p.dt
                       * (F_node_Z[i][x + dx][y + dy][(z + dz) % p.Lz]
                          - F_node_Z[i][x][y][z])
                       / p.dx);
        }
      }
    }
  }

  halo_field(f.Z[0], p);
  halo_field(f.Z[1], p);
  halo_field(f.Z[2], p);

  free_vector(p, F_node_Z);
  free_vector(p, delta);

#endif
}

/** Advection for `E` in all directions.
 *
 *  Note that this is not what the Grant and Wilson
 * suggest, rather they state the update should be split into two half
 * steps, and the order should be reversed for each half.
 */
void advect_E(hydro_fields f, hydro_params p, int adv_order)
{

  int directions[18] = { 0, 1, 2, 2, 1, 0, 2, 0, 1, 1, 0, 2, 1, 2, 0, 0, 2, 1 };

#ifdef TRANSPORT
  transport_E_dir(f, p, directions[(adv_order % 6) * 3]);
  transport_E_dir(f, p, directions[(adv_order % 6) * 3 + 1]);
  transport_E_dir(f, p, directions[(adv_order % 6) * 3 + 2]);
#else
  donor_E_dir(f, p, directions[(adv_order % 6) * 3]);
  donor_E_dir(f, p, directions[(adv_order % 6) * 3 + 1]);
  donor_E_dir(f, p, directions[(adv_order % 6) * 3 + 2]);
#endif
}

/** Advection for `f.Z` in all directions.
 *
 * Note that this is not what the Grant and Wilson
 * suggest, rather they state the update should be split into two half
 * steps, and the order of advection should be reversed for each half.
 */
void advect_Z(hydro_fields f, hydro_params p, int adv_order)
{

  int directions[18] = { 0, 1, 2, 2, 1, 0, 2, 0, 1, 1, 0, 2, 1, 2, 0, 0, 2, 1 };

#ifdef TRANSPORT
#ifdef WMMOMADVECT
  transport_Z_WM_dir(f, p, directions[(adv_order % 6) * 3]);
  transport_Z_WM_dir(f, p, directions[(adv_order % 6) * 3 + 1]);
  transport_Z_WM_dir(f, p, directions[(adv_order % 6) * 3 + 2]);
#else
  transport_Z_dir(f, p, directions[(adv_order % 6) * 3]);
  transport_Z_dir(f, p, directions[(adv_order % 6) * 3 + 1]);
  transport_Z_dir(f, p, directions[(adv_order % 6) * 3 + 2]);
#endif // WMMOMADVECT
#else
#ifdef WMMOMADVECT
#error "WMMOMADVECT not implemented yet for donor cell"
#endif // WMMOMADVECT
  donor_Z_dir(f, p, directions[(adv_order % 6) * 3]);
  donor_Z_dir(f, p, directions[(adv_order % 6) * 3 + 1]);
  donor_Z_dir(f, p, directions[(adv_order % 6) * 3 + 2]);
#endif
}
