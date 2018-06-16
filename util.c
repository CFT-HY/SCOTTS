/* util.c
 *
 * Utility functions, not currently used.
 * May be used later if we introduce more sophisticated advection.
 * They are mostly logical selection of one number from three subject
 * to certain constraints.
 */
#include "hydro.h"

Real minof3(Real a, Real b, Real c) {
  if (a < b) {
    if (b < c) {
      return a;
    } else {
      if (c < a) {
	return c;
      } else {
	return a;
      }
    }
  } else if (b < c) {
    return b;
  }

  return c;  
}


int minof3_int(int a, int b, int c) {
  if (a < b) {
    if (b < c) {
      return a;
    } else {
      if (c < a) {
	return c;
      } else {
	return a;
      }
    }
  } else if (b < c) {
    return b;
  }

  return c;  
}


Real maxof3(Real a, Real b, Real c) {
  if (a > b) {
    if (b > c) {
      return a;
    } else {
      if (c > a) {
	return c;
      } else {
	return a;
      }
    }
  } else if (b > c) {
    return b;
  }

  return c;  
}

Real minof2(Real a, Real b) {
  if (a < b)
    return a;
  return b;
}
