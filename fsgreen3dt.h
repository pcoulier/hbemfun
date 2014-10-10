#ifndef _FSGREEN3DT_
#define _FSGREEN3DT_

void fsgreen3dt(const double& Cs, const double& Cp, const double& rho,
                const int& ftyp, const double& delt, const double* const r, 
                const double* const z, const double * const t, const int& nrRec, 
                const int& nzRec, const int& nTime, double* const Ug, 
                double* const Sg, const bool calcUg, const bool calcSg);
/*

  computes the displacements and stresses for the homogeneous fullspace Stokes solution 
  using analytical expressions. The Stokes solution is the convolution of 
  the Green's function with a modulation function f. The Stokes solution 
  is returned in cylindrical coordinates.
  
  Cs    Shear wave velocity.
  Cp    Dilatational wave velocity.
  rho   Density.
  ftyp  modulation function type.
        0: Unit shape function from -delt to 0.
        1: Triangular shape function from -delt to delt.
        2: Triangular shape function from -delt to 0.
  delt  Load time
  r     Receiver locations (radial coordinate) (nrRec * 1).
  z     Receiver locations (vertical coordinate) (nzRec * 1).
  t     Time [s].
  nrRec number of r-receivers
  nzRec number of z-receivers
  nTime number of Times
  Ug    Green's displacements (5 * nrRec * nzRec * nTime).
  Sg    Green's stresses (10 * nrRec * nzRec * nTime).
  calcUg Compute Green's displacements
  calcSg Compute Green's tractions
 */
#endif
