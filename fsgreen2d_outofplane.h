#ifndef _FSGREEN2D_OUTOFPLANE_H_
#define _FSGREEN2D_OUTOFPLANE_H_

void fsgreen2d_outofplane(const double Cs, const double Ds, const double rho,
                          const double* const x, const double* const z,
                          const double* const omega, const int nxRec, 
                          const int nzRec, const int nFreq, 
                          std::complex<double>* const Ug, std::complex<double>* const Sg,
                          const bool calcUg, const bool calcSg);
                       
/*   Twodimensional Green's function of a homogeneous fullspace.
 *   Cs    Shear wave velocity.
 *   Cp    Dilatational wave velocity.
 *   Ds    Shear damping ratio.
 *   Dp    Dilatational damping ratio.
 *   rho   Density.
 *   x     Receiver locations (x-coordinate) (nxRec).
 *   z     Receiver locations (z-coordinate) (nzRec).
 *   omega Circular frequency (nFreq).
 *   nxRec Number of x-receivers.
 *   nzRec Number of z-receivers.
 *   ug    Green's displacements (1 * 1 * nxRec * nzRec * nFreq).
 *   sg    Green's stresses (1 * 3 * nxRec * nzRec * nFreq).
 */
#endif
