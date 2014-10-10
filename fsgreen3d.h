#ifndef _FSGREEN3D_
#define _FSGREEN3D_
void fsgreen3d(const double Cs, const double Cp,
               const double Ds, const double Dp, const double rho,
               const double* const r,
               const double* const z,
               const double* const omega, const int& nrRec,
               const int& nzRec, const int& nFreq,
               std::complex<double>* const Ug, std::complex<double>* const Sg,
               const bool calcUg, const bool calcSg);
/*   Green's function of a homogeneous fullspace in the (x,y,z) domain.
 *   Cs     Shear wave velocity.
 *   Cp     Dilatational wave velocity.
 *   Ds     Shear damping ratio.
 *   Dp     Dilatational damping ratio.
 *   rho    Density.
 *   r      Receiver location (r-coordinate) (nrRec).
 *   z      Receiver locations (z-coordinate) (nzRec).
 *   omega  Circular frequency (nFreq).
 *   nrRec  Number of r-receivers.
 *   nzRec  Number of z-receivers.
 *   nFreq  Number of frequencies.
 *   Ug     Green's displacements (3 * 3 * nxRec * nzRec * nWavey * nFreq).
 *   Sg     Green's stresses (3 * 6 * nxRec * nzRec * nWavey * nFreq).
 *   calcUg Flag to compute Ug.
 *   calcSg Flag to compute Sg.
 */
#endif
