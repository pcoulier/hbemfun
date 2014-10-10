#ifndef _FSGREENF_H_
#define _FSGREENF_H_

void fsgreenf(const double Cs, const double Cp,
              const double Ds, const double Dp, const double rho,
              const double* const x,
              const double* const py,
              const double* const z,
              const double* const omega,
              const int nxRec, const int nWave,
              const int nzRec, const int nFreq,
              std::complex<double>* const Ug, std::complex<double>* const Sg,
              const bool calcUg, const bool calcSg);
/*   Green's function of a homogeneous fullspace in the (x,ky)-domain.
 *   Cs    Shear wave velocity.
 *   Cp    Dilatational wave velocity.
 *   Ds    Shear damping ratio.
 *   Dp    Dilatational damping ratio.
 *   rho   Density.
 *   x     Receiver locations (x-coordinate) (nxRec).
 *   py    Slowness (y-coordinate) (nWave).
 *         If omega ~= 0, the wavenumber sampling is given by ky = omega * py.
 *         If omega == 0, the wavenumber sampling is given by ky = py.
 *   z     Receiver locations (z-coordinate) (nzRec).
 *   omega Circular frequency (nFreq).
 *   nxRec Number of x-receivers.
 *   nWave Number of wavenumbers ky.
 *   nzRec Number of z-receivers.
 *   ug    Green's displacements (3 * 3 * nxRec * nzRec * nyWave * nFreq).
 *   sg    Green's stresses (3 * 6 * nxRec * nzRec * nyWave * nFreq).
 */
#endif
