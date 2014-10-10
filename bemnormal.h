#ifndef _BEMNORMAL_
#define _BEMNORMAL_
void bemnormal(const double* const a, const unsigned int& nXi, const unsigned int& EltDim, double* const normal);
/* Boundary element normals in integration points.
 *
 * a       Element natural basis (6 * nXi) for 3D problems or (2 * nXi) for 2D problems.
 * nXi     Number of integration points.
 * EltDim  Element dimension, 2 for 3D problems, 1 for 2D problems.
 * normal  Element normal in integration points (3 * nXi).
 *
 */

#endif
