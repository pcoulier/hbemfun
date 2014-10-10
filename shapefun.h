#ifndef _SHAPEFUN_
#define _SHAPEFUN_
void shapefun(const unsigned int& ShapeType, const unsigned int& nXi,
              const double* const xi ,double* const N);
/*  Shape functions.
 *
 *  Shapetype  Shape function type index.
 *  nXi        Number of integration points.
 *  xi         Local coordinates of the integration points (nXi * eltDim).
 *  N          Resulting shape functions (nNod * nXi).
 *
 */
#endif

#ifndef _SHAPEDERIV_
#define _SHAPEDERIV_
void shapederiv(const unsigned int& ShapeType, const unsigned int& nXi,
                const double* const xi ,double* const N);
/*  Shape function derivatives
 *
 */
#endif

#ifndef _SHAPENATCOORD_
#define _SHAPENATCOORD_
void shapenatcoord(const double* const dN, const unsigned int& nNod,
                   const unsigned int& nXi, const double* const NodCoord,
                   double* const a, const unsigned int& EltDim);
/* Element Natural basis a in the integration points. This is used in furter computations
 * of the element Jacobian and the element normal in the integration points.
 *
 * dN        Derivative of shape function.
 * nNod      Number of nodes.
 * nXi       Number of integration points
 * NodCoord  Nodal Coordinates.
 * a         Element natural basis, (6 * nXi) for 3D problems or (2 * nXi) for 2D problems.
 * EltDim    Element dimension.
 */
#endif

#ifndef _JACOBIAN_
#define _JACOBIAN_
void jacobian(const double* const a, const unsigned int& nXi,double* const Jac,
              const unsigned int& EltDim);
#endif

#ifndef _TRIANGDIV_
#define _TRIANGDIV_
void triangdiv(const double* const xiSing, const unsigned int& Parent, unsigned int& nDiv,
               double* const am, double* const a1, double* const a2,
               double* const rhom, double* const rho1, double* const rho2);
#endif
