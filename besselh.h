#ifndef _BESSELH_H_
#define _BESSELH_H_

int zbesh_(double *zr,double *zi,double *fnu,int *kode,
           int *m,int *n,double *cyr,double * cyi,int *nz,int *ierr);

//         Compute a sequence of the Hankel functions H(m,a,z)
//         for superscript m=1 or 2, real nonnegative orders a=b,
//         b+1,... where b>0, and nonzero complex argument z.  A
//         scaling option is available to help avoid overflow.
//
//         On KODE=1, ZBESH computes an N member sequence of complex
//         Hankel (Bessel) functions CY(L)=H(M,FNU+L-1,Z) for super-
//         script M=1 or 2, real nonnegative orders FNU+L-1, L=1,...,
//         N, and complex nonzero Z in the cut plane -pi<arg(Z)<=pi
//         where Z=ZR+i*ZI.  On KODE=2, CBESH returns the scaled
//         functions
//
//            CY(L) = H(M,FNU+L-1,Z)*exp(-(3-2*M)*Z*i),  i**2=-1
//
//         which removes the exponential behavior in both the upper
//         and lower half planes.  Definitions and notation are found
//         in the NBS Handbook of Mathematical Functions (Ref. 1).
//
//         Input
//           ZR     - DOUBLE PRECISION real part of nonzero argument Z
//           ZI     - DOUBLE PRECISION imag part of nonzero argument Z
//           FNU    - DOUBLE PRECISION initial order, FNU>=0
//           KODE   - A parameter to indicate the scaling option
//                    KODE=1  returns
//                            CY(L)=H(M,FNU+L-1,Z), L=1,...,N
//                        =2  returns
//                            CY(L)=H(M,FNU+L-1,Z)*exp(-(3-2M)*Z*i),
//                            L=1,...,N
//           M      - Superscript of Hankel function, M=1 or 2
//           N      - Number of terms in the sequence, N>=1
//
//         Output
//           CYR    - DOUBLE PRECISION real part of result vector
//           CYI    - DOUBLE PRECISION imag part of result vector
//           NZ     - Number of underflows set to zero
//                    NZ=0    Normal return
//                    NZ>0    CY(L)=0 for NZ values of L (if M=1 and
//                            Im(Z)>0 or if M=2 and Im(Z)<0, then
//                            CY(L)=0 for L=1,...,NZ; in the com-
//                            plementary half planes, the underflows
//                            may not be in an uninterrupted sequence)
//           IERR   - Error flag
//                    IERR=0  Normal return     - COMPUTATION COMPLETED
//                    IERR=1  Input error       - NO COMPUTATION
//                    IERR=2  Overflow          - NO COMPUTATION
//                            (abs(Z) too small and/or FNU+N-1
//                            too large)
//                    IERR=3  Precision warning - COMPUTATION COMPLETED
//                            (Result has half precision or less
//                            because abs(Z) or FNU+N-1 is large)
//                    IERR=4  Precision error   - NO COMPUTATION
//                            (Result has no precision because
//                            abs(Z) or FNU+N-1 is too large)
//                    IERR=5  Algorithmic error - NO COMPUTATION
//                            (Termination condition not met)
//
// *Long Description:
//
//         The computation is carried out by the formula
//
//            H(m,a,z) = (1/t)*exp(-a*t)*K(a,z*exp(-t))
//                   t = (3-2*m)*i*pi/2
//
//         where the K Bessel function is computed as described in the
//         prologue to CBESK.
//
//         Exponential decay of H(m,a,z) occurs in the upper half z
//         plane for m=1 and the lower half z plane for m=2.  Exponential
//         growth occurs in the complementary half planes.  Scaling
//         by exp(-(3-2*m)*z*i) removes the exponential behavior in the
//         whole z plane as z goes to infinity.
//
//         For negative orders, the formula
//
//            H(m,-a,z) = H(m,a,z)*exp((3-2*m)*a*pi*i)
//
//         can be used.
//
//         In most complex variable computation, one must evaluate ele-
//         mentary functions.  When the magnitude of Z or FNU+N-1 is
//         large, losses of significance by argument reduction occur.
//         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
//         losses exceeding half precision are likely and an error flag
//         IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is double
//         precision unit roundoff limited to 18 digits precision.  Also,
//         if either is larger than U2=0.5/UR, then all significance is
//         lost and IERR=4.  In order to use the INT function, arguments
//         must be further restricted not to exceed the largest machine
//         integer, U3=I1MACH(9).  Thus, the magnitude of Z and FNU+N-1
//         is restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, and
//         U3 approximate 2.0E+3, 4.2E+6, 2.1E+9 in single precision
//         and 4.7E+7, 2.3E+15 and 2.1E+9 in double precision.  This
//         makes U2 limiting in single precision and U3 limiting in
//         double precision.  This means that one can expect to retain,
//         in the worst cases on IEEE machines, no digits in single pre-
//         cision and only 6 digits in double precision.  Similar con-
//         siderations hold for other machines.
//
//         The approximate relative error in the magnitude of a complex
//         Bessel function can be expressed as P*10**S where P=MAX(UNIT
//         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre-
//         sents the increase in error due to argument reduction in the
//         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))),
//         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF
//         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may
//         have only absolute accuracy.  This is most likely to occur
//         when one component (in magnitude) is larger than the other by
//         several orders of magnitude.  If one component is 10**K larger
//         than the other, then one can expect only MAX(ABS(LOG10(P))-K,
//         0) significant digits; or, stated another way, when K exceeds
//         the exponent of P, no significant digits remain in the smaller
//         component.  However, the phase angle retains absolute accuracy
//         because, in complex arithmetic with precision P, the smaller
//         component will not (as a rule) decrease below P times the
//         magnitude of the larger component.  In these extreme cases,
//         the principal phase angle is on the order of +P, -P, PI/2-P,
//         or -PI/2+P.

#endif
