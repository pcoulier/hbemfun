%BEMTIMECONV   Time convolution of Green's function.
%
%   u = BEMTIMECONV(t,ug,tBem,delt,type) computes the convolution of the
%   Green's function ug and the boundary element shape function psi. The
%   convolution is performed over the last dimension of ugt. The convolution
%   is defined as:
%
%       /+inf
%   u = |       ug(tau) * psi(tBem-tau) d tau
%       /-inf
%
%   t     Time sampling of the Green's function (1 * nTime) or (nTime * 1).
%         The vector should be monotonically increasing.
%   ug    Green's functions (... * nTime). This is usually a Green's
%         displacement or traction. The convolution is performed over
%         the last dimension of ug. The function is assumed to be zero
%         outside the sampling interval.
%   tBem  Time sampling (1 * nTimeBem) or (nTimeBem * 1) for which the
%         convolution is evaluated.
%   delt  Time step that is used in the definition of the shape function
%   type  Shape function type.
%         1: Constant shape function from -delt to 0. Used for displacements.
%         2: Triangular shape function from -delt to delt. Used for tractions
%         3: Modified triangular shape function from -delt to 0. Used for
%         tractions
%   u     Convolution (... * nTimeBem) evaluated at tBem.
