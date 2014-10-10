function Ktil = newmarkstiff(a,K,C,M)
%NEWMARKSTIFF Equivalent stiffness matrix.
% 
%   KTIL = NEWMARKSTIFF(TIMPAR,K,C,M) Computes the equivalent stiffness
%   matrix for the Newmark analysis. K,C and M are the stiffness, damping
%   and mass matrices in the analysis and TIMPAR is a vector containing
%   coefficients.
%   
%   For the definition of this matrix, refer to:
%   K.J. Bathe,"Finite element procedures", 1996, Prentice-Hall,
%   Englewood Cliffs, New Jersey.
%   
%   This routine allows to write a customized time stepping procedure
%   in Matlab. Typically, this is performed as follows:
%   1. Choose the time step dt and the Newmark parameters alpha, delta.
%      Typically, alpha=1/4,delta=1/2 for the constant average acceleration
%      method which is unconditionally stable for linear systems.
%   2. Compute the coefficients:
%        timpar = newmarkcoef(alpha,delta,dt);
%   3. For a linear problem, the stiffness matrix is constant, so we can
%        Ktil = newmarkstiff(timpar,K,C,M);
%   4. Define initial displacements u0 and velocities v0 for time t=0. This
%      allows to compute the initial acceleration of the system through the
%      dynamic equilibrium equation.
%   5. Time stepping procedure: from time dt to time Nt*dt:
%      - Compute the equivalent force vector:
%           Ftil = newmarkforce(timpar,M,C,u0,v0,a0);
%      - Compute the displacements:
%           u1=Ktil\Ftil;
%      - Update the displacement,velocity and acceleration history:
%           [v1,a1] = newmarkupdate(a,u0,v0,a0,u1)
%   
%   See Also NEWMARKUPDATE, NEWMARFORCE, NEWMARKCOEF.

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

% Stijn François, April 2006
ndof=size(K,1);
nzk = nnz(K);
nzc = nnz(M);
nzm = nnz(C);

Ktil=spalloc(ndof,ndof,max([nzk,nzc,nzm]));
Ktil(:,:)= a(1)*M+a(2)*C+K;
