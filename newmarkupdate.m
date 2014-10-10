function [v1,a1] = newmarkupdate(a,u0,v0,a0,u1)
%NEWMARKUPDATE update velocity and acceleration history.
%
%   [v1,a1] = newmarkupdate(a,u0,v0,a0,u1) updates the velocity and acceleration
%   history based on the previous values and the current displacement solution.
%   
%   For the definition of the Newmark scheme, refer to:
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
%   3. For a linear problem, the stiffness matrix is constant:
%        Ktil = newmarkstiff(timpar,K,C,M);
%   4. Define initial displacements u0 and velocities v0 for time t=0. This
%      allows to compute the initial acceleration of the system through the
%      dynamic equilibrium equation.
%   5. Time stepping procedure from time dt to time Nt*dt:
%      - Compute the equivalent force vector:
%           Ftil = newmarkforce(timpar,M,C,u0,v0,a0);
%      - Compute the displacements:
%           u1=Ktil\Ftil;
%      - Update the displacement,velocity and acceleration history:
%           [v1,a1] = newmarkupdate(a,u0,v0,a0,u1)
%   
%   See Also NEWMARKFORCE, NEWMARSTIFF, NEWMARKCOEF.

% Stijn François, April 2006

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

a1=a(1)*(u1-u0) -a(3)*v0 -a(4)*a0;
v1=v0+ a(7)*a0 + a(8)*a1;
