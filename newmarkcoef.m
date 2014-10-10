function timpar = newmarkcoef(alpha,delta,dt)
%NEWMARKCOEF   Newmark parameters.
%
%   TIMPAR = NEWMARKCOEF(ALPHA,DELTA,DT) Computes the coefficients
%   a0 to a7 required for the computation of the equivalent stiffness
%   matrix and equivalent force vector in the Newmark time stepping
%   procedure. ALPHA and DELTA are the Newmark parameters and DT is
%   the time step in the analysis.
%
%   For the definition of the coefficients, refer to:
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
%           Ftil = newmarkforce(timpar,K,C,M,u0,v0,a0);
%      - Compute the displacements:
%           u1=Ktil\Ftil;
%      - Update the displacement,velocity and acceleration history:
%           [v1,a1] = newmarkupdate(a,u0,v0,a0,u1)
%
%   See Also NEWMARKUPDATE, NEWMARFORCE, NEWMARKSTIFF.

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

%  Stijn François, april 2006
timpar=zeros(8,1);
timpar(1)=1/alpha/dt^2;
timpar(2)=delta/alpha/dt;
timpar(3)=1/alpha/dt;
timpar(4)=1/2/alpha-1;
timpar(5)=delta/alpha-1;
timpar(6)=dt/2*(delta/alpha-2);
timpar(7)=dt*(1-delta);
timpar(8)=delta*dt;


% EXAMPLE
%
%pars=newmarkcoef(alpha,delta,dt);
%Ktil=newmarkstiff(pars,K,C,M);
%u=zeros(length(DOF),N);
%v=zeros(length(DOF),N);
%a=zeros(length(DOF),N);
%
%u(:,1)=0;
%v(:,1)=0;
%a(:,1)=M\(-K*u1);
%
%for it=2:N
%  progbar(2,N,it);
%  Ftil=P(:,it)+newmarkforce(pars,M,C,u(:,it-1),v(:,it-1),a(:,it-1));
%  u(:,it)=Ktil\Ftil;
%  [v(:,it),a(:,it)]=newmarkupdate(pars,u(:,it-1),v(:,it-1),a(:,it-1),u(:,it));
%end