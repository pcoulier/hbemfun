function [L,I]=unselectdof(DOF,seldof)

%UNSELECTDOF   Unselect degrees of freedom.
%
%   [L,I]=unselectdof(dof,seldof) creates the matrix to unselect degrees of 
%   freedom from the global degrees of freedom. 
%
%   dof        Degrees of freedom  (nDof * 1)
%   seldof     Unselected dof labels (ndof * 1)
%   L          Selection matrix ((nDof-ndof) * nDof)
%   I          Index vector ((nDof-ndof) * 1)

% David Dooms
% March 2008

% PREPROCESSING
bemfunlicense('VerifyOnce');

DOF=DOF(:);
seldof=seldof(:);

if ~ isempty(find(seldof==0.00))
    error('The wild card 0.00 is not allowed')
end

nDOF=length(DOF);
ndof=length(seldof);
indj=ones(1,nDOF);
for idof=1:ndof
    if floor(seldof(idof,1))==0       % wild cards 0.0X
        indjdof=find(abs(rem(DOF,1)-rem(seldof(idof,1),1))<0.0001);
    elseif rem(seldof(idof,1),1)==0   % wild cards X.00
        indjdof=find(abs(floor(DOF)-floor(seldof(idof,1)))<0.0001);
    else                                 % standard case
        indjdof=find(abs(DOF-seldof(idof,1))<0.0001);
    end
    indj(indjdof)=0;
end 
I=find(indj);
ndof=length(I);
indi=1:ndof;
s=ones(1,ndof);
L=sparse(indi,I,s,ndof,nDOF);
I=I(:);