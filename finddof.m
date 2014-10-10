function ind=finddof(DOF,seldof)

%FINDDOF   Find indices of specified degrees of freedom.
%  ind=finddof(DOF,seldof) returns the indices of the degrees of freedom
%  seldof in the degrees of freedom DOF. The degrees of freedom follow the
%  node.dof convention, where nod.01 to nod.03 represent the 
%  
%
%   DOF        Degrees of freedom  (nDOF * 1)
%   seldof     Selected DOF labels (kDOF * 1)
%   ind        Index vector (kDOF * 1)
%
bemfunlicense('VerifyOnce');

[L,ind]=selectdof(DOF,seldof)