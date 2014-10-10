function [nod,elt,nod2ID,elt2ID]=bemmeshrep(nod0,elt0,N,dir)
%BEMMESHREP Replicate a boundary element mesh.
%
%   [nod,elt]=BEMMESHREP(nod0,elt0,N,dir) replicates a boundary element mesh 
%   N times along the distance dir. Node and element renumbering is performed to 
%   ensure unique node and element numbers.
%
%   nod      Node array.
%   elt      Element array. 
%   N        Number of copies.
%   dir      Vector [x y z] over which copies are made.
%
bemfunlicense('VerifyOnce');
nod=nod0;
elt=elt0;
for iCopy=1:N
  nod1=nod0;
  nod1(:,2)=nod0(:,2)+iCopy/N*dir(1);
  nod1(:,3)=nod0(:,3)+iCopy/N*dir(2);
  nod1(:,4)=nod0(:,4)+iCopy/N*dir(3);
  [nod,elt]=bemmeshcat(nod,elt,nod1,elt0);
end
