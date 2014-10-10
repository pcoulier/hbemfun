function [nod1,elt1,eltID]=bemmesharc(x0,R,theta,nElt,typ,typID,nod0,elt0)
%BEMMESHARC  2D arc mesh.
%   [nod,elt,eltID] = BEMMESHARC(x0,R,theta,nElt,typ,typID,nod0,elt0) meshes an
%   arc defined by its center, radius and angles with equal-sized 2D boundary 
%   elements of the specified element type.
%   
%   x0       Vector [x0 0 z0] with coordinates of the circle center.
%   R        Circle radius.
%   theta    Vector [theta0 theta1] of the start and end angles in radians.
%            The angle is defined as a rotation from the z-axis to the 
%            x-axis.
%   nElt     Number of elements.
%   typ      Element types.
%   typID    Element type ID.
%   nod0     Node array to which new nodes are added. Defaults [].
%   elt0     Element array to which new elements are added. Defaults [].
%   nod      Node array.
%   elt      Element array. 
%   eltID  Element ID's of the meshed elements.
%
%
%   See also bemmeshcircle

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

% INPUT ARGUMENT PROCESSING
if nargin<7, nod0=zeros(0,4); end
if nargin<8, elt0=zeros(0,2); end
if isempty(nod0), nod0=zeros(0,4); end
if isempty(elt0), elt0=zeros(0,2); end

[Parent,nEltNod,nCol,TypeN,TypeM,NodDef,EltDim]=bemeltdef(typID,typ);
if ~(EltDim==1), error('The element type specified is not 2D.'); end

nNod=nElt*(nEltNod-1)+1;
phi=linspace(theta(1),theta(2),nNod).';
nod=[[1:nNod].'  R*sin(phi)+x0(1) zeros(nNod,1)+x0(2) R*cos(phi)+x0(3)];

elt=zeros(nElt,nEltNod+2);
elt(:,1)=1:nElt;
elt(:,2)=typID;
for iEltNod=1:nEltNod
  elt(:,2+iEltNod)=iEltNod:(nEltNod-1):(nNod)-(nEltNod-iEltNod);
end

% CHANGE ID'S TO AVAILABLE ID'S
[nod1,elt1,dum,eltID]=bemmeshcat(nod0,elt0,nod,elt);
