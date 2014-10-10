function [nod1,elt1,eltID]=bemmeshcircle(x0,R,nElt,typ,typID,orient,nod0,elt0)
%BEMMESHCIRCLE  2D circular mesh.
%   [nod,elt,eltID] = BEMMESHCIRCLE(x0,R,nElt,typ,typID,orient,nod0,elt0) meshes 
%   a circle defined by its center and radius with equal-sized 2D boundary 
%   elements of the specified element type.
%   
%   x0       Vector [x0 0 z0] with coordinates of the circle center.
%   R        Circle radius.
%   nElt     Number of elements.
%   typ      Element types.
%   typID    Element type ID.
%   orient   Orientation of the element normal, either 1 (default) for outward 
%            element normal or -1 for inward element normal.
%   nod0     Node array to which new nodes are added. Defaults [].
%   elt0     Element array to which new elements are added. Defaults [].
%   nod      Node array.
%   elt      Element array. 
%   eltID    Element ID's of the meshed elements.

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

% INPUT ARGUMENT PROCESSING
if nargin<6, orient=1; end
if nargin<7, nod0=zeros(0,4); end
if nargin<8, elt0=zeros(0,2); end
if isempty(nod0), nod0=zeros(0,4); end
if isempty(elt0), elt0=zeros(0,2); end
if ~((orient==1)||(orient==-1)), error('Input argument ''orient'' should be either 1 or -1.'); end

[Parent,nEltNod,nCol,TypeN,TypeM,NodDef,EltDim]=bemeltdef(typID,typ);
if ~(EltDim==1), error('The element type specified is not 2D.'); end

nNod=nElt*(nEltNod-1);
theta=2*pi*linspace(0,-orient,nNod+1);
theta=theta(1:nNod).';
nod=[[1:nNod].'  R*cos(theta)+x0(1) zeros(nNod,1)+x0(2) R*sin(theta)+x0(3)];

elt=zeros(nElt,nEltNod+2);
elt(:,1)=1:nElt;
elt(:,2)=typID;
for iEltNod=1:nEltNod
  elt(:,2+iEltNod)=iEltNod:(nEltNod-1):(nNod+1)-(nEltNod-iEltNod);
end
elt(end)=elt(1,3); % close loop

[nod1,elt1,dum,eltID]=bemmeshcat(nod0,elt0,nod,elt);

