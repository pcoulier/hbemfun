function [nod1,elt1,eltID] = bemmeshline(x1,x2,nElt,typ,typID,nod0,elt0)
%BEMMESHLINE   2D line meshing.
%
%   [nod,elt,eltID] = BEMMESHLINE(x1,x2,nElt,typ,typID,nod0,elt0) meshes a line
%   defined by its end points with equal-sized 2D boundary elements of the
%   specified element type.
%
%   x1       Vector [x0 0 z0] with coordinates of the starting point of the line.
%   x2       Vector [x1 0 z1] with coordinates of the end point of the line.
%   nElt     Number of elements on the line.
%   typ      Element types.
%   typID    Element type ID.
%   nod0     Node array to which new nodes are added. Defaults [].
%   elt0     Element array to which new elements are added. Defaults [].
%   nod      Node array.
%   elt      Element array. 
%   eltID    Element ID's of the meshed elements.

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

% INPUT ARGUMENT PROCESSING
if nargin<6, nod0=zeros(0,4); end
if nargin<7, elt0=zeros(0,2); end
if isempty(nod0), nod0=zeros(0,4); end
if isempty(elt0), elt0=zeros(0,2); end

[Parent,nEltNod,nCol,TypeN,TypeM,NodDef,EltDim]=bemeltdef(typID,typ);
if ~(EltDim==1), error('The element type specified is not 2D.'); end

nNod=(nEltNod-1)*nElt+1;
Nod=zeros(nNod,4);
Nod(:,1)=1:nNod;
Nod(:,2)=linspace(x1(1),x2(1),nNod);
Nod(:,4)=linspace(x1(3),x2(3),nNod);

Elt=zeros(nElt,nEltNod+2);
Elt(:,1)=1:nElt;
Elt(:,2)=typID;
for iEltNod=1:nEltNod
  Elt(:,2+iEltNod)=iEltNod:(nEltNod-1):nNod-(nEltNod-iEltNod);
end

[nod1,elt1,dum,eltID]=bemmeshcat(nod0,elt0,Nod,Elt);