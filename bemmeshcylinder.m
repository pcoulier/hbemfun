function [nod,elt,eltID]=bemmeshcylinder(x0,normal,R,H,nEltR,nEltH,typ,typID,orient,nod0,elt0)
%BEMMESHCYLINDER 3D cylinder mesh
%   [A,AIND] = BEMMESHCYLINDER(x0,normal,R,H,nEltR,nEltH,typ,typID,orient,nod0,elt0) 
%   meshes a cylindrical surface defined by a center node, its radius, height 
%   and direction.
%   
%   x0       Vector [x0 y0 z0] with coordinates of the origin of the cylinder.
%   normal   Vector parallel to cylinder axis [nx ny nz].
%   R        Cylinder radius.
%   H        Cylinder height.
%   nElt     Desired number of elements. The number of elements is rounded
%            to the nearest multiple of 12.
%   typ      Element types.
%   typID    Element type ID.
%   orient   Orientation of the element normal, either 1 (default) for outward 
%            element normal or -1 for inward element normal.
%   nod0     Node array to which new nodes are added. Defaults [].
%   elt0     Element array to which new elements are added. Defaults [].
%   nod      Node array.
%   elt      Element array. 
%   eltID    Element ID's of the meshed elements.

% Stijn Francois
% August 2008

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

% INPUT ARGUMENT PROCESSING
if nargin<9, orient=1; end
if nargin<10, nod0=zeros(0,4); end
if nargin<11, elt0=zeros(0,2); end
if isempty(nod0), nod0=zeros(0,4); end
if isempty(elt0), elt0=zeros(0,2); end

[Parent,nEltNod,nCol,TypeN,TypeM,NodDef,EltDim]=bemeltdef(typID,typ);
if ~(EltDim==2), error('The element type specified is not 3D.'); end

% CONSTRUCT UNIT CYLINDER
[nod,elt]=bemmesharea([0,0,0;2*pi,0,0;2*pi,0,H;0,0,H],nEltR,nEltH,typ,typID);
nNod=size(nod,1);
nElt=size(elt,1);
theta=nod(:,2);
nod(:,2)=R*cos(theta);
nod(:,3)=R*sin(theta);

% REMOVE COINCIDENT NODES
inod=0;
coinc=zeros(nNod,2);
for inod1=1:nNod
  if coinc(inod1)==0
    inod=inod+1;
    for inod2=inod1:nNod
      if norm(nod(inod1,2:end)-nod(inod2,2:end))<1e-10
        coinc(inod2,1)=nod(inod1,1);
        coinc(inod2,2)=inod;
      end
    end
  end
end

nEltNod=size(elt,2)-2;
for ielt=1:nElt
  for inod=3:nEltNod+2
    elt(ielt,inod)=coinc(elt(ielt,inod),2);
  end
end
[util,utilind]=unique(coinc(:,1));
nod=[coinc(utilind,2) nod(utilind,2:end)];
nNod=size(nod,1);

% CORRECT FOR CENTER
nod(:,2)=nod(:,2)+x0(1);
nod(:,3)=nod(:,3)+x0(2);
nod(:,4)=nod(:,4)+x0(3);

% ROTATE MESH
normal=normal/norm(normal);
util1=cross(normal,[1 0 0]);
if norm(util1)==0
   util1=cross(normal,[0 1 0]);
   util1=util1./norm(util1);
end
util2=cross(util1,normal);
util2=util2./norm(util2);
nod=bemmeshrotate(nod,x0,[util2;util1;normal].');

% CHANGE ELEMENT ORIENTATION
if orient==-1
  elt=bemeltreverse(elt,typ);
end

% CHANGE ID'S TO AVAILABLE ID'S
[nod,elt,dum,eltID]=bemmeshcat(nod0,elt0,nod,elt);


%-------------------------------------------------------------------------------
function [nod]=bemmeshrotate(nod0,x0,A)
%BEMMESHROTATE  Rotate finite element mesh.
%   [nod] = BEMMESHROTATE(nod0,x0,angles) rotates the nodes of a mesh.
%   
%   nod0     Unrotated nodes.
%   x0       Center of rotation [x0 y0 z0].
%   A        Direction cosines [n1;n2;n3].
%   nod      Rotated nodes.

nNod=size(nod0,1);
nod=nod0;
for iNod=1:nNod
  vec0=nod0(iNod,2:4);
  vec=A*(vec0(:)-x0(:))+x0(:);
  nod(iNod,2:4)=vec;
end