function [nod,elt,eltID]=bemmeshdisk(x0,normal,R,nElt,typ,typID,nod0,elt0)
%BEMMESHDISK  3D disk mesh.
%   [nod,elt,eltID] = BEMMESHDISK(x0,normal,R,nElt,typ,typID,nod0,elt0) meshes
%   a disk defined by its center and radius with 3D boundary elements 
%   of the specified element type.
%   
%   x0       Vector [x0 y0 z0] with coordinates of the center of the disk.
%   R        Radius.
%   normal   Vector normal to disk [nx ny nz].
%   nElt     Desired number of elements. The number of elements equals 
%            nElt = 12 * n^2 where n is an integer.
%   typ      Element types.
%   typID    Element type ID.
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
if nargin<7, nod0=zeros(0,4); end
if nargin<8, elt0=zeros(0,2); end
if isempty(nod0), nod0=zeros(0,4); end
if isempty(elt0), elt0=zeros(0,2); end

[Parent,nEltNod,nCol,TypeN,TypeM,NodDef,EltDim]=bemeltdef(typID,typ);
if ~(EltDim==2), error('The element type specified is not 3D.'); end

nDiv=round(sqrt(nElt/12));
nDiv=max(nDiv,1);

[nod,elt]=bemmesharea([0.5,pi/4,0;0,pi/4,0;0,-pi/4,0;0.5,-pi/4,0],nDiv,2*nDiv,typ,typID);
nNod=size(nod,1);
nElt=size(elt,1);
theta=atan(4/pi*nod(:,3)).*(1-1/0.5*nod(:,2))  + nod(:,3).*nod(:,2)/0.5;
radius=(0.500/0.5*nod(:,2)+1./cos(theta).*(0.500-0.500/0.5*nod(:,2))) + nod(:,2);

nod1=nod;
nod1(:,2)=radius.*cos(theta);
nod1(:,3)=radius.*sin(theta);
nod1(:,4)=0;
elt1=elt;

nod2(:,1)=nod(:,1)+1*nNod;
nod2(:,2)=radius.*cos(theta+1*pi/2);
nod2(:,3)=radius.*sin(theta+1*pi/2);
nod2(:,4)=0;
elt2=elt;
elt2(:,1)=elt(:,1)+1*nElt;
elt2(:,3:end)=elt(:,3:end)+1*nNod;

nod3(:,1)=nod(:,1)+2*nNod;
nod3(:,2)=radius.*cos(theta+2*pi/2);
nod3(:,3)=radius.*sin(theta+2*pi/2);
nod3(:,4)=0;
elt3=elt;
elt3(:,1)=elt(:,1)+2*nElt;
elt3(:,3:end)=elt(:,3:end)+2*nNod;

nod4(:,1)=nod(:,1)+3*nNod;
nod4(:,2)=radius.*cos(theta+3*pi/2);
nod4(:,3)=radius.*sin(theta+3*pi/2);
nod4(:,4)=0;
elt4=elt;
elt4(:,1)=elt(:,1)+3*nElt;
elt4(:,3:end)=elt(:,3:end)+3*nNod;

[nod,elt]=bemmesharea([0.500,0.500,0;-0.500,0.500,0;-0.500,-0.500,0;0.500,-0.500,0],2*nDiv,2*nDiv,typ,typID,[nod1;nod2;nod3;nod4],[elt1;elt2;elt3;elt4]);
nNod=size(nod,1);
nElt=size(elt,1);

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

% CORRECT FOR RADIUS AND CENTER
nod(:,2:3)=R*nod(:,2:3);
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

