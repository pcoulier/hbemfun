function [nod,elt,eltID]=bemmeshsphere(x0,R,nElt,typ,typID,orient,nod0,elt0)
%BEMMESHSPHERE  3D sphere mesh.
%   [nod,elt,eltID] = BEMMESHSPHERE(x0,R,nElt,typ,typID,orient,nod0,elt0) meshes
%   a sphere defined by its center and radius with 3D boundary elements 
%   of the specified element type.
%   
%   x0       Vector [x0 y0 z0] with coordinates of the center of the sphere.
%   R        Radius.
%   nElt     Desired number of elements.
%   typ      Element types.
%   typID    Element type ID.
%   orient   Orientation of the element normal, either 1 (default) for outward 
%            element normal or -1 for inward element normal.
%   nod0     Node array to which new nodes are added. Defaults [].
%   elt0     Element array to which new elements are added. Defaults [].
%   nod      Node array.
%   elt      Element array. 
%   eltID  Element ID's of the meshed elements.

% Stijn Francois
% August 2008

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
if ~(EltDim==2), error('The element type specified is not 2D.'); end

nDiv=floor(sqrt(nElt/6));
nDiv=max(nDiv,1);

%% MESH A CUBE
[nod1,elt1]=bemmesharea([1,1,1;-1,1,1;-1,-1,1;1,-1,1],nDiv,nDiv,typ,typID);
nNod1=size(nod1,1);
nElt1=size(elt1,1);

nod2=nod1;
nod2(:,1)=nod1(:,1)+nNod1;
nod2(:,4)=-1;
elt2=elt1;
elt2(:,1)=elt1(:,1)+nElt1;
elt2(:,3:end)=elt1(:,3:end)+nNod1;
elt2=bemeltreverse(elt2,typ);

nod3=nod1;
nod3(:,1)=nod1(:,1)+2*nNod1;
nod3(:,2)=nod1(:,4);
nod3(:,4)=nod1(:,2);
elt3=elt1;
elt3(:,1)=elt1(:,1)+2*nElt1;
elt3(:,3:end)=elt1(:,3:end)+2*nNod1;
elt3=bemeltreverse(elt3,typ);

nod4=nod1;
nod4(:,1)=nod1(:,1)+3*nNod1;
nod4(:,2)=-nod1(:,4);
nod4(:,4)=nod1(:,2);
elt4=elt1;
elt4(:,1)=elt1(:,1)+3*nElt1;
elt4(:,3:end)=elt1(:,3:end)+3*nNod1;

nod5=nod1;
nod5(:,1)=nod1(:,1)+4*nNod1;
nod5(:,3)=nod1(:,4);
nod5(:,4)=nod1(:,3);
elt5=elt1;
elt5(:,1)=elt1(:,1)+4*nElt1;
elt5(:,3:end)=elt1(:,3:end)+4*nNod1;
elt5=bemeltreverse(elt5,typ);

nod6=nod1;
nod6(:,1)=nod1(:,1)+5*nNod1;
nod6(:,3)=-nod1(:,4);
nod6(:,4)=nod1(:,3);
elt6=elt1;
elt6(:,1)=elt1(:,1)+5*nElt1;
elt6(:,3:end)=elt1(:,3:end)+5*nNod1;

nod=[nod1;nod2;nod3;nod4;nod5;nod6];
elt=[elt1;elt2;elt3;elt4;elt5;elt6];
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
r=sqrt(sum(nod(:,2:4).^2,2));
nod(:,2:4)=R*nod(:,2:4).*[1./r 1./r 1./r];
nod(:,2)=nod(:,2)+x0(1);
nod(:,3)=nod(:,3)+x0(2);
nod(:,4)=nod(:,4)+x0(3);

% REVERSE ELEMENT ORIENTATION IF DESIRED
if (orient==-1)
  elt=bemeltreverse(elt,typ);
end

% CHANGE ID'S TO AVAILABLE ID'S
[nod,elt,dum,eltID]=bemmeshcat(nod0,elt0,nod,elt);