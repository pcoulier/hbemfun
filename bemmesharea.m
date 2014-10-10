function [nod1,elt1,eltID] = bemmesharea(corner,n1,n2,typ,typID,nod0,elt0)
%BEMMESHAREA   Quadrilateral area meshing.
%
%   [nod,elt,eltID] = BEMMESHAREA(corner,nDiv1,nDiv2,typ,typID,nod0,elt0) 
%   meshes a quadrilateral area with the elements of a specified type.
%
%   corner Coordinates of 4 corner nodes of the quadrilateral (4 * 3).
%   nDiv1  Number of elements on the line from the first to the
%          second node.
%   nDiv2  Number of elements on the line from the second to the
%          third node.
%   typ    Element types.
%   typID  element type ID.
%   nod0     Node array to which new nodes are added. Defaults [].
%   elt0     Element array to which new elements are added. Defaults [].
%   nod    Nodes.
%   elt    Elements.
%   eltID  Element ID's of the meshed elements.

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

% INPUT ARGUMENT PROCESSING
if ~(size(corner,1)==4 & size(corner,2)==3) 
 error('Input argument ''corner'' should have dimension (4,3)')
end

if nargin<6, nod0=zeros(0,4); end
if nargin<7, elt0=zeros(0,2); end
if isempty(nod0), nod0=zeros(0,4); end
if isempty(elt0), elt0=zeros(0,2); end


% CHECK ELEMENT TYPE
typInd=[];
for iTyp=1:size(typ,1)
  if typ{iTyp,1}==typID, typInd=iTyp; end
end
if isempty(typInd), error('type ID not found in input argument ''typ'''); end
eltType = typ{typInd,2};

% PATTERN IN NATURAL COORDINATE SYSTEM
if strcmpi(eltType,'tria3')
  cellnod=[1 0 0;
           2 1 0;
           3 1 1;
           4 0 1];
  cellelt=[1 2 4;
           2 3 4];
elseif strcmpi(eltType,'quad4')
  cellnod=[1 0 0;
           2 1 0;
           3 1 1;
           4 0 1];
  cellelt=[1 2 3 4];
elseif strcmpi(eltType,'tria6')
  cellnod=[1 0   0;
           2 0.5 0;
           3 1   0;
           4 0   0.5;
           5 0.5 0.5;
           6 1   0.5;
           7 0   1;
           8 0.5 1
           9 1   1];
  cellelt=[1 3 7 2 5 4;
           3 9 7 6 8 5];
elseif strcmpi(eltType,'quad8')
  cellnod=[1 0   0;
           2 1   0;
           3 1   1;
           4 0   1;
           5 0.5 0;
           6 1   0.5;
           7 0.5 1;
           8 0   0.5];
  cellelt=[1 2 3 4 5 6 7 8];
elseif strcmpi(eltType,'quad9')
  cellnod=[1   0.0   0.0;
           2   1.0   0.0;
           3   1.0   1.0;
           4   0.0   1.0;
           5   0.5   0.0;
           6   1.0   0.5;
           7   0.5   1.0;
           8   0.0   0.5;
           9   0.5   0.5];
  cellelt=[1 2 3 4 5 6 7 8 9];
else
  error('Element type unsupported')
end

ncellelt=size(cellelt,1);
ncellnod=size(cellnod,1);
nelt=ncellelt*n1*n2;
neltnod=size(cellelt,2);
locnode=zeros(n1*n2*ncellnod,3);
locelt=zeros(nelt,neltnod);
for idiv1 = 1:n1
  for idiv2 = 1:n2
    icell=n2*(idiv1-1)+idiv2;
    locnode(ncellnod*(icell-1)+1:ncellnod*icell,1)=cellnod(:,1)+ncellnod*(icell-1);
    locnode(ncellnod*(icell-1)+1:ncellnod*icell,2)=(cellnod(:,2)+idiv1-1)/n1;
    locnode(ncellnod*(icell-1)+1:ncellnod*icell,3)=(cellnod(:,3)+idiv2-1)/n2;
    locelt(ncellelt*(icell-1)+1:ncellelt*icell,:)=cellelt+(icell-1)*ncellnod;
  end
end

% REMOVE COINCIDENT NODES
nlocnode =size(locnode,1);
coinc=zeros(nlocnode,2);
inod=0;
for inod1=1:nlocnode
  if coinc(inod1)==0
    inod=inod+1;
    for inod2=inod1:nlocnode
      if norm(locnode(inod1,2:end)-locnode(inod2,2:end))<1e-10
        coinc(inod2,1)=locnode(inod1,1);
        coinc(inod2,2)=inod;
      end
    end
  end
end
for ielt=1:nelt
  for inod=1:neltnod
    locelt(ielt,inod)=coinc(locelt(ielt,inod),2);
  end
end
[util,utilind]=unique(coinc(:,1));
locnode=[coinc(utilind,2) locnode(utilind,2:end)];

xi  = locnode(:,2);
eta = locnode(:,3);
N = [(1-xi).*(1-eta) xi.*(1-eta)  xi.*eta  (1-xi).*eta];

globnode = N * corner;
nnod=size(N,1);
nodeid=(1:nnod)';
nod =[nodeid globnode];

eltutil = locelt;
for inod = 1:neltnod
  for ielt =1:nelt
    eltutil(ielt,inod) = nodeid(eltutil(ielt,inod));
  end
end

elt=[(1:nelt)' repmat(typID,nelt,1) eltutil];

[nod1,elt1,dum,eltID]=bemmeshcat(nod0,elt0,nod,elt);
