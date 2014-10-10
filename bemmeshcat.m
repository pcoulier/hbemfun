function [nod,elt,nod2ID,elt2ID]=bemmeshcat(varargin)
%BEMMESHCAT concatenate boundary element meshes.
%
%  [nod,elt]=BEMMESHCAT(nod1,elt1,nod2,elt2) concatenates the boundary 
%  element meshes 1 and 2. Node and element renumbering is performed to
%  ensure unique node and element numbers.
%  
%  [nod,elt,nodID2,eltID2]=BEMMESHCAT(nod1,elt1,nod2,elt2) returns the node 
%  and element ID's of nod2 and elt2 in the resulting node and element arrays
%  nod and elt.
%
%  [nod,elt]=BEMMESHCAT(nod1,elt1,nod2,elt2,...,nodN,eltN)
%  concatenates multiple meshes at once.

bemfunlicense('VerifyOnce');

if nargin<4, error('minimum 4 input arguments required'); end
if rem(nargin,2), error('number of input arguments should be even'); end

if nargin>4
   if (nargout>2), error('Too many output arguments.'); end
  [nod3,elt3]=bemmeshcat(varargin{1:end-2});
  [nod,elt]=bemmeshcat(nod3,elt3,varargin{end-1:end});
else
  nod1=varargin{1};
  elt1=varargin{2};
  nod2=varargin{3};
  elt2=varargin{4};
  nNod1=size(nod1,1);
  nNod2=size(nod2,1);
  nElt1=size(elt1,1);
  nElt2=size(elt2,1);
  
  % CHANGE ID'S TO AVAILABLE ID'S
  nEltColumn1=size(elt1,2);
  nEltColumn2=size(elt2,2);
  colDiff=nEltColumn2-nEltColumn1;
  nod2ID=availID(nod1(:,1),nNod2);
  elt2ID=availID(elt1(:,1),nElt2);
  nod2(:,1)=nod2ID;
  elt2(:,1)=elt2ID;
  for iElt=1:nElt2
   for jElt=3:nEltColumn2
     elt2(iElt,jElt)= nod2ID(elt2(iElt,jElt));
   end
  end
  
  % CONCATENATE NODE AND ELT ARRAYS
  nod=[nod1;nod2];
  elt=[elt1 zeros(nElt1,max(0,colDiff));
       elt2 zeros(nElt2,max(0,-colDiff))];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newID=availID(ID,nNewID)
% Get the available ID's from a list of ID's in use

nID=numel(ID);
newID=zeros(nNewID,1);

iNewID=0;
iID=0;
while iNewID<nNewID
  iID=iID+1;
  if ~ismember(iID,ID)
    iNewID=iNewID+1;
    newID(iNewID)=iID;
  end
end