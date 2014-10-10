function t=bempressure(nod,elt,typ,p,nColDof)
%
% BEMPRESSURE Apply tractions normal to elements.
%
% t=BEMPRESSURE(nod,elt,typ,p) constructs the boundary element traction
% vector corresponding to a pressure at the boundary element collocation points.
%
%   nod   Nodes (nNod * 4). Each row has the layout [nodID x y z] where
%         nodID is the node number and x, y, and z are the nodal
%         coordinates.
%   elt   Elements (nElt * nColumn). Each row has the layout
%         [eltID typID n1 n2 n3 ... ] where eltID is the element number,
%         typID is the element type number and n1, n2, n3, ... are the node
%         numbers representing the nodal connectivity of the element.
%   typ   Element type definitions. Cell array with the layout
%         {{typID type keyOpts} ... } where typID is the element type number,
%         type is the element type (string) and keyOpts is a cell array of
%         strings with key options.
%   p     Pressure (1 * 1).
%   t     Boundary element tractions (nDof * 1).


bemfunlicense('VerifyOnce');

if nargin<5
  dim=bemdimension(elt,typ);
  if dim==2, nColDof=2; end
  if dim==3, nColDof=3; end
end
if ~(nColDof==2 | nColDof==3)
  error('Number of degrees of freedom per collocation point should be 2 or 3.');
end

[col,colTyp,colID]=bemcollpoints(nod,elt,typ);

nCol=size(col,1);
t=zeros([nColDof*nCol 1]);

for iCol=1:nCol
  if colTyp(iCol)==1  % Element collocation
    
    ielt=id2sel(elt,colID(iCol));
    [parent,nEltNod,nEltCol,TypeN,TypeM,NodXi]=bemeltdef(ielt(2),typ);
    
    if parent==0         % 1D element
      colXi=[0 0];
    elseif parent==1     % 2D trianglular element
      colXi=[1/3 1/3];
    elseif parent==2     % 2D quadrilateral element
      colXi=[0 0];
    end
    normal = bemnormal(nod,ielt,typ,colXi);
    if nColDof==2
      t(2*(iCol-1)+1)=-p*normal(1);
      t(2*(iCol-1)+2)=-p*normal(3);
    elseif nColDof==3
      t(3*(iCol-1)+1)=-p*normal(1);
      t(3*(iCol-1)+2)=-p*normal(2);
      t(3*(iCol-1)+3)=-p*normal(3);
    end
  else % Nodal collocation
    iNod=colID(iCol);
    [iEltPos,iEltNod]=find(elt(:,3:end)==iNod);
    nElt=numel(iEltPos);
    for iPos=1:nElt
      ielt=elt(iEltPos(iPos),:);
      [parent,nEltNod,nEltCol,TypeN,TypeM,NodXi]=bemeltdef(ielt(2),typ);
      normal = bemnormal(nod,ielt,typ,NodXi(iEltNod(iPos),:));
      if nColDof==2
        t(2*(iCol-1)+1)=t(2*(iCol-1)+1)-p*normal(1)/nElt;
        t(2*(iCol-1)+2)=t(2*(iCol-1)+2)-p*normal(3)/nElt;
      elseif nColDof==3
        t(3*(iCol-1)+1)=t(3*(iCol-1)+1)-p*normal(1)/nElt;
        t(3*(iCol-1)+2)=t(3*(iCol-1)+2)-p*normal(2)/nElt;
        t(3*(iCol-1)+3)=t(3*(iCol-1)+3)-p*normal(3)/nElt;
      end
    end
  end
end

%===============================================================================
function [prop,propind]=id2sel(item,itemID)
%ID2PROP item properties from ID.
%   PROP=ID2PROP(ITEM,ITEMID)

nID=numel(itemID);
propind=[];
for iID=1:nID
  ind=find(item(:,1)==itemID(iID),1);
  propind=[propind,ind];
end
prop=item(propind,1:end);
