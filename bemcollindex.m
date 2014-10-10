function ind=bemcollindex(nod,elt,typ,eltSel)
%
%BEMCOLLINDEX   Collocation point indices.
% 
%  ind=BEMCOLLINDEX(nod,elt,typ,eltsel) returns the indices of the collocation
%  points of the boundary element mesh defined by the nod, elt and typ arrays 
%  corresponding to the elements eltsel. 
%
%   nod      Node array.
%   elt      Element array. 
%   typ      Element types.
%   eltsel   Element ID's of elements for which corresponding collocation
%            points are returned.
%   ind      Collocation point indices.
%
bemfunlicense('VerifyOnce');

% LOOKUP INDICES OF NEW COLLOCATION POINTS
[col,colTyp,colID]=bemcollpoints(nod,elt,typ);
nCol=size(col,1);

% FIND NODES OF ELTSEL
nEltSel=numel(eltSel);
nodSel=[];
for iElt=1:nEltSel
  [dum,eltind]=id2prop(elt,eltSel(iElt));
  [parent,nEltNod]=bemeltdef(elt(eltind,2),typ);
  nodSel=[nodSel;elt(eltind,3:nEltNod+2)];
end
nodSel=unique(nodSel);

% GO TROUGH COLLOCATION POINTS
ind=[];
for iCol=1:nCol
  switch colTyp(iCol)
  case 1 % element collocation
    k=find(colID(iCol)==eltSel,1);
    if ~isempty(k)
      ind=[ind;iCol];
    end
  case 2 % nodal collocation
    k=find(colID(iCol)==nodSel,1);
    if ~isempty(k)
      ind=[ind;iCol];
    end
  end
end


%===============================================================================
function [prop,propind]=id2prop(item,itemID)
%ID2PROP item properties from ID.
%   PROP=ID2PROP(ITEM,ITEMID)

nID=numel(itemID);
propind=[];
for iID=1:nID
  ind=find(item(:,1)==itemID(iID),1);
  propind=[propind,ind];
end
prop=item(propind,2:end);