function elt=bemeltreverse(elt0,typ)
%BEMELTREVERSE   Reverse boundary element normals.
%
%   elt = BEMELTREVERSE(elt0,typ) reverses the normals of all boundary elements 
%   defined in elt0.
%
%   elt0 Elements.
%   typ  Element types.
%   elt  Reversed elements.

bemfunlicense('VerifyOnce');

nElt=size(elt0,1);
elt=elt0;
for iElt=1:nElt
  typID=elt0(iElt,2);
  typInd=[];
  for iTyp=1:size(typ,1)
    if typ{iTyp,1}==typID, typInd=iTyp; end
  end
  if isempty(typInd), error('type ID not found in input argument ''typ'''); end
  eltType = typ{typInd,2};
  
  if strcmp(eltType,'tria3')
    ReverseInd=[1 3 2];
  elseif strcmp(eltType,'quad4')
    ReverseInd=[1 4 3 2];
  elseif strcmp(eltType,'tria6')
    ReverseInd=[1 3 2 6 5 4];
  elseif strcmp(eltType,'quad8')
    ReverseInd=[1 4 3 2 8 7 6 5];
  elseif strcmp(eltType,'quad9')
    ReverseInd=[1 4 3 2 8 7 6 5 9];
  elseif strcmp(eltType,'line2')
    ReverseInd=[2 1];
  elseif strcmp(eltType,'line3')
    ReverseInd=[3 2 1];
  elseif strcmp(eltType,'line4')
    ReverseInd=[4 3 2 1];
  else
    error(['Element type ' num2str(eltType) ' not supported for element reversing']);
  end
  nEltNod=length(ReverseInd);
  elt(iElt,3:2+nEltNod)=elt0(iElt,2+ReverseInd);
end
