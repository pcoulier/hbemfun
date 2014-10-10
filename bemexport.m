function bemexport(nod,elt,typ,exportformat,varargin)
%BEMEXPORT   Boundary element mesh exporting.
%
%   BEMEXPORT(nod,elt,typ,'miss',file) exports a mesh for the use in MISS.
%
%   nod   Nodes.
%   elt   Elements.
%   typ   Element types.
%   file  Resulting file for the use in MISS. Commonly named "mail.txt".

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

switch lower(exportformat)
   case 'miss'
      exportfile=varargin{1};
      nNod = size(nod,1);
      nElt = size(elt,1);
      nEltCol=size(elt,2);

      % Check element type
      typID=unique(elt(:,2));
      nTypID=length(typID);
      for iTypID=1:nTypID
        typInd=[];
        for iTyp=1:size(typ,1)
          if typ{iTyp,1}==typID, typInd=iTyp; end
        end    
        if isempty(typInd), error('type ID not found in input argument ''typ'''); end
        eltType = typ{typInd,2};
        if ~(strcmpi(eltType,'quad4')), error('export routine only supports quad4 elements'); end
      end

      % Change nodeID to nodeindex in the element-array
      [nodcoord,nodind]=id2prop(nod,elt(:,3:nEltCol));
      elt(:,3:nEltCol)=reshape(nodind,nElt,nEltCol-2);

      fid = fopen(exportfile,'w');
      fprintf(fid,'MAILLAGE BEMFUNMESH: %g noeuds %g elements Q4\n %g %g\n(3E21.13)\n', ...
              nNod,nElt,nNod,nElt);
      fprintf(fid,'%20.13f %20.13f %20.13f\n',nod(:,2:4).');
      fprintf(fid,'%d 0 %d 0 %d 0 %d 0 GR 1\n',elt(:,3:6).');
      fprintf(fid,'*');
      fclose(fid);
   otherwise
      error('Export format unknown');
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