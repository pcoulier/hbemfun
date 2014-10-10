function Tu=bemtu(nod,elt,typ,probDim)
%BEMTU   Boundary element displacement transfer matrix.
%
%   Tu = BEMTU(nod,elt,typ) computes the boundary element displacment transfer
%   matrix defined as:
%
%   u=Tu*unod                                          
%
%   where u are the boundary element degrees of freedom and unod are the
%   nodal degrees of freedom. u and unod are equal for a nodal collocated
%   boundary element formulation, in which case Tu is a unity matrix.
%   
%   In the case of a centroid collocated boundary element formulation,
%   the elements of the matrix Tu are derived from the boundary element
%   shape functions.
%   
%   nod  Nodes.
%   elt  Elements.
%   typ  Element types.
%   Tu   Sparse displacement transfer matrix.

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

if (nargin<4), probDim=3; end;
[col,colType,colID]=bemcollpoints(nod,elt,typ);
nCol=size(col,1);
nNod=size(nod,1);
nElt=size(elt,1);

nColDof=probDim*nCol;
nNodDof=probDim*nNod;

Tu=sparse(nColDof,nNodDof);

for iElt=1:nElt
  thisElt=elt(iElt,:);
  [parent,nEltNod,nEltCol,TypeN,TypeM,NodXi]=bemeltdef(thisElt(2),typ);
  
  if nEltCol==1
    if parent==0         % 1D element
      colXi=[0 0];
    elseif parent==1     % 2D trianglular element
      colXi=[1/3 1/3];
    elseif parent==2     % 2D quadrilateral element
      colXi=[0 0];
    end
  else
    colXi=NodXi;
  end
  N=bemshape(TypeN,colXi);
  
  Tuelt=zeros(probDim*nEltCol,probDim*nEltNod);
  for iDim=1:probDim
    Tuelt(iDim:probDim:probDim*nEltCol,iDim:probDim:probDim*nEltNod)=N;
  end
  
  % LOOKUP NODES FOR THIS ELEMENT
  eltNodes=zeros(1,nEltNod);
  for iEltNod=1:nEltNod
    for iNod=1:nNod
      if (nod(iNod,1)==thisElt(2+iEltNod));
       eltNodes(iEltNod)=iNod;
      end
    end
  end
  
  % LOOKUP COLLOCATION POINTS INDICES
  % Element collocation points
   if nEltCol==1
     for iCol=1:nCol
       if (colType(iCol)==1 & colID(iCol)==elt(iElt,1))
         eltColIndex=iCol;
       end
     end
   %  Nodal collocation points
   else
     eltColIndex=zeros(1,nEltNod);
     for iNod=1:nEltNod
       for iCol=1:nCol
         if (colType(iCol)==2 & colID(iCol)==elt(iElt,3+iNod))
           eltColIndex(iNod)=iCol;
         end
       end
     end
   end  
   
   ind0=zeros(1,probDim*nEltNod);
   for iDim=1:probDim
     ind0(iDim:probDim:probDim*nEltNod)=probDim*(eltNodes-1)+iDim;
   end
   
   ind1=zeros(1,probDim*nEltCol);
   for iDim=1:probDim
     ind1(iDim:probDim:probDim*nEltCol)=probDim*(eltColIndex-1)+iDim;
   end
   
   Tu(ind1,ind0)=Tuelt;
end
