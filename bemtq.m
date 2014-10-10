function Tq=bemtq(nod,elt,typ,probDim)
%BEMTQ   Boundary element stress transfer matrix.
%
%   Tq = BEMTQ(nod,elt,typ,probDim) computes the boundary element stress 
%   transfer matrix Tq defined as the integral:
%  
%           /
%   Tq_ij = |  M_i * N_j  dS
%           /
%         Gamma
%  
%   of the product  M_i * N_j of the boundary element interpolation function
%   M_i and the finite element interpolation function N_j. In the case of a
%   conforming boundary element-finite element coupling, Tq relates the
%   boundary element tractions t and the finite element
%   forces Q as: 
%   
%       Q=Tq*t
%   
%   This relationship is derived from the principle of virtual work along the 
%   interface: the virtual work performed by the load vector Q and the interface
%   traction should be equal under any virtual displacement field. 
%   
%   The traction interpolation function is constant in the case of a centroid
%   collocated boundary element. In this case, the finite element interpolation
%   function is derived from the equivalent nodal collocated boundary element.
%   Hence, the function N_j is at least linear.
%   
%   nod  Nodes.
%   elt  Elements.
%   typ  Element types.
%   Tq   Sparse stress transfer matrix (nDof * nDof).


%  UNSUPPORTED SYNTAX
%   Tq = BEMTQ(NOD,ELT,PHI) computes the boundary stress transfer matrix Tq with
%   respect to the mode PHI. In this case, the matrix TQ relates the boundary
%   element tractions t to the modal force fm as fm=Tq*t.

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

if (nargin<4), probDim=3; end;
[col,colType,colID]=bemcollpoints(nod,elt,typ);
nCol=size(col,1);
nElt=size(elt,1);
nNod=size(nod,1);

nDof=probDim*nCol;
Tq=sparse(nDof,nDof);
  
for iElt=1:nElt
  thisElt=elt(iElt,:);
  [parent,nEltNod,nEltCol,TypeN,TypeM,NodXi]=bemeltdef(thisElt(2),typ);
  if (nEltCol==1) % Centroid collocated element
    EltCollocation=1;
    nEltCol=1;
    tbasis=repmat(eye(probDim),nEltNod,1);
  else            % Nodal collocated element
    EltCollocation=0;
    nEltCol=nEltNod;
    tbasis=eye(probDim*nEltCol);
  end
  ubasis=eye(probDim*nEltNod);
  
  [thisEltCol,thisEltColTyp,thisEltColID]=bemcollpoints(nod,thisElt,typ);
  TqElt=bemint(nod,thisElt,typ,tbasis,ubasis);
  
  
% LOOKUP COLLOCATION POINTS INDICES
  % Element collocation points
   if EltCollocation
     for iCol=1:nCol
       if (colType(iCol)==1 & colID(iCol)==thisEltColID)
         eltColIndex=iCol;
       end
     end
   % Nodal collocation points
   else
     eltColIndex=zeros(1,nEltNod);
     for iEltNod=1:nEltNod
       for iCol=1:nCol
         if (colType(iCol)==2 & colID(iCol)==thisEltColID(iEltNod))
           eltColIndex(iEltNod)=iCol;
         end
       end
     end
   end
   
   ind=zeros(1,probDim*nEltCol);
   for iDim=1:probDim
     ind(iDim:probDim:probDim*nEltCol)=probDim*(eltColIndex-1)+iDim;
   end
   Tq(ind,ind)=Tq(ind,ind)+TqElt;
end
