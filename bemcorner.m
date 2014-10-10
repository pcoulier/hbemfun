function [Tc,ind1,ind2]=bemcorner(nod,elt,typ,E,nu,nColDof,CoincEps,ky)
%BEMCORNER Boundary element constraint equations at corners.
%
%   [Tc,ind1,ind2] = BEMCORNER(nod,elt,typ,E,nu,nColDof,tol,ky)
%   generates constraint equations at corners and edges of a boundary element
%   mesh. A distinction is made between traction degrees of freedom t1=t(ind1)
%   not related to a corner and traction degrees of freedom t2=t(ind2) related
%   to corners or edges. The traction degrees of freedom t2 are related to the
%   boundary element displacments as:
%
%          t2=Tc*u
%
%   By means of this equation, Dirichlet boundary conditions are converted
%   to Neumann conditions on corners and edges.
%
%   nod      Nodes (nNod * 4). Each row has the layout [nodID x y z] where
%            nodID is the node number and x, y, and z are the nodal
%            coordinates.
%   elt      Elements (nElt * nColumn). Each row has the layout
%            [eltID typID n1 n2 n3 ... ] where eltID is the element number,
%            typID is the element type number and n1, n2, n3, ... are the node
%            numbers representing the nodal connectivity of the element.
%   typ      Element type definitions. Cell array with the layout
%            {{typID type keyOpts} ... } where typID is the element type number,
%            type is the element type (string) and keyOpts is a cell array of
%            strings with key options.
%   E        Young's modulus (1 * 1).
%   nu       Poisson coefficient (1 * 1).
%   nColDof  Number of DOFs per collocation point. 1 for 2D in-plane problems,
%            2 for 2D out-of-plane problems, 3 for 2.5D problems and 
%            3 for 3D problems. Defaults 3 for 3D problems and 2 for 2D problems.
%   tol      Tolerance for node coincidence (1 * 1). If two nodes are located less
%            than the tolerance from another, they are considered to be at a corner.
%            Defaults 1e-6.
%   ky       Wavenumber for 2.5D problems (2D mesh, nColDof=3).
%   Tc       Constraint equation coefficients for displacements (nDof2 * nDof)
%   ind1     Indices of traction degrees of freedom not on a corner (nDof1 * 1).
%   ind2     Indices of traction degrees of freedom on a corner (nDof2 * 1).
%

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

dim=bemdimension(elt,typ);
if nargin<6  nColDof=dim; end
if nargin<7, CoincEps=1e-6; end
if nargin<8, ky=0; end

% MATERIAL PROPERTIES
lambda=E*nu/(1+nu)/(1-2*nu);
mu=E/2/(1+nu);

% FIND SLAVE NODES
[col,colTyp,ID] = bemcollpoints(nod,elt,typ);
ID(find(colTyp==1))=NaN; % REMOVE ELEMENT COLLOCATION POINTS FROM ID.

nCol=size(col,1);
slaves=zeros(nCol,1);  % master=0, slave=masterindex
for iCol=1:nCol
  for jCol=iCol+1:nCol
    if (norm(col(iCol,:)-col(jCol,:)) < CoincEps) && ~slaves(jCol)
      slaves(jCol)=iCol;
    end
  end
end
slaveInd=find(slaves);
nSlave=length(slaveInd);

% CORNER INDICES
cornerInd=unique([slaves(slaveInd) slaveInd]);
notcornerInd=setdiff([1:nCol].',cornerInd);
if nColDof==2
  ind1=sort([nColDof*notcornerInd-1; nColDof*notcornerInd]);
  ind2=sort([nColDof*cornerInd-1;    nColDof*cornerInd]);
elseif nColDof==3
  ind1=sort([nColDof*notcornerInd-2; nColDof*notcornerInd-1; nColDof*notcornerInd]);
  ind2=sort([nColDof*cornerInd-2;   nColDof*cornerInd-1;     nColDof*cornerInd]);
end
ind1=ind1(:);
ind2=ind2(:);

ID1=nod(cornerInd,:);
ID2=nod(notcornerInd,:);

nCornerNod=length(cornerInd);

Tc=zeros(nColDof*nCornerNod,nColDof*nCol);
for iSlave=1:nSlave
  masterInd=slaves(slaveInd(iSlave));

  % ELEMENT OF MASTER NODE
  [iElt,iEltNod]=find(elt(:,3:end)==ID(masterInd),1); % Take normal of first element found
  [mEltMap,mEltNod,mEltCol,mNShape,mMShape,mNodXi]=bemeltdef(elt(iElt,2),typ);
  mXi=mNodXi(iEltNod,:);
  n1 = bemnormal(nod,elt(iElt,:),typ,mXi);
  tan1 = bemtangent(nod,elt(iElt,:),typ,mXi);
  mInd=[];
  for iEltNod=1:mEltNod
    mInd(iEltNod)=find((ID==elt(iElt,2+iEltNod)));
  end
  mCorner=find(cornerInd==masterInd);
  if dim==2, mXi=mXi(1); end
  dN1 = bemshapederiv(mMShape,mXi);

  % ELEMENT OF SLAVE NODE
  [iElt,iEltNod]=find(elt(:,3:end)==ID(slaveInd(iSlave)),1); % Take normal of first element found
  [nEltMap,nEltNod,nEltCol,nNShape,nMShape,nNodXi]=bemeltdef(elt(iElt,2),typ);
  nXi=nNodXi(iEltNod,:);
  n2 = bemnormal(nod,elt(iElt,:),typ,nXi);
  tan2 = bemtangent(nod,elt(iElt,:),typ,nXi);
  nInd=[];
  for iEltNod=1:nEltNod
    nInd(iEltNod)=find((ID==elt(iElt,2+iEltNod)));
  end
  nCorner=find(cornerInd==slaveInd(iSlave));
  if dim==2, nXi=nXi(1); end
  dN2 = bemshapederiv(nMShape,nXi);

  if dim==2   % 2D ROUTINE
    dxdxi=[tan1(1) tan1(3);tan2(1) tan2(3)];
    dxidx=inv(dxdxi);
    if nColDof==1 % 2D OUT-OF-PLANE
      T=dxidx;
      
      B=zeros(2,nColDof*nCol);
      for iNod=1:mEltNod
        B(1,nColDof*(mInd(iNod)-1)+1)=dN1(iNod);
      end
      for iNod=1:nEltNod
        B(2,nColDof*(nInd(iNod)-1)+1)=dN2(iNod);
      end
      
      C=[mu 0         % Constitutive equation
         0  mu];
         
      sigma=C*T*B;
      
       % Tractions
      N1=[n1(1) n1(3)];
      N2=[n2(1) n2(3)];
      Tc(nColDof*(mCorner-1)+1,:)=N1*sigma;
      Tc(nColDof*(nCorner-1)+1,:)=N2*sigma;
      
    elseif nColDof==2 % 2D IN-PLANE
      T=zeros(4,4);
      T(1:2,1:2)=dxidx; % Transformation matrix from xi to x.
      T(3:4,3:4)=dxidx;
      A=[1 0 0 0; 0 0 0 1; 0 1 1 0]; % Gradient to strain

      B=zeros(4,nColDof*nCol);
      for iNod=1:mEltNod
        B(1,nColDof*(mInd(iNod)-1)+1)=dN1(iNod);
        B(3,nColDof*(mInd(iNod)-1)+2)=dN1(iNod);
      end
      for iNod=1:nEltNod
        B(2,nColDof*(nInd(iNod)-1)+1)=dN2(iNod);
        B(4,nColDof*(nInd(iNod)-1)+2)=dN2(iNod);
      end

      C=[lambda+2*mu lambda      0    % Constitutive equation
        lambda       lambda+2*mu 0
        0            0           mu];
      sigma=C*A*T*B;
      
      % Tractions
      N1=[n1(1) 0     n1(3);
          0     n1(3) n1(1)];
      N2=[n2(1) 0     n2(3);
          0     n2(3) n2(1)];
      Tc(nColDof*(mCorner-1)+[1 2],:)=N1*sigma;
      Tc(nColDof*(nCorner-1)+[1 2],:)=N2*sigma;

    elseif nColDof==3  % 2.5D ROUTINE
      T=zeros(6,6);
      T(1:2,1:2)=dxidx; % Transformation matrix from xi to x.
      T(3:4,3:4)=dxidx;
      T(5:6,5:6)=dxidx;
      
      A1=[1 0 0 0 0 0
          0 0 0 0 0 0
          0 0 0 0 0 1
          0 0 1 0 0 0
          0 0 0 1 0 0
          0 1 0 0 1 0];
      A2=zeros(6,nColDof*nCol);
      for iNod=1:mEltNod
        A2(2,nColDof*(mInd(iNod)-1)+2)=-i*ky;
        A2(4,nColDof*(mInd(iNod)-1)+1)=-i*ky;
        A2(5,nColDof*(mInd(iNod)-1)+3)=-i*ky;
      end
      for iNod=1:nEltNod
        A2(2,nColDof*(nInd(iNod)-1)+2)=-i*ky;
        A2(4,nColDof*(nInd(iNod)-1)+1)=-i*ky;
        A2(5,nColDof*(nInd(iNod)-1)+3)=-i*ky;
      end      

      C=[lambda+2*mu  lambda       lambda        0   0   0    % Constitutive equation
         lambda       lambda+2*mu  lambda        0   0   0
         lambda       lambda       lambda+2*mu   0   0   0
         0            0            0             mu  0   0
         0            0            0             0   mu  0
         0            0            0             0   0   mu];

      B=zeros(6,nColDof*nCol);
      for iNod=1:mEltNod
        B(1,nColDof*(mInd(iNod)-1)+1)=dN1(iNod);
        B(3,nColDof*(mInd(iNod)-1)+2)=dN1(iNod);
        B(5,nColDof*(mInd(iNod)-1)+3)=dN1(iNod);
      end
      for iNod=1:nEltNod
        B(2,nColDof*(nInd(iNod)-1)+1)=dN2(iNod);
        B(4,nColDof*(nInd(iNod)-1)+2)=dN2(iNod);
        B(6,nColDof*(nInd(iNod)-1)+3)=dN2(iNod);
      end
      
      sigma=C*(A1*T*B+A2);
      
      % Tractions
      N1=[n1(1) 0  0     0     0      n1(3);
          0     0  0     n1(3) n1(1)  0;
          0     0  n1(3) 0     0      n1(1)];
      N2=[n2(1) 0  0     0     0      n2(3);
          0     0  0     n2(3) n2(1)  0;
          0     0  n2(3) 0     0      n2(1)];
      Tc(nColDof*(mCorner-1)+[1 2 3],:)=N1*sigma;
      Tc(nColDof*(nCorner-1)+[1 2 3],:)=N2*sigma;
    end

  elseif dim==3 % 3D ROUTINE
    for iDir=1:2
      for jDir=1:2
        cp=sum(cross(tan1(:,iDir),tan2(:,jDir)));
        if abs(cp)<1e-10
          commonEdge=jDir;
        end
      end
    end
    if commonEdge==1, xi3ind=2; end
    if commonEdge==2, xi3ind=1; end

    dxdxi=[tan1(:,1).';tan1(:,2).';tan2(:,xi3ind).'];
    dxidx=inv(dxdxi);
    T=zeros(9,9);    % Transformation matrix from xi to x.
    T(1:3,1:3)=dxidx;
    T(4:6,4:6)=dxidx;
    T(7:9,7:9)=dxidx;
    A=[1 0 0 0 0 0 0 0 0; % Gradient to strain
       0 0 0 0 1 0 0 0 0;
       0 0 0 0 0 0 0 0 1;
       0 1 0 1 0 0 0 0 0;
       0 0 0 0 0 1 0 1 0;
       0 0 1 0 0 0 1 0 0];
    B=zeros(9,nColDof*nCol);
    for iNod=1:mEltNod
      B(1,nColDof*(mInd(iNod)-1)+1)=dN1(iNod,1);
      B(4,nColDof*(mInd(iNod)-1)+2)=dN1(iNod,1);
      B(7,nColDof*(mInd(iNod)-1)+3)=dN1(iNod,1);
      B(2,nColDof*(mInd(iNod)-1)+1)=dN1(iNod,2);
      B(5,nColDof*(mInd(iNod)-1)+2)=dN1(iNod,2);
      B(8,nColDof*(mInd(iNod)-1)+3)=dN1(iNod,2);
    end
    for iNod=1:nEltNod
      B(3,nColDof*(nInd(iNod)-1)+1)=dN2(iNod,xi3ind);
      B(6,nColDof*(nInd(iNod)-1)+2)=dN2(iNod,xi3ind);
      B(9,nColDof*(nInd(iNod)-1)+3)=dN2(iNod,xi3ind);
    end

    C=[lambda+2*mu  lambda       lambda        0   0   0    % Constitutive equation
       lambda       lambda+2*mu  lambda        0   0   0
       lambda       lambda       lambda+2*mu   0   0   0
       0            0            0             mu  0   0
       0            0            0             0   mu  0
       0            0            0             0   0   mu];
    sigma=C*A*T*B;

    N1=[n1(1) 0     0      n1(2)  0     n1(3)
        0     n1(2) 0      n1(1)  n1(3) 0
        0     0     n1(3)  0      n1(2) n1(1)];
    N2=[n2(1) 0     0      n2(2)  0     n2(3)
        0     n2(2) 0      n2(1)  n2(3) 0
        0     0     n2(3)  0      n2(2) n2(1)];
    Tc(nColDof*(mCorner-1)+[1 2 3],:)=N1*sigma;
    Tc(nColDof*(nCorner-1)+[1 2 3],:)=N2*sigma;
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