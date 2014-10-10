function B=bemcompress(A,tresh)
% BEMCOMPRESS   Compress time boundary element matrices.
%   
%   B = bemcompress(A,tresh) compresses the boundary element
%   system matrix A by replacing elements below a treshhold value by zero. 
%   The result is stored as the sparse matrix B.
%
%   A     boundary element matrix (nDof * nDof).
%   tresh treshold value.
%   B     compressed, sparse boundary element matrix (nDof * nDof).
%
%   This routine may only be used in the case of time domain boundary elements,
%   where a large part of the system matrices is zero due to a finite wave
%   propagation speed and causality.

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

Asize=size(A);
tresh=abs(tresh);

m=Asize(1);
n=Asize(2);

% Count number of non-zero values
refval=max(max(abs(A)));
refval=max(refval,10*realmin);
minval=refval*tresh;
count=0;
for iA=1:m
  for jA=1:n
    val=A(iA,jA);
    if (abs(val) >  minval )
      count=count+1;
    end
  end
end

for iA=1:m
  for jA=1:n
    val=A(iA,jA);
    if (abs(val) >  minval )
      count=count+1;
      sval(count)=val;
      iindex(count)=iA;
      jindex(count)=jA;
    end
  end
end

B=sparse(iindex(1:count),jindex(1:count),sval(1:count),m,n,count);
