function [L,I]=selectdof(DOF,seldof,varargin)

%SELECTDOF   Select degrees of freedom.
%
%   [L,I]=selectdof(DOF,seldof,'Ordering',ordering) creates the matrix to 
%   extract degrees of freedom from the global degrees of freedom by matrix 
%   multiplication.
%
%   DOF        Degrees of freedom  (nDOF * 1).
%   seldof     Selected DOF labels (kDOF * 1).
%   L          Selection matrix (kDOF * nDOF).
%   I          Index vector (kDOF * 1).
%   ordering   'seldof','DOF' or 'sorted'. Default: 'seldof'
%              Ordering of L and I similar as seldof, DOF or sorted.

% David Dooms
% March 2008

bemfunlicense('VerifyOnce');

% PREPROCESSING
DOF=DOF(:);
seldof=seldof(:);
if nargin<3
    varargin={};
end
paramlist=varargin;
[ordering,paramlist]=cutparam('Ordering','seldof',paramlist,{'seldof','DOF','sorted'});
if ~isempty(paramlist) error(['Undefined parameter ''' paramlist{1} '''.']); end
DOF=DOF(:);
seldof=seldof(:);

if ~ isempty(find(seldof==0.00))
    error('The wild card 0.00 is not allowed')
end

nDOF=length(DOF);
kDOF=length(seldof);
indj=zeros(1,kDOF);
nterm=0;
for idof=1:kDOF
    if floor(seldof(idof,1))==0       % wild cards 0.0X
        indjdof=find(abs(rem(DOF,1)-rem(seldof(idof,1),1))<0.0001);
    elseif rem(seldof(idof,1),1)==0   % wild cards X.00
        indjdof=find(abs(floor(DOF)-floor(seldof(idof,1)))<0.0001);
    else                                 % standard case
        indjdof=find(abs(DOF-seldof(idof,1))<0.0001);
    end
    if isempty(indjdof)
        indj(1,nterm+1)=NaN;
        nterm=nterm+1;
    else
        indj(1,nterm+1:nterm+length(indjdof))=indjdof.';
        nterm=nterm+length(indjdof);
    end
end  

indi=find(~isnan(indj));

if strcmpi(ordering,'seldof')                    % default: ordered as in seldof
    I=indj.';
    indj=indj(indi);
    s=ones(size(indi));
    kDOF=nterm;
elseif strcmpi(ordering,'DOF')                   % ordered as in DOF
    indj=indj(indi);
    indj=unique(indj);
    I=indj.';
    kDOF=length(indj);
    indi=1:kDOF;
    s=ones(1,kDOF);
elseif strcmpi(ordering,'sorted')                % sorted
    indj=indj(indi);
    [sortseldof,indDOF]=unique(DOF(indj));
    indj=indj(indDOF);
    I=indj.';
    kDOF=length(indj);
    indi=1:kDOF;
    s=ones(1,kDOF);
end

L=sparse(indi,indj,s,kDOF,nDOF);
I=I(:);

%-------------------------------------------------------------------------------
% CUT PARAMETER FROM LIST
function [value,paramlist]=cutparam(name,default,paramlist,allowed)
value=default;
for iarg=length(paramlist)-1:-1:1
  if strcmpi(name,paramlist{iarg})
    value=paramlist{iarg+1};
    paramlist=paramlist([1:iarg-1 iarg+2:end]);
    break
  end
end
if nargin==4
  if isempty(find(strcmpi(value,allowed)))
    error(['Parameter ' name ' should be one of ' sprintf('''%s'' ',allowed{:}) '.']);
  end
end
