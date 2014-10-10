function bemfunmake(outdir)

%BEMFUNMAKE   Build the MEX-files in the ElastoDynamics Toolbox.
%
%   BEMFUNMAKE(outdir) compiles and links the MEX-files in the
%   Boundary Elements Toolbox.  The resulting files are moved to the specified
%   target directory.

% Stijn Francois
% Mattias Schevenels
% June 2008

try
  if ~exist(outdir,'dir'), mkdir(outdir); end
%   pmake('bemcollindex',outdir); 
%   pmake('bemcompress',outdir); 
%   pmake('bemcorner',outdir); 
%   pmake('bemeltreverse',outdir); 
%   pmake('bemexport',outdir); 
%   pmake('bemfuncheck',outdir); 
%   pmake('bemgrid_cyl',outdir); 
%   pmake('bemgrid_rec',outdir); 
%   pmake('bemgrid_pol',outdir); 
%   pmake('bemgridplot_cyl',outdir); 
%   pmake('bemgridplot_rec',outdir); 
%   pmake('bemgridplot_pol',outdir); 
%   pmake('bemimport',outdir); 
%   pmake('bemmesharc',outdir); 
%   pmake('bemmesharea',outdir); 
%   pmake('bemmeshcat',outdir); 
%   pmake('bemmeshcircle',outdir); 
%   pmake('bemmeshcylinder',outdir); 
%   pmake('bemmeshdisk',outdir); 
%   pmake('bemmeshline',outdir);
%   pmake('bemmeshrep',outdir);  
%   pmake('bemmeshsphere',outdir); 
%   pmake('bemplot',outdir); 
%   pmake('bempressure',outdir); 
%   pmake('bemrigid',outdir); 
%   pmake('bemtq',outdir); 
%   pmake('bemtu',outdir); 
%   pmake('newmarkcoef',outdir); 
%   pmake('newmarkforce',outdir); 
%   pmake('newmarkstiff',outdir); 
%   pmake('newmarkupdate',outdir); 
%   pmake('timflexmat',outdir); 
%   pmake('timflexmat2',outdir); 
%   pmake('selectdof',outdir); 
%   pmake('unselectdof',outdir); 

  try bemfunlicense('reset'); end
  clear bemfunlicense;
  try if exist(which('bemfunlicense'),'file'), delete(which('bemfunlicense')); end; end
  pause(0.1);
  try bemfunlicense('reset'); end
  clear bemfunlicense;
  compile('bigint.cpp');
  compile('ripemd128.cpp');
  compile('rsa.cpp');
  compile('getmac.cpp');
  compile('bemfunlicense.cpp')
  compile('bemfunlicense_mex.cpp');
  compile('hashfile_mex.cpp');
  link(sprintf('%s/bemfunlicense',outdir),'bemfunlicense_mex.o','bemfunlicense.o','bigint.o','ripemd128.o','rsa.o','getmac.o');
  link('hashfile','hashfile_mex.o','ripemd128.o');

  fprintf('Computing hash code for bemfunlicense.%s ...\n',mexext);
  hash=hashfile(sprintf('%s/bemfunlicense.%s',outdir,mexext));
  fprintf('Updating checklicense.cpp ...\n');
  lic=txtread('checklicense.cpp');
  lic=regexprep(lic,'(#define BEMFUNLICENSE_HASH) "\w*"',sprintf('$1 "%s"',hash));
  txtwrite('checklicense.cpp',lic);
  fprintf('Deleting hashfile.%s and hashfile.m ...\n',mexext);
  try
    clear('hashfile');
    delete('hashfile.m');
    delete(['hashfile.' mexext]);
  end

  
  compile('bemmat_mex.cpp');
  compile('s2coll.cpp');
  compile('uniquecoll.cpp');
  compile('bemmat.cpp');
  compile('greeneval3d.cpp');
  compile('bemtangent_mex.cpp');
  compile('bemshapederiv_mex.cpp');
  compile('bemcollpoints.cpp');
  compile('bemcollpoints_mex.cpp');
  compile('bemdimension.cpp');
  compile('bemdimension_mex.cpp');
  compile('bemeltdef_mex.cpp');
  compile('bemint.cpp');
  compile('bemint_mex.cpp');
  compile('bemintpoints_mex.cpp');
  compile('bemintreg2d.cpp');
  compile('bemintreg3d.cpp');
  compile('bemintreg3dnodiag.cpp');
  compile('bemintreg3ddiag.cpp');
  compile('bemintreg3dperiodic.cpp');
  compile('bemintregaxi.cpp');
  compile('bemintsing2d.cpp');
  compile('bemintsing3d.cpp');
  compile('bemintsing3dperiodic.cpp');
  compile('bemintsingaxi.cpp');
  compile('bemisaxisym.cpp');
  compile('bemisaxisym_mex.cpp');
  compile('bemisperiodic.cpp');
  compile('bemisperiodic_mex.cpp');
  compile('bemmatconv_mex.cpp');
  compile('bemnormal.cpp');
  compile('bemnormal_mex.cpp');
  compile('bemshape_mex.cpp');
  compile('bemtimeconv_mex.cpp');
  compile('bemxfer2d.cpp');
  compile('bemxfer3d.cpp');
  compile('bemxfer3dperiodic.cpp');
  compile('bemxfer_mex.cpp');
  compile('bemxferaxi.cpp');
  compile('besselh.cpp');
  compile('boundaryrec2d.cpp');
  compile('boundaryrec3d.cpp');
  compile('checklicense.cpp');
  compile('eltdef.cpp');
  compile('fminstep.cpp');
  compile('fsgreen2d_inplane.cpp');
  compile('fsgreen2d_outofplane.cpp');
  compile('fsgreen3d.cpp');
  compile('fsgreen3dt.cpp');
  compile('fsgreenf.cpp');
  compile('gausspw.cpp');
  compile('gausspw1d_mex.cpp');
  compile('gausspw2d_mex.cpp');
  compile('greeneval2d.cpp');
  compile('greenrotate2d.cpp');
  compile('greenrotate3d.cpp');
  compile('search1.cpp');
  compile('shapefun.cpp');
  
  link(sprintf('%s/bemcollpoints',outdir),'bemcollpoints_mex.o','eltdef.o','bemcollpoints.o','shapefun.o','bemdimension.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemdimension',outdir),'bemdimension_mex.o','bemdimension.o','eltdef.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemeltdef',outdir),'bemeltdef_mex.o','eltdef.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemint',outdir),'bemint_mex.o','eltdef.o','bemcollpoints.o','shapefun.o','gausspw.o','bemint.o','bemisaxisym.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemintpoints',outdir),'bemintpoints_mex.o','eltdef.o','gausspw.o','bemcollpoints.o','shapefun.o','bemdimension.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemisaxisym',outdir),'bemisaxisym_mex.o','bemisaxisym.o','eltdef.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemisperiodic',outdir),'bemisperiodic_mex.o','bemisperiodic.o','eltdef.o','checklicense.o','ripemd128.o');
% %   link(sprintf('%s/bemmat',outdir),'bemmat_mex.o','bemmat.o','s2coll.o','uniquecoll.o','eltdef.o','bemcollpoints.o','shapefun.o','bemintreg3d.o','bemintreg3dperiodic.o','bemintreg2d.o','bemintregaxi.o','bemintsing3d.o','bemintsing3dperiodic.o','bemintsing2d.o','bemintsingaxi.o','gausspw.o','search1.o','bemnormal.o','bemdimension.o','bemisaxisym.o','bemisperiodic.o','greeneval2d.o','greeneval3d.o','greenrotate2d.o','greenrotate3d.o','fsgreenf.o','fsgreen2d_inplane.o','fsgreen2d_outofplane.o','besselh.o','fsgreen3d.o','fsgreen3dt.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemmat',outdir),'bemmat_mex.o','bemmat.o','s2coll.o','uniquecoll.o','eltdef.o','bemcollpoints.o','shapefun.o','bemintreg3dnodiag.o','bemintreg3ddiag.o','bemintreg3dperiodic.o','bemintreg2d.o','bemintregaxi.o','bemintsing3d.o','bemintsing3dperiodic.o','bemintsing2d.o','bemintsingaxi.o','gausspw.o','search1.o','bemnormal.o','bemdimension.o','bemisaxisym.o','bemisperiodic.o','greeneval2d.o','greeneval3d.o','greenrotate2d.o','greenrotate3d.o','fsgreenf.o','fsgreen2d_inplane.o','fsgreen2d_outofplane.o','besselh.o','fsgreen3d.o','fsgreen3dt.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemmatconv',outdir),'bemmatconv_mex.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemnormal',outdir),'bemnormal_mex.o','eltdef.o','shapefun.o','bemnormal.o','bemcollpoints.o','bemdimension.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemtangent',outdir),'bemtangent_mex.o','eltdef.o','shapefun.o','bemnormal.o','bemcollpoints.o','bemdimension.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemshape',outdir),'bemshape_mex.o','shapefun.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemshapederiv',outdir),'bemshapederiv_mex.o','shapefun.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemtimeconv',outdir),'bemtimeconv_mex.o','search1.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/bemxfer',outdir),'bemxfer_mex.o','eltdef.o','bemcollpoints.o','shapefun.o','bemnormal.o','gausspw.o','search1.o','bemxfer3d.o','bemxfer3dperiodic.o','bemxfer2d.o','bemxferaxi.o','bemdimension.o','bemisaxisym.o','bemisperiodic.o','greeneval2d.o','fsgreenf.o','fsgreen3d.o','fsgreen3dt.o','fsgreen2d_inplane.o','fsgreen2d_outofplane.o','besselh.o','greeneval3d.o','greenrotate2d.o','boundaryrec2d.o','boundaryrec3d.o','fminstep.o','greenrotate3d.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/gausspw1d',outdir),'gausspw1d_mex.o','gausspw.o','checklicense.o','ripemd128.o');
  link(sprintf('%s/gausspw2d',outdir),'gausspw2d_mex.o','gausspw.o','checklicense.o','ripemd128.o');
  
  finalize
catch
  finalize
  error(lasterr);
end
%-------------------------------------------------------------------------------

function pmake(fun,outdir)
fprintf('Compiling %s.m ...\n',fun);
mfile=txtread(sprintf('%s.m',fun));
mfile=strrep(mfile,char([13 10]),char(10));
mfile=mfile(find(mfile=='%',1):end);
k=regexp(mfile,[10 '[^\%]']);
mfile=mfile(1:k(2)-1);
mfile=strtrim(mfile);
mfile=[mfile char(10)];
txtwrite(sprintf('%s/%s.m',outdir,fun),mfile);
pcode(fun);
movefile(sprintf('%s.p',fun),outdir);

%-------------------------------------------------------------------------------

% % % % function compile(cppFile)
% % % % fprintf('Compiling %s ...\n',cppFile);
% % % % % if ~exist('verLessThan','file') || verLessThan('matlab','7.2')
% % % % if ~exist('verLessThan','file') || verLessThan('matlab','7.2')
% % % %     %   mex('-O','-c',cppFile);
% % % %     mex('-O','-c','CFLAGS=-O2 -g -DDEBUG',cppFile);
% % % % %     mex('-O','-c','-g -DDEBUG',cppFile);
% % % % else
% % % %     if strcmp(cppFile,'bemmat_mex.cpp') || strcmp(cppFile,'bemmat.cpp')
% % % % %     if  strcmp(cppFile,'bemmat.cpp')
% % % %         %   mex('-O','-c','-largeArrayDims',cppFile);
% % % % %         mex('-O','-c','-g','-f "/home/zuhal/pieter/NewDirectory/source/mexopts.bat"','-largeArrayDims',cppFile);
% % % % % mex('-O','-c','-g','-f mexopts.sh','-largeArrayDims',cppFile);
% % % %         mex('-O','-c','CXXOPTIMFLAGS="$CXXOPTIMFLAGS -O"','CXXDEBUGFLAGS="$CXXDEBUGFLAGS -g -DDEBUG"','-largeArrayDims',cppFile);
% % % % %         mex('-O','-c','-g -DDEBUG','-largeArrayDims',cppFile);
% % % %     else
% % % %         mex('-O','-c','-largeArrayDims',cppFile);
% % % %     end
% % % %         
% % % % end
% % % % [dum,file]=fileparts(cppFile);
% % % % objFile=strcat(file,'.obj');
% % % % oFile=strcat(file,'.o');
% % % % if exist(objFile,'file')
% % % %   movefile(objFile,oFile);
% % % % end

function compile(cppFile)
fprintf('Compiling %s ...\n',cppFile);
if ~exist('verLessThan','file') || verLessThan('matlab','7.2')
  mex('-O','-c',cppFile);
else
  mex('-O','-c','-largeArrayDims',cppFile);
end
[dum,file]=fileparts(cppFile);
objFile=strcat(file,'.obj');
oFile=strcat(file,'.o');
if exist(objFile,'file')
  movefile(objFile,oFile);
end

%-------------------------------------------------------------------------------
function link(varargin)

% ACTUAL LINKING
[outDir,funName]=fileparts(varargin{1});
fprintf('Linking %s.%s ...\n',funName,mexext);
if ispc
  if strcmpi(funName,'bemfunlicense')
    if ~exist('verLessThan','file') || verLessThan('matlab','7.2')
      mex -O -g -output bemfunlicense GCCLIBS="$GCCLIBS -liphlpapi" bemfunlicense_mex.o bemfunlicense.o bigint.o ripemd128.o rsa.o getmac.o
    else
      mex -largeArrayDims -O -g -output bemfunlicense GCCLIBS="$GCCLIBS -liphlpapi" bemfunlicense_mex.o bemfunlicense.o bigint.o ripemd128.o rsa.o getmac.o
    end
    movefile(['bemfunlicense.' mexext],outDir);
  else
    if ~exist('verLessThan','file') || verLessThan('matlab','7.2')
      mex('-O','-output','-g',varargin{:});
    else
%       mex('-O','-largeArrayDims','-output',varargin{:});
        mex('-O','CXXOPTIMFLAGS="$CXXOPTIMFLAGS -O"','CXXDEBUGFLAGS="$CXXDEBUGFLAGS -g -DDEBUG"','-largeArrayDims','-output',varargin{:});
%         mex('-O','-g -DDEBUG','-largeArrayDims','-output',varargin{:});
    end
  end
else
  if ~exist('verLessThan','file') || verLessThan('matlab','7.2')
    mex('-cxx','-O','-output','-g',varargin{:});
  else
%     mex('-cxx','-O','-largeArrayDims','-output',varargin{:});
      mex('-cxx','-O','CXXOPTIMFLAGS="$CXXOPTIMFLAGS -O"','CXXDEBUGFLAGS="$CXXDEBUGFLAGS -g -DDEBUG"','-largeArrayDims','-output',varargin{:});
%     mex('-cxx','-O','-g -DDEBUG','-largeArrayDims','-output',varargin{:});
  end
end

% EXTRACT HELP SECTION
srcFile=varargin{~cellfun('isempty',strfind(varargin,'_mex.o'))};
srcFile=strrep(srcFile,'.o','.cpp');
mFile=[varargin{1} '.m'];
srcFid = fopen(srcFile);
mFid = fopen(mFile,'w');
help = 0;
while 1
  line = fgetl(srcFid);
  if ~ischar(line), break; end;
  line = strtrim(line);
  if strcmp(line(1:min(2,end)),'/*') && (help == 0), help = 1; end;
  if strcmp(line(1:min(2,end)),'*/'), help = 2; end;
  if help == 1
    if strcmp(line(1:min(2,end)),'/*'), line = line(2:end); end;
    if line(1:min(1,end))=='*', line(1) = '%'; end;
    fprintf(mFid,'%s\n',line);
  end
end
fclose(srcFid);
fclose(mFid);

%-------------------------------------------------------------------------------
function finalize
fprintf('Deleting object files ...\n');
try
  delete('*.o');
end
try
  delete('*.obj');
end

%-------------------------------------------------------------------------------
function s=txtread(file)
fid=fopen(file,'r');
s=char(fread(fid,'uchar')');
fclose(fid);

%-------------------------------------------------------------------------------
function txtwrite(file,s)
fid=fopen(file,'w');
fwrite(fid,s,'uchar');
fclose(fid);
