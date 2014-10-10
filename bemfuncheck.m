function bemfuncheck

%BEMFUNCHECK   Check whether the MEX files in BEMFUN are properly compiled.
%
%   BEMFUNCHECK attempts to run one of the MEX files in BEMFUN to check whether
%   these files are available in compiled form for the current platform.
%   The result is displayed in the MATLAB command window.

% Stijn François
% July 2008

try
  bemeltdef(1,{1 'line2' {'nodcol'}});
  fprintf('The MEX files in BEMFUN are properly compiled for the current platform.\nBEMFUN is ready to be used.\n')
catch
  fprintf('The MEX files in BEMFUN are not yet compiled for the current platform.\nPlease contact the authors of BEMFUN on http://www.kuleuven.be/bwm/bemfun.\n');
end