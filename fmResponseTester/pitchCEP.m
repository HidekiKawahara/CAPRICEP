function output = pitchCEP(xa, fs)
% Interface function for pitch extractor of Noll's CEPSTRUM in MATLAB
%   output = pitchCEP(xa, fs)
%     Use the function name "@pitchCEP" for the argument of the evaluator
% Augment
%   xa   : test signal with CAPRICEP FM and simulatee /a/ spectrum
%   fs   : sampling frequency (Hz)
% Output
%   output : structure varialbe with the following fields
%      fo   : extracted fundamental frequency (Hz)
%      tt   : discrete temporal locations of fo measurement (s)
%      titleStr    : string for the first item of the figure title
%      filePrefix  : string for the beggining of the output files

% LICENSE: refer to LICENSE in this folder

output = struct;
[f0, loc] = pitch(xa, fs,"Method","CEP","Range",[70 450]);
output.fo = f0;
output.tt = loc/fs-0.028;
output.titleStr = "CEP ";
output.filePrefix = "pCEP";
end