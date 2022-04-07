function output = pitchPEF(xa, fs)
% Interface function for pitch function of MATLAB with SRH option
%   output = pitchPEF(xa, fs)
%     Use the function name "@pitchPEF" for the argument of the evaluator
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
[f0,loc] = pitch(xa, fs, "Method","PEF");
output.fo = f0;
output.tt = loc/fs;
output.titleStr = "PEF ";
output.filePrefix = "pPEF";
end