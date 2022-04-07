function output = pitchXSX(xa, fs)
% Interface function for pitch extractor of TANDEM-STRAIGHT
%   output = pitchXSX(xa, fs)
%     Use the function name "@pitchXSX" for the argument of the evaluator
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
r = exF0candidatesTSTRAIGHTGB(xa,fs);
output.fo = r.f0;
output.tt = r.temporalPositions;
output.titleStr = "XSX ";
output.filePrefix = "pXSX";
end