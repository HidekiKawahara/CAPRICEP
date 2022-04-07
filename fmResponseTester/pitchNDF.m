function output = pitchNDF(xa, fs)
% Interface function for pitch extractor of legacy-STRAIGHT option
%   output = pitchNDF(xa, fs)
%     Use the function name "@pitchNDF" for the argument of the evaluator
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
[f0raw,vuv,auxouts]=MulticueF0v14(xa,fs,70,450);
output.fo = f0raw;
output.tt = (1:length(f0raw))'/1000;
output.titleStr = "NDF ";
output.filePrefix = "pNDF";
end