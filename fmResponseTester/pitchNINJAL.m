function output = pitchNINJAL(xa, fs)
% Interface function for pitch extractor designed for CSJ-corpus analysis
%   output = pitchNINJAL(xa, fs)
%     Use the function name "@pitchNINJAL" for the argument of the evaluator
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
foOut = sourceAttributesAnalysis(xa, fs);
output.fo = foOut.fixed_points_freq(:,1);
output.tt = foOut.time_axis_wavelet(:);
output.titleStr = "NINJAL ";
output.filePrefix = "pNINJAL";
end