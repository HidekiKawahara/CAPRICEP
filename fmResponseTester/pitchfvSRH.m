function output = pitchfvSRH(xa, fs)
% Interface function for pitch extractor of SRH in COVREP
%   output = pitchYANG(xa, fs)
%     Use the function name "@pitchYANG" for the argument of the evaluator
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
%sourceStr = AnalyzeSpeechSource(xa, fs);
[srh_f0,~,~,srh_time] = pitch_srh(xa,fs,50,500,5);
output.fo = srh_f0(:);
output.tt = srh_time(:);
output.titleStr = "cSRH ";
output.filePrefix = "pcSRHG";
end