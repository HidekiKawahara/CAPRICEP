function output = pitchHarvest(xa, fs)
% Interface function for pitch extractor of WORLD
%   output = pitchHarvest(xa, fs)
%     Use the function name "@pitchHarvest" for the argument of the evaluator
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
option_harvest.f0_floor = 40;
f0_parameter = Harvest(xa, fs, option_harvest);
output.fo = f0_parameter.f0;
output.tt = f0_parameter.temporal_positions;
output.titleStr = "Harvest ";
output.filePrefix = "pHarvest";
end