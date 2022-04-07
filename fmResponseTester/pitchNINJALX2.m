function output = pitchNINJALX2(xa, fs)
% Interface function for pitch extractor designed for CSJ-corpus analysis
%   output = pitchNINJALX2(xa, fs)
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

% sourceAttributesAnalysis(x, fs, range, low_frequency, high_freuency, ...
%     channels_in_octave, dsOn, stretching_factor, wintype, ...
%     integration_time)

% LICENSE: refer to LICENSE in this folder

output = struct;
range = [1 length(xa)];
low_frequency = 55;
high_freuency = 1200;
channels_in_octave = 24;
dsOn = 1;
stretching_factor = 1.0;
wintype = 'sixterm';
integration_time = 10; % This is changed from the default 40
foOut = sourceAttributesAnalysis(xa, fs, range, low_frequency, high_freuency, ...
     channels_in_octave, dsOn, stretching_factor, wintype, integration_time);
output.fo = foOut.fixed_points_freq(:,1);
output.tt = foOut.time_axis_wavelet(:);
output.titleStr = "NINJALX2 ";
output.filePrefix = "pNINJALX2";
end