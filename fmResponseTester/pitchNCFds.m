function output = pitchNCFds(xa, fs)
% Interface function for pitch function of MATLAB using down sampling
% This shows how to convert the test signal (sampling at 44100 Hz) to
% a test signal with lower sampling frequency. This example uses
% 44100/4 Hz.
%
%   output = pitchNCFds(xa, fs)
%     Use the function name "@pitchNCF" for the argument of the evaluator
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
xad = decimate(xa, 4);
[f0, loc] = pitch(xad, fs/4);
output.fo = f0;
output.tt = loc/(fs/4);
output.titleStr = "NCFdf ";
output.filePrefix = "pNCFdf";
end