function output = pitchOpenSMILE(xa, fs)
% Interface function for pitch function of MATLAB
%   output = pitchNCF(xa, fs)
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
%[f0, loc] = pitch(xa, fs);
audiowrite("testSignal.wav", xa/max(abs(xa))*0.8, fs,"BitsPerSample",16);
!~/Downloads/opensmile-3.0.1-macos-x64/bin/SMILExtract -C ~/Downloads/opensmile-3.0.1-macos-x64/config/prosody/prosodyAcf.conf -I testSignal.wav -O testSignal.htk
htkData = v_readhtk('testSignal.htk');
output.fo = htkData(:,2);
output.tt = (1:length(output.fo))'/100;
output.titleStr = "openSMILE ";
output.filePrefix = "pOSMILE";
end