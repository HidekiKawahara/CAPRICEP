function output = pitchPraatDefault(xa, fs)
% Interface function for the default pitch extractor of Praat
% This is for macOS. Plase referr Praat's help for other machines
% This function also uses a function "praatFoFileRead.m" 
%
%   output = pitchPraatDefault(xa, fs)
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
audiowrite("testSignal.wav", xa/max(abs(xa))*0.8, 44100, "BitsPerSample", 24);
fNamePrScript = "cppMat.praat";
fid = fopen(fNamePrScript,"W");
cmdStr1 = ['Read from file: "' pwd '/testSignal.wav"'];
cmdStr2 = 'To Pitch: 0, 75, 600';
cmdStr3 = ['Save as text file: "' pwd '/testSignalFo.txt"'];
fprintf(fid, '%s\n', cmdStr1);
fprintf(fid, '%s\n', cmdStr2);
fprintf(fid, '%s\n', cmdStr3);
fclose(fid);
%type cppMat.praat (for debugging)
% Following line is for batch processing of Praat on macOS
!/Applications/Praat.app/Contents/MacOS/Praat --run cppMat.praat >> praatLog.txt
foo = praatFoFileRead([ pwd '/testSignalFo.txt']);
output.fs = fs;
output.fo = foo.fo;
output.tt = (1:foo.n)'*foo.dtx;
output.titleStr = "praatDef ";
output.filePrefix = "praatDef";
end
