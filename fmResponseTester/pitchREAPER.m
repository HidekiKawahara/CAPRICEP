function output = pitchREAPER(xa, fs)
% Interface function for pitch function of Google REAPER
%   output = pitchREAPER(xa, fs)
%     Use the function name "@pitchREAPER" for the argument of the evaluator
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
audiowrite("testSignal.wav", xa/max(abs(xa))*0.8, fs,"BitsPerSample",16);
!~/Downloads/REAPER-master/build/reaper -i testSignal.wav -f ttt.fo -a
fid = fopen('ttt.fo','r');
tline = fgetl(fid);
while string(tline) ~= "EST_Header_End"
    tline = fgetl(fid);
end
nmax = 1000000;
tt = zeros(nmax, 1);
fo = zeros(nmax, 1);
nid = 0;
while nid < nmax
    tmp = fscanf(fid, '%f',3);
    if length(tmp) <3
        break
    end
    nid = nid + 1;
    tt(nid) = tmp(1);
    fo(nid) = tmp(3);
end
fclose(fid);
tt = tt(1:nid);
fo = fo(1:nid);
%
output.fo = fo;
output.tt = tt;
output.titleStr = "REAPER ";
output.filePrefix = "pREAPER";
end