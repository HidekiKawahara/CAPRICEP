%% Test script for two channel input

close all
clear variables

%%
fs = 44100;
nBits = 24;
inChannel = 2;
outChannel = 'L-ch';
nRepetition = 44;
tResponse = 400;
option = struct;

%% Note that the sound pressure level was 80 dB using A-weight

recTest2chStr = capricepResponseTest(fs, tResponse, nRepetition, outChannel, ...
    inChannel);
option.calibrationConst = 80 - recTest2chStr.lAeq;
recTest2chStr.lAeq = [80 80];

%%
if 1 == 2
save recTest2chStr recTest2chStr
end

%%

tmp = load('recTest2chStr');

analysisStr = tmp.recTest2chStr;

%%
specStr = capricepResponseAnalysis(analysisStr);

%%

repStr = capricepResponseReport(specStr);