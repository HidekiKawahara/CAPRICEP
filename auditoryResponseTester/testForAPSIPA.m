%%   fundamental component cancelling test


close all
clear variables
%%

baseDir = '/Volumes/HD-CD-A/TAFlog/';
fullList = {'/Volumes/HD-CD-A/TAFlog/audResp20210702T184558.wav', ...
    '/Volumes/HD-CD-A/TAFlog/audResp20210702T184658.wav', ...
    '/Volumes/HD-CD-A/TAFlog/audResp20210711T035454.wav', ...
    '/Volumes/HD-CD-A/TAFlog/audResp20210717T202331.wav', ...
    '/Volumes/HD-CD-A/TAFlog/audResp20210717T202520.wav', ...
    '/Volumes/HD-CD-A/TAFlog/audResp20210717T202606.wav', ...
    '/Volumes/HD-CD-A/TAFlog/audResp20210717T202744.wav', ...
    '/Volumes/HD-CD-A/TAFlog/audResp20210717T202831.wav'};

%% for SINE

audioinfo(fullList{4})
[x, fs] = audioread(fullList{4});
tx = (1:length(x))'/fs;
figure;
subplot(211)
plot(tx, x(:, 1));grid on;
subplot(212)
plot(tx, x(:, 2));grid on;

%%
idx1s = (1:fs);
initP = fs * 18;
scanL = round(fs/260*1.5);
diffV = zeros(scanL, 1);
for ii = 1:scanL
    diffV(ii) = sum(abs(x(initP + idx1s, 1) - x(initP + idx1s + ii, 2)));
end
[~, minIdx] = min(diffV);

%%
lvl1 = std(x(initP + idx1s, 1));
lvl2 = std(x(initP + idx1s + minIdx, 2));

modSig = x(:, 2) / lvl2 * lvl1;
testSig = x(fs + (1:18*fs), 1) - modSig(fs + (1:18*fs) + minIdx, 1);

%%

scanL2 = round(0.02*fs);
scanIdx = (-scanL2:scanL2);
diffV2 = zeros(length(scanIdx), 1);
for ii = 1:length(scanIdx)
    diffV2(ii) = std(x(fs*15 + (1:3*fs), 1) - modSig(fs*15 + (1:3*fs)+scanIdx(ii), 1));
end

%%
[~, minIdx2] = min(diffV2);
realIdx = scanIdx(minIdx2);

%%
magList = 0.9:0.01:1.2;
diffV3 = zeros(length(magList), 1);
for ii = 1:length(magList)
    diffV3(ii) = std(x(fs*15 + (1:3*fs), 1) - magList(ii) * modSig(fs*15 + (1:3*fs)+realIdx, 1));
end
%%
[~, magIdx] = min(diffV3);
fixedSignal = x(fs + (1:18*fs), 1) - magList(magIdx) * modSig(fs + (1:18*fs)+realIdx, 1);
fixedSignal2 = x(fs + (1:18*fs), 1) - modSig(fs + (1:18*fs)+realIdx, 1);
