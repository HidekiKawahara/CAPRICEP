function output = foResponseAnalysis(testSignalData, fopi, trueF0, fo, fs)
%baseIdx = (1:length(fopi))';
%[r, lags] = xcorr(fopi-fo, trueF0Org(baseIdx)-fo);
%[~, maxLag] = max(r);
%trueF0 = trueF0Org(baseIdx-lags(maxLag));
fopi = double(fopi);
%% select matching TSP
Bk = [1 1 1 1 1 1 1 1; ...
    1 -1 1 -1 1 -1 1 -1; ...
    1 1 -1 -1 1 1 -1 -1; ...
    1 1 1 1 -1 -1 -1 -1];
fftlSeg = testSignalData.fftlSeg;
setIds = testSignalData.cOrdered(1, :);
baseList = 1:10;
nonMember = baseList(~ismember(baseList, setIds));
xTSPset = testSignalData.xTSPSel(:, [setIds nonMember(1)]);
filtrdSig = fftfilt(flip(xTSPset, 1), 1200*log2(trueF0/fo));
filtrdSigMes = fftfilt(flip(xTSPset, 1), 1200*log2(fopi/fo));
idsWhole = 1:length(trueF0);
orthSig = filtrdSig;
orthSigMes = filtrdSigMes;
for ii = 2:8
    orthSig = orthSig + filtrdSig(circshift(idsWhole, (ii-1)*fftlSeg), :)*diag(Bk(:,ii));
    orthSigMes = orthSigMes + filtrdSigMes(circshift(idsWhole, (ii-1)*fftlSeg), :)*diag(Bk(:,ii));
end
impulseLong = orthSig * [1 -1 2 0]'/32;
impulseLongMes = orthSigMes * [1 -1 2 0]'/32;
orthSigMes = orthSigMes / 8;
orthSigRef = orthSig / 8;
maxinpulseLevel = max(impulseLong);
pulseIdx = 1:length(impulseLong);
pulseLocations = pulseIdx(impulseLong >= impulseLong([1 1:end-1]) ...
    & impulseLong > impulseLong([2:end end]));
pulseLocations = pulseLocations(impulseLong(pulseLocations) > 0.99 * maxinpulseLevel);
longPeriod = fftlSeg*4;
shortPeriod = fftlSeg;
baseLongRespIdx = (1:longPeriod) - round(fs*0.15);
baseShortRespIdx = (1:shortPeriod) - round(fs*0.15);
averageLongRef = baseLongRespIdx(:) * 0;
averageRandRef = baseLongRespIdx(:) * 0;
averageLongMes = baseLongRespIdx(:) * 0;
averageRandMes = baseLongRespIdx(:) * 0;
indivShortRef = baseShortRespIdx(:) * [0 0 0 0];
indivShortMes = baseShortRespIdx(:) * [0 0 0 0];
for ii = 2:length(pulseLocations) - 2
    averageLongRef = averageLongRef + impulseLong(baseLongRespIdx + pulseLocations(1) ...
        + (ii - 1) * longPeriod);
    averageRandRef = averageRandRef + orthSig(baseLongRespIdx + pulseLocations(1) ...
        + (ii - 1) * longPeriod, 4).^2;
    averageLongMes = averageLongMes + impulseLongMes(baseLongRespIdx + pulseLocations(1) ...
        + (ii - 1) * longPeriod);
    averageRandMes = averageRandMes + orthSigMes(baseLongRespIdx + pulseLocations(1) ...
        + (ii - 1) * longPeriod, 4).^2;
    indivShortRef = indivShortRef + orthSigRef(baseShortRespIdx + pulseLocations(1) ...
        + (ii - 1) * longPeriod, :);
    indivShortMes = indivShortMes + orthSigMes(baseShortRespIdx + pulseLocations(1) ...
        + (ii - 1) * longPeriod, :);
end
averageLongRef = averageLongRef / (length(pulseLocations) - 3);
averageRandRef = averageRandRef / (length(pulseLocations) - 3);
averageLongMes = averageLongMes / (length(pulseLocations) - 3);
averageRandMes = averageRandMes / (length(pulseLocations) - 3);
indivShortMes = indivShortMes / (length(pulseLocations) - 3);
indivShortRef = indivShortRef / (length(pulseLocations) - 3);
output.fx = (0:4*fftlSeg-1)/4/fftlSeg*fs;
dd = diag([1 -1 1]');
output.refGain = 20*log10(abs(fft(averageLongRef)));
output.refNoiseGain = 20*log10(abs(fft(averageRandRef)));
output.mesNoiseGain = 20*log10(abs(fft(averageRandMes)));
output.mesGain = 10*log10(mean(abs(fft(indivShortMes(:,1:3),fftlSeg*4)).^2,2));
output.mesNlGain = 10*log10(mean(abs(fft(indivShortMes(:,1:3)*dd - mean(indivShortMes(:,1:3)*dd,2),fftlSeg*4)).^2,2));
output.averageLongRef = averageLongRef;
output.averageRandRef = averageRandRef;
output.averageLongMes = averageLongMes;
output.averageRandMes = averageRandMes;
output.indivShortRef = indivShortRef;
output.indivShortMes = indivShortMes;
output.timeAxisLong = baseLongRespIdx(:)/fs;
output.timeAxisShort = baseShortRespIdx(:)/fs;
end