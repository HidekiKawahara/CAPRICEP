function output = fmTransfTestGaussN(fo, fMcent, outDirName, foWrapper)
%   Response measurement of pitch extractors to FM test signal
%   using narrow Gaussian pulse　(fmTransfTestGauss is obsolate)
%      output = fmTransfTestGaussN(fo, fMcent, outDirName, foWrapper)
%
%  Argument
%      fo         : fundamental frequency of the test signal (Hz)
%      fMcent     : frequency modulation depth (musical cent)
%      outDirname : path to the file output directory
%      foWrapper  : pointer to wrapper function of the pitch extractor
%                   for example "@pitchNCF" makes this program to call
%                   pitchNCF.m 
%
%  Output
%      output : structure with the following fields
%        foOut          : output of the tested pitch extractor
%        testCondStr    : structure consisting of the test conditions
%        testSignalData : structure consisting of used CAPRICEP set
%        testSigOut     : structure consisting of all internal data and
%                         test data
%        xa             : test signal: simulated vowel /a/ with FM
%        elapsedTime    : elapsed time of whole test (s)


% LICENSE: refer to LICENSE in this folder

startTic = tic;
testDuration = 20; % in second
testSignalData = load("testSigBase400ms.mat");
fs = testSignalData.fs;
tt = (0:1/fs:testDuration)' - testDuration;
carrierTypeList = {'SINE';'SINES';'MFND';'MFNDH'};
tmp = load("vaTmplt.mat");
a1Lpc = tmp.vaTmplt.a1Lpc;
combiID = 1;
carrierID = 2;
carrierType = carrierTypeList{carrierID};
testSigOut = zgenerateTestFMSignalG(testSignalData, fo, ...
    carrierType, combiID, fMcent, length(tt));
x = testSigOut.testSignal;
trueF0 = testSigOut.fundamentalFreq;
xa = filter(1, a1Lpc, x);
foOut = foWrapper(xa, fs);
testCondStr.outDirName = outDirName;
testCondStr.fo = fo;
testCondStr.fs = fs;
testCondStr.fMcent = fMcent;
testCondStr.carrierType = carrierType;
testCondStr.titleStr = foOut.titleStr;
testCondStr.filePrefix = foOut.filePrefix;
zdisplayFmTfncTest(xa, fs, trueF0, foOut.fo, foOut.tt, testCondStr, testSignalData);
output.foOut = foOut;
output.testCondStr = testCondStr;
output.testSignalData = testSignalData;
output.testSigOut = testSigOut;
output.xa = xa;
output.elaspedTime = toc(startTic);
end

%--- private functions 
function output = zgenerateTestFMSignalG(baseStr, fo, sigType, combId, fmCent, nSignal)
% output = generateTestFMSignalG(baseStr, fo, sigType, combId, fmCent, nSignal)
%   Argument

startTic = tic;
output = struct;
fs = baseStr.fs;
xTSPSel = baseStr.xTSPSel;
fftl = size(xTSPSel, 1);
fftlSeg = fftl / 16;
orthSeq = zeros(nSignal, 3);
%--- generate orthogonal sequences
bMat = [1 1 1 1;1 -1 1 -1;1 1 -1 -1]';
c = baseStr.cOrdered;
baseIdx = (1:fftl)';
headLoc = 0;
capriTmp(:, 1) = xTSPSel(:, c(combId, 1));
capriTmp(:, 2) = xTSPSel(:, c(combId, 2));
capriTmp(:, 3) = xTSPSel(:, c(combId, 3));
nSegment = floor((nSignal - fftl - fftlSeg) / (fftlSeg));
for ii = 1:nSegment
    kk = rem(ii - 1, 4) + 1;
    orthSeq(headLoc + baseIdx, :) = orthSeq(headLoc + baseIdx, :) ...
        +[capriTmp(:, 1), bMat(kk, 2) * capriTmp(:, 2), bMat(kk, 3) * capriTmp(:, 3)];
    headLoc = headLoc + fftlSeg;
end
output.orthSeq = orthSeq;
output.combination = c(combId, :);
%--- make filtered modulation signal
tmp = sum(orthSeq, 2);
fftlNgauss = fftl / 64;
ttw = ((1:fftlNgauss)' - fftlNgauss/2)/(fftlNgauss/2);
smoother = exp(-(6*ttw).^2); % Gaussian smoother
y = fftfilt(smoother, tmp);
yTrim = y(fftl + (1:fftl));
y = y / std(yTrim);
output.modulationSignal = y;
%--- make FM test signal
fi = 2 .^ (log2(fo) + y * fmCent / 1200);
output.fundamentalFreq = fi;
foBase = sin(cumsum(2*pi*fi/fs));
switch sigType
    case 'SINE'
        disp('SINE is disabled')
        return;
    case 'SINES'
        testSignal = foBase;
        for ii = 2:5000/fo % old 20
            testSignal = testSignal + sin(cumsum(ii*2*pi*fi/fs));
        end
    case 'MFND'
        testSignal = foBase * 0;
        for ii = 2:20
            testSignal = testSignal + sin(cumsum(ii*2*pi*fi/fs));
        end
    case 'MFNDH'
        disp('MFNDH is disabled')
        return;
end
w = hanning(4410);
testSignal(1:2205) = testSignal(1:2205) .* w(1:2205);
testSignal((1:2205) + end - 2205) = testSignal((1:2205) + end - 2205) ...
    .* w((1:2205) + end - 2205);
output.testSignal = testSignal;
output.xTSPSel = xTSPSel;
output.smoother = smoother;
output.fxs = (0:fftlSeg-1)/fftlSeg*fs;
output.fftlSeg = fftlSeg;
output.elapsedTime = toc(startTic);
end

function output = zdisplayFmTfncTest(x, fs, trueF0, estF0, estTime, testCondStr, testSignalData)
%%
startTic = tic;
titleStr = testCondStr.titleStr;
fo = testCondStr.fo;
fMcent = testCondStr.fMcent;
outDirName = testCondStr.outDirName;
filePrefix = testCondStr.filePrefix;
carrierType = testCondStr.carrierType;
tputG120 = zfoFixAndAligner(x, fs, fo, trueF0, estF0, estTime);
utputG120 = zfoResponseAnalysis(testSignalData, tputG120.rawF0i, tputG120.trueF0, fo, fs);
figure;
set(gcf, 'position', [681   700   671   277])
semilogx(utputG120.fx, utputG120.mesGain-utputG120.refGain, 'k-', 'linewidth', 2);
hold all;
semilogx(utputG120.fx, utputG120.mesNoiseGain-utputG120.refGain, 'k:', 'linewidth', 2);
semilogx(utputG120.fx, utputG120.mesNlGain-utputG120.refGain, 'k-.', 'linewidth', 2);
grid on;
set(gca, 'fontsize', 15, 'linewidth', 2)
axis([1 fo/2 -100 10])
title(titleStr + " fo:" + num2str(fo, '%7.1f') + " Hz, FMdpth:" ...
     + num2str(fMcent, '%7.1f') + " cent, narrow Gauss:");
xlabel('frequency (Hz)')
ylabel('gain (dB)')
legend('LTI', 'TV-rand', 'non-LTI', 'location', 'southeast', 'fontsize', 15);
fNameOut = filePrefix + datestr(now,30) + "NG";
drawnow;
print("-depsc", outDirName + "/" + fNameOut + ".eps");
save(outDirName + "/" + fNameOut + ".mat", "testCondStr");
disp(titleStr + "  at:" + datestr(now));
output.elapsedTime = toc(startTic);
end

function output = zfoResponseAnalysis(testSignalData, fopi, trueF0, fo, fs)

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

function output = zfoFixAndAligner(x, fs, fo, trueF0Org, rawF0, frameTime)
timeMargin = round(fs*0.4);
rawF0 = double(rawF0);
rawF0i = interp1(frameTime, rawF0, ...
    (timeMargin:length(x)-timeMargin)/fs,"linear","extrap");
rawF0i = rawF0i(:);
baseIdx = (1:length(rawF0i))';
[r, lags] = xcorr(rawF0i-fo, trueF0Org(baseIdx)-fo);
[~, maxLag] = max(r);
nData = length(trueF0Org);
trueF0 = trueF0Org(max(1, min(nData, baseIdx-lags(maxLag))));
timeAxis = (baseIdx-lags(maxLag))/fs;
output.rawF0i = rawF0i;
output.trueF0 = trueF0;
output.timeAxis = timeAxis(:);
end
