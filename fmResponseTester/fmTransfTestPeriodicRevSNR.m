function output = fmTransfTestPeriodicRevSNR(targetFo, fMcent, foWrapper, noisetype, snr)
%   Response measurement of pitch extractors to FM test signal
%   using narrow Gaussian pulse　(fmTransfTestGauss is obsolate)
%   with periodic modification with duration extension and periodicity fix
%      output = fmTransfTestPeriodicRev(fo, fMcent, outDirName, foWrapper)
%
%  Argument
%      targetFo   : fundamental frequency of the test signal (Hz)
%      fMcent     : frequency modulation depth (musical cent)
%      foWrapper  : pointer to wrapper function of the pitch extractor
%                   for example "@pitchNCF" makes this program to call
%                   pitchNCF.m 
%      noisetype  : noise type: default pink noise
%      snr        : signal to noise ratio (more than 1000 no noise)
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
testSigOut = zgenerateTestFMSignalG(testSignalData, targetFo, ...
    carrierType, combiID, fMcent, length(tt));
x = testSigOut.testSignal;
%trueF0 = testSigOut.fundamentalFreq;
xa = filter(1, a1Lpc, x);
switch noisetype
    case 'pink'
        if snr < 1000
            snrLvel = 10^(snr/20);
            xNoise = pinknoise(length(xa), 1);
            xa = xa + xNoise/std(xNoise)*std(xa)/snrLvel;
        end
end
ticExtract = tic;
foOut = foWrapper(xa, fs);
elapsedTimeExtraction = toc(ticExtract);

foMes = foOut.fo;
foMes(abs(foMes-targetFo) > targetFo/10) = targetFo;
txFrame = foOut.tt;
foRef = testSigOut.fundamentalFreq;
txSig = (1:length(foRef))'/fs;
foMesi = interp1(txFrame, foMes, txSig, "linear","extrap");
errorInd = 0;
if sum(abs(foMes - targetFo)) < 1
    foMesi = testSigOut.fundamentalFreq(:) + randn(length(foMesi) ,1);
    errorInd = 1;
end
fftl = testSignalData.fftl;
fftlSeg = testSignalData.fftlSeg;
xTsp = testSignalData.xTSPSel;
combination = testSigOut.combination;
xTspSel = xTsp(:, combination);
foMesiCent = 1200*log2(foMesi/targetFo);
foMesiCent(isnan(foMesiCent)) = 0;
foRefCent = 1200*log2(foRef/targetFo);
foRefCent(isnan(foRefCent)) = 0;
foRefX = fftfilt(xTspSel(end:-1:1,:), [foRefCent;zeros(fftl,1)]);
foRefX = foRefX(fftl/2 + (1:length(foRefCent)), :);
foMesX = fftfilt(xTspSel(end:-1:1,:), [foMesiCent;zeros(fftl,1)]);
foMesX = foMesX(fftl/2 + (1:length(foRefCent)), :);
lPeriod = fftlSeg*4;
nSeg = floor(length(foRef)/lPeriod);
segDiff = zeros(nSeg-1, 1);
baseIdx = 1:lPeriod;
for ii = 1:nSeg-1
    segDiff(ii) = std(foRefCent((ii-1)*lPeriod+baseIdx) ...
        - foRefCent(ii*lPeriod+baseIdx));
end
safeTh = -150;
segId = 1:nSeg-1;
goodSeg = segId(20*log10(segDiff) < safeTh);
%goodSeg = [goodSeg goodSeg(end)+1];
% use signal as a whole
nGoodSeg = length(goodSeg);
sRespInd = zeros(fftlSeg, nGoodSeg,3);
sRefInd = zeros(lPeriod, nGoodSeg,3);
lRefInd = zeros(lPeriod, nGoodSeg);
lMesInd = zeros(lPeriod, nGoodSeg);
lRespInd = zeros(lPeriod, nGoodSeg);
headShaper = (1:lPeriod)'/lPeriod;
postShaper = (lPeriod:-1:1)'/lPeriod;
baseSegIdx = 1:fftlSeg;
bmat = [1 1 1 1;1 -1 1 -1;1 1 -1 -1];
for ii = 1:nGoodSeg
    idg = goodSeg(ii);
    foRefTmp = foRefX((idg-1)*lPeriod+baseIdx,:) .* headShaper ...
        + foRefX(idg*lPeriod+baseIdx,:) .* postShaper;
    foMesTmp = foMesX((idg-1)*lPeriod+baseIdx,:) .* headShaper ...
        + foMesX(idg*lPeriod+baseIdx,:) .* postShaper;
    ref1 = zeros(lPeriod, 1);
    ref2 = zeros(lPeriod, 1);
    ref3 = zeros(lPeriod, 1);
    mes1 = zeros(lPeriod, 1);
    mes2 = zeros(lPeriod, 1);
    mes3 = zeros(lPeriod, 1);
    for jj = 1:4
        ref1 = ref1 + foRefTmp(circshift(baseIdx, -(jj-1)*fftlSeg),1);
        ref2 = ref2 + bmat(2,jj)*foRefTmp(circshift(baseIdx, -(jj-1)*fftlSeg),2);
        ref3 = ref3 + bmat(3,jj)*foRefTmp(circshift(baseIdx, -(jj-1)*fftlSeg),3);
        mes1 = mes1 + foMesTmp(circshift(baseIdx, -(jj-1)*fftlSeg),1);
        mes2 = mes2 + bmat(2,jj)*foMesTmp(circshift(baseIdx, -(jj-1)*fftlSeg),2);
        mes3 = mes3 + bmat(3,jj)*foMesTmp(circshift(baseIdx, -(jj-1)*fftlSeg),3);
    end
    sRespInd(:, ii, 1) = fft(mes1(baseSegIdx)) ./ fft(ref1(baseSegIdx));
    sRespInd(:, ii, 2) = fft(mes2(baseSegIdx)) ./ fft(ref2(baseSegIdx));
    sRespInd(:, ii, 3) = fft(mes3(baseSegIdx)) ./ fft(ref3(baseSegIdx));
    sRefInd(:, ii, 1) = fft(ref1);
    sRefInd(:, ii, 2) = fft(ref2);
    sRefInd(:, ii, 3) = fft(ref3);
    lRefF = fft(ref1+ref2+2*ref3)/4;
    lMesF = fft(mes1+mes2+2*mes3)/4;
    lRespInd(:, ii) = lMesF ./ lRefF;
    lRefInd(:, ii) = lRefF;
    lMesInd(:, ii) = lMesF;
end
varNL = zeros(fftlSeg, 1);
meanNL = zeros(fftlSeg, 1);
for ii = 1:3
    tmpResp = sRespInd(baseSegIdx,:,ii);
    %semilogx(fxSeg, 20*log10(abs(tmpResp)));grid on;        
    %hold all
    tmpAvResp = mean(tmpResp,2);
    meanNL = meanNL + tmpAvResp;
end
meanNL = meanNL/3;
meanSq = zeros(fftlSeg, 3);
for ii = 1:3
    tmpResp = sRespInd(baseSegIdx,:,ii);
    %semilogx(fxSeg, 20*log10(abs(tmpResp)));grid on;        
    %hold all
    tmpAvResp = mean(tmpResp,2);
    meanSq(:, ii) = tmpAvResp;
    varNL = varNL + abs(tmpAvResp - meanNL) .^2 ;
end
varNL = varNL /2;
varTV = zeros(fftlSeg, 3);
averageLResp = mean(lRespInd,2);
varTVL = zeros(lPeriod, 1);
for ii = 1:nGoodSeg
    tmpResp = sRespInd(baseSegIdx,ii,1);
    varTV(:, 1) = varTV(baseSegIdx, 1) + abs(meanSq(:, 1) - tmpResp) .^2;
    tmpResp = sRespInd(baseSegIdx,ii,2);
    varTV(:, 2) = varTV(:, 2) + abs(meanSq(:, 2) - tmpResp) .^2;
    tmpResp = sRespInd(baseSegIdx,ii,1);
    varTV(:, 3) = varTV(:, 3) + abs(meanSq(:, 3) - tmpResp) .^2;
    varTVL = varTVL + abs(lRespInd(:,ii)-averageLResp) .^2;
end
varTV = varTV/(nGoodSeg-1);
varTVL = varTVL/(nGoodSeg-1);
avVarTV = mean(varTV,2);

%testCondStr.outDirName = outDirName;
testCondStr.targetFo = targetFo;
testCondStr.fs = fs;
testCondStr.fMcent = fMcent;
testCondStr.carrierType = carrierType;
testCondStr.titleStr = foOut.titleStr;
testCondStr.filePrefix = foOut.filePrefix;
output.foOut = foOut;
output.testCondStr = testCondStr;
output.testSignalData = testSignalData;
output.testSigOut = testSigOut;

output.xa = xa;
output.sRefInd = sRefInd;
output.lRefInd = lRefInd;
output.lMesInd = lMesInd;
output.lRespInd = lRespInd;
output.averageResponse = abs(meanNL);
output.varNLTI = varNL;
output.varTV = avVarTV / 0.8164965^2; % compensation for weighting
output.varTVL = varTVL / 0.8164965^2; % compensation for weighting
output.fxSeg = (0:fftlSeg-1)'/fftlSeg*fs;
output.fxlResp = (0:lPeriod-1)'/lPeriod*fs;
%---- figure of merit
tmpVarNL = varNL;
tmpVarNL(1) = tmpVarNL(2); % ignore bias effect
varNLTIL = interp1(output.fxSeg,tmpVarNL,output.fxlResp);
totalErr = varNLTIL + output.varTVL;
snr = 20*log10(abs(mean(output.lRespInd,2)))-10*log10(totalErr);
snr(isnan(snr)) = 0;
snr(1) = snr(2);
fxIntLim = min(output.fxlResp(snr < 0));
pwspec = abs(mean(output.lRespInd,2)) .^2;
bww = sqrt(sum(pwspec(output.fxlResp<fxIntLim) .* output.fxlResp(output.fxlResp<fxIntLim).^2) ...
    / sum(pwspec(output.fxlResp<fxIntLim)));
effSnr = 10*log10(sum(pwspec(output.fxlResp<bww))) ...
    -10*log10(sum(totalErr(output.fxlResp<bww)));
output.minFqLimit = fxIntLim;
output.bandWidth = bww;
output.respSNR = effSnr;
output.avTotalErr = 10*log10(mean(totalErr(output.fxlResp<bww)));
%output.randomErr = 10*log10(mean(output.varTVL(output.fxlResp<bww)));
%-----
output.segInfo = 20*log10(segDiff);
output.elapsedTimeExtraction = elapsedTimeExtraction;
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
