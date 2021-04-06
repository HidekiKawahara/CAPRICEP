function output = generateTestFMSignal(baseStr, fo, sigType, combId, fmCent, xin)
% output = generateTestFMSignal(baseStr, fo, sigType, combId, fmCent, xin)
%   Argument

startTic = tic;
output = struct;
fs = baseStr.fs;
xTSPSel = baseStr.xTSPSel;
fftl = size(xTSPSel, 1);
fftlSeg = fftl / 16;
orthSeq = zeros(length(xin), 3);
%--- generate orthogonal sequences
bMat = [1 1 1 1;1 -1 1 -1;1 1 -1 -1]';
c = baseStr.cOrdered;
baseIdx = (1:fftl)';
headLoc = 0;
capriTmp(:, 1) = xTSPSel(:, c(combId, 1));
capriTmp(:, 2) = xTSPSel(:, c(combId, 2));
capriTmp(:, 3) = xTSPSel(:, c(combId, 3));
nSegment = floor((length(xin) - fftl - fftlSeg) / (fftlSeg));
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
y = fftfilt(sixtermCos(fftlSeg/2), tmp);
yTrim = y(fftl * 2 + (1:fftl));
y = y / std(yTrim);
output.modulationSignal = y;
%--- make FM test signal
fi = 2 .^ (log2(fo) + y * fmCent / 1200);
output.fundamentalFreq = fi;
foBase = sin(cumsum(2*pi*fi/fs));
switch sigType
    case 'SINE'
        testSignal = foBase;
    case 'SINES'
        testSignal = foBase;
        for ii = 2:20
            testSignal = testSignal + sin(cumsum(ii*2*pi*fi/fs));
        end
    case 'MFND'
        testSignal = foBase * 0;
        for ii = 2:20
            testSignal = testSignal + sin(cumsum(ii*2*pi*fi/fs));
        end
    case 'MFNDH'
        testSignal = foBase * 0;
        for ii = 8:20
            testSignal = testSignal + sin(cumsum(ii*2*pi*fi/fs));
        end
end
w = hanning(4410);
testSignal(1:2205) = testSignal(1:2205) .* w(1:2205);
testSignal((1:2205) + end - 2205) = testSignal((1:2205) + end - 2205) ...
    .* w((1:2205) + end - 2205);
output.testSignal = testSignal;
output.xTSPSel = xTSPSel;
output.elapsedTime = toc(startTic);
end

