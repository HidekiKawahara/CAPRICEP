function simulationOut = acousticSimulator(xTestPink, option)

nto = option.nto;
fftl = option.fftl;
%nRepetition = option.nRepetition;
fs = option.fs;
snr = option.snr;
testMode = option.testMode;

%---- pink noise shaper
fLow = 40;
fx = (0:fftl - 1)' / fftl * fs;
fx(fx > fs/2) = fx(fx > fs/2) - fs;
g = 1.0 ./ abs(fx / fLow);
g(g > 1) = 1;
rawAutoCorr = real(ifft(g));
rawAutoCorr = rawAutoCorr / rawAutoCorr(1);
np = 50;
[pinkLPC, ~] = levinson(rawAutoCorr, np);

% ---- test body
switch testMode
    case 'SNR'
        sigStd = std(xTestPink(round(length(xTestPink) /2) + (-nto:nto)));
        pingNoise = randn(length(xTestPink), 1);
        pingNoise = filter(1, pinkLPC, pingNoise);
        
        simulationOut.y = xTestPink ...
            + pingNoise / std(pingNoise) * sigStd * 10 ^ (-snr / 20);
        simulationOut.xTestPink = xTestPink;
        simulationOut.noise = pingNoise / std(pingNoise) * sigStd * 10 ^ (-snr / 20);
    case 'nonLinear'
        xTestPink = xTestPink * 10 ^ (option.level/20);
        sigStd = std(xTestPink(round(length(xTestPink) /2) + (-nto:nto)));
        pingNoise = randn(length(xTestPink), 1);
        pingNoise = filter(1, pinkLPC, pingNoise);
        xOrg = filter(pinkLPC, 1, xTestPink);
        nlXorg = expNonlinear(xOrg, 0.3);
        
        simulationOut.y = filter(1, pinkLPC, nlXorg) ...
            + pingNoise / std(pingNoise) * sigStd * 10 ^ (-snr / 20);
        simulationOut.xTestPink = xTestPink;
        simulationOut.noise = pingNoise / std(pingNoise) * sigStd * 10 ^ (-snr / 20);
        simulationOut.y = simulationOut.y - mean(simulationOut.y(round(length(xTestPink) /2) + (-nto:nto)));
        lSeg = xOrg(round(length(xTestPink) /2) + (-nto:nto));
        nlSeg = nlXorg(round(length(xTestPink) /2) + (-nto:nto));
        nlSegN = nlSeg - mean(nlSeg);
        lSegN = lSeg - mean(lSeg);
        areg = sum(lSegN .* nlSegN)/sum(lSegN.*lSegN);
        ee = nlSegN - areg * lSegN;
        simulationOut.nlLeveldB = 10*log10(sum(ee .* ee)/sum(lSegN.*lSegN));
end
end