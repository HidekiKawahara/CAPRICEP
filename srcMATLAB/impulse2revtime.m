function output = impulse2revtime(specStr)

strtTic = tic;
fs = specStr.fs;
aweightfilter = weightingFilter("A-weighting", fs);
strtTic2 = tic;
xa = aweightfilter(specStr.averageShortResponse);
strtTic3 = tic;
tx = specStr.timeAxisShort;
revPower = cumsum(abs(xa).^2, "reverse");
xnStd = std(xa(tx > tx(end)/2));
simNoise = ones(length(tx), 1) * xnStd;
revNoiseSimPower = cumsum(abs(simNoise).^2, "reverse");
fixedCumResp = 10*log10(abs(revPower - revNoiseSimPower));
noiseResp = 10*log10(revNoiseSimPower);
initMarginList = 0.020:0.001:0.035;
biasExt = 0:0.3:5;
reverberationTimeDistribution = zeros(length(initMarginList)*length(biasExt), 1);
iterId = 0;
for jj = 1:length(initMarginList)
    initialMargin = initMarginList(jj);
    for kk = 1:length(biasExt)
        biasL = biasExt(kk);
        iterId = iterId + 1;
        [~, bestIndex10] = min(abs(fixedCumResp(tx < tx(end)/2) - noiseResp(tx < tx(end)/2) + biasL));
        [~, best10ms] = min(abs(tx - initialMargin));
        tTrim = tx(best10ms:bestIndex10);
        H = [tTrim, ones(length(tTrim), 1)];
        regCoeff = (H'*H) \ (H'*fixedCumResp(best10ms:bestIndex10));
        reverberationTimeDistribution(iterId) = -60/regCoeff(1);
    end
end
tbound = 0.015;
directRevDiffDB = 10*log10(sum(xa(tx < tbound).^2))-10*log10(sum(xa(tx > tbound).^2));
magFactor = 10^(directRevDiffDB/20);
octStr = designOctabeFilger(specStr);
octRevStr = octaveReverb(specStr, octStr);
output.octStr = octStr;
output.octRevStr = octRevStr;
output.magFactor = magFactor;
output.reverberationTimeDistribution = reverberationTimeDistribution;
output.reverberationTime = median(reverberationTimeDistribution);
output.totalElapsedTime = toc(strtTic);
output.elapsedTime1 = toc(strtTic2);
output.elapsedTime2 = toc(strtTic3);
end

function octRevStr = octaveReverb(specStr, octStr)
fs = specStr.fs;
fcList = octStr.fcList;
filtOut = octStr.filtOut;
nte = round(fs*specStr.tResponse/1000*4);
tx = (1:nte)'/fs;
initMarginList = 0.020:0.001:0.035;
biasExt = 0:0.3:5;
revTimeList = zeros(length(2:3:length(fcList)),length(initMarginList)*length(biasExt));
for ii = 1:length(2:3:length(fcList))
    xnStd = sqrt(mean(filtOut(round(nte/2):nte, ii) .^2));
    simNoise = ones(length(tx), 1) * xnStd;
    revPower = cumsum(abs(filtOut(1:nte, ii)).^2, "reverse");
    revNoiseSimPower = cumsum(abs(simNoise).^2, "reverse");
    fixedCumResp = 10*log10(abs(revPower - revNoiseSimPower));
    noiseResp = 10*log10(revNoiseSimPower);
    iterId = 0;
    for jj = 1:length(initMarginList)
        initialMargin = initMarginList(jj);
        for kk = 1:length(biasExt)
            biasL = biasExt(kk);
            iterId = iterId + 1;
            [~, bestIndex10] = min(abs(fixedCumResp(tx < tx(end)/2) - noiseResp(tx < tx(end)/2) + biasL));
            [~, best35ms] = min(abs(tx - initialMargin));
            tTrim = tx(best35ms:bestIndex10);
            H = [tTrim, ones(length(tTrim), 1)];
            regCoeff = (H'*H) \ (H'*fixedCumResp(best35ms:bestIndex10));
            revTimeList(ii, iterId) = -60/regCoeff(1);
        end
    end
end
octRevStr.revTimeList = revTimeList;
octRevStr.fcOctFilter = octStr.fcOctFilter;

end

function octStr = designOctabeFilger(specStr)
fs = specStr.fs;
fftlw = 2^ceil(log2(length(specStr.averageLongResponse)*2));
fcList = 50 * 2 .^(0:1/3:log2(fs/2/50));
fxw = (0:fftlw-1)/fftlw*fs;
octFilter = zeros(fftlw, length(fcList));
for ii = 1:length(fcList)
    fc = fcList(ii);
    ntc = round(3*fs/fc);
    w = blackman(ntc*2+1);
    tmt = ((1:2*ntc+1)' - ntc)/fs;
    idx = fftlw/2 + (-ntc:ntc)';
    octFilter(idx, ii) = w/sum(w) .* exp(2*1i*pi*fc*tmt);
end

octFilterSet = zeros(fftlw, length(2:3:length(fcList)));
fcOctFilter = zeros(length(2:3:length(fcList)), 1);
id = 0;
for ii = 2:3:length(fcList)
    id = id + 1;
    fcOctFilter(id) = fcList(ii);
    octFilterSet(:, id) = sum(octFilter(:, ii+(-1:1)),2);
end

tmpOcFilter = fftshift(octFilterSet,1);
octFilterFx = fft(tmpOcFilter);
xTmp = zeros(fftlw,1);
xTmp(1:length(specStr.averageLongResponse)) = specStr.averageLongResponse;
xFx = fft(xTmp);
filtOut = zeros(fftlw, length(2:3:length(fcList)));
for ii = 1:length(2:3:length(fcList))
    filtOut(:, ii) = abs(ifft(xFx .* octFilterFx(:, ii)));
end
octStr.filtOut = filtOut;
octStr.fcList = fcList;
octStr.fcOctFilter = fcOctFilter;
end