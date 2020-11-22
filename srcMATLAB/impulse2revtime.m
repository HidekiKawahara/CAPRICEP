function output = impulse2revtime(specStr)
%  impulse response to reverberation time converter
%  output = impulse2revtime(specStr)
%
%  Argument
%    specStr      : structure variable of spectrum analysis
%    
% Output
%    output    : structure variable with the following fields
%           specStr          : copy of the analysis results
%           reportStr        : copy of the report results
%
%
% Licence 
% Copyright 2020 Hideki Kawahara
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

strtTic = tic;
fs = specStr.fs;
aweightfilter = weightingFilter("A-weighting", fs);
strtTic2 = tic;
xa = aweightfilter(specStr.averageLongResponse);
xn = aweightfilter(specStr.averageLongRTVcomp);
strtTic3 = tic;
tx = specStr.timeAxisLong;
[~, idxHlf] = min(abs(tx - tx(end)/2));
revPower = cumsum(abs(xa).^2, "reverse");
revRTVPower = cumsum(abs(xn).^2, "reverse");
dBmismatch = 10*log10(revRTVPower(idxHlf) / revPower(idxHlf));
revRTVPower = revRTVPower / revRTVPower(idxHlf) * revPower(idxHlf);
fixedCumResp = 10*log10(abs(revPower - revRTVPower));
noiseResp = 10*log10(revRTVPower);
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
        tTrim = tx(best10ms:max(best10ms + 441, bestIndex10));
        H = [tTrim, ones(length(tTrim), 1)];
        regCoeff = (H'* H) \ (H'*fixedCumResp(best10ms:max(best10ms + 441, bestIndex10)));
        reverberationTimeDistribution(iterId) = -60/regCoeff(1);
    end
end
tbound = 0.015;
directRevDiffDB = 10*log10(sum(xa(tx < tbound).^2))-10*log10(sum(xa(tx > tbound).^2));
magFactor = 10^(directRevDiffDB/20);
octStr = designOctabeFilger(specStr);
octRevStr = octaveReverb(specStr, octStr);
%%
output.octStr = octStr;
output.octRevStr = octRevStr;
output.magFactor = magFactor;
output.dBmismatch = dBmismatch;
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
filtOutRV = octStr.filtOutRV;
nte = round(fs*specStr.tResponse*4/1000);
tx = ((1:nte)' - octStr.maxNtc)/fs;
[~, idxHlf] = min(abs(tx - tx(end)/2));
initMarginList = 0.020:0.001:0.035;
biasExt = 0:0.3:5;
revTimeList = zeros(length(2:3:length(fcList)),length(initMarginList)*length(biasExt));
for ii = 1:length(2:3:length(fcList))
    revPower = cumsum(abs(filtOut(1:nte, ii)).^2, "reverse");
    revRVPower = cumsum(abs(filtOutRV(1:nte, ii)).^2, "reverse")/2;
    revRVPower = revRVPower / revRVPower(idxHlf) * revPower(idxHlf);
    fixedCumResp = 10*log10(abs(revPower - revRVPower));
    noiseResp = 10*log10(revRVPower);
    iterId = 0;
    for jj = 1:length(initMarginList)
        initialMargin = initMarginList(jj);
        for kk = 1:length(biasExt)
            biasL = biasExt(kk);
            iterId = iterId + 1;
            [~, bestIndex10] = min(abs(fixedCumResp(tx < tx(end)/2) - noiseResp(tx < tx(end)/2) + biasL));
            [~, best35ms] = min(abs(tx - initialMargin));
            %tTrim = tx(best35ms:bestIndex10);
            tTrim = tx(best35ms:max(best35ms + 441, bestIndex10));
            H = [tTrim, ones(length(tTrim), 1)];
            regCoeff = (H'*H) \ (H'*fixedCumResp(best35ms:max(best35ms + 441, bestIndex10)));
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
%fxw = (0:fftlw-1)/fftlw*fs;
octFilter = zeros(fftlw, length(fcList));
maxNtc = 0;
for ii = 1:length(fcList)
    fc = fcList(ii);
    ntc = round(3*fs/fc);
    w = blackman(ntc*2+1);
    tmt = ((1:2*ntc+1)' - ntc)/fs;
    idx = fftlw/2 + (-ntc:ntc)';
    octFilter(idx, ii) = w/sum(w) .* exp(2*1i*pi*fc*tmt);
    if maxNtc < ntc
        maxNtc = ntc;
    end
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
tmpOcFilter = circshift(tmpOcFilter,maxNtc);
octFilterFx = fft(tmpOcFilter);
xTmp = zeros(fftlw,1);
xRVTmp = zeros(fftlw,1);
xTmp(1:length(specStr.averageLongResponse)) = specStr.averageLongResponse;
xRVTmp(1:length(specStr.averageLongResponse)) = specStr.averageLongRTVcomp;
xFx = fft(xTmp);
xRVFx = fft(xRVTmp);
filtOut = zeros(fftlw, length(2:3:length(fcList)));
filtOutRV = zeros(fftlw, length(2:3:length(fcList)));
for ii = 1:length(2:3:length(fcList))
    filtOut(:, ii) = abs(ifft(xFx .* octFilterFx(:, ii)));
    filtOutRV(:, ii) = abs(ifft(xRVFx .* octFilterFx(:, ii)));
end
octStr.filtOut = filtOut;
octStr.filtOutRV = filtOutRV;
octStr.maxNtc = maxNtc;
octStr.fcList = fcList;
octStr.fcOctFilter = fcOctFilter;
end