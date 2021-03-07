function output = capricepResponseAnalysis(analysisStr)
% Analyse the acquired sound data
% output = capricepResponseAnalysis(analysisStr)
% 
%
% Argument
%   analysisStr  : structure variable with the following fields
%                : Below shows an example from acquisition function
%   (necessary fields)
%              yRecorded: [241937×1 double] % recorded signal
%                     fs: 44100
%                pinkLPC: [1×51 double] % paramet for IIR pink noise 
%              tResponse: 100
%             outChannel: 'L-ch'
%            numChannels: 1    % number of channels assigned to A/D
%            nRepetition: 30
%                   lAeq: -30.6956  % calibrated sound pressure level in
%                                     A-weighted dB rel. 20 micro Pa
%   (The following fields are ignored)
%          simulationOut: [1×1 struct] % only appears 'simulator mode'
%                tspName: "setCAPRICEP100msRc.mat" % file name of the unit
%       calibrationConst: 0   % simulator did not used this
%                   fLow: 40  % low frequency limit for pink noise design
%              xTestPink: [197837×1 double] % used output test signal
%     testSignalDuration: 4.4861
%             yRecovered: [241937×1 double] % pulse commpressed signal
%       orthogonalSignal: [241937×4 double] % orthogonalized signal
%    orthogonalSignalRef: [241937×4 double] % for debug. tentative
%            elapsedTime: 0.3602  % total elapsed time (s)
%     elapsedTimeRecover: 0.1838  % elapsed time for post processing (s)
%
% Output
%  output    : structure variable with the following fields
%                  shortSpec: [65536×1 double]  % raw LTI spectrum
%                   longSpec: [65536×1 double]  % raw extended LIT spectrum
%         deviationPowerSpec: [65536×1 double]  % raw nonlinear TI pw spec
%            randomPowerSpec: [65536×1 double]  % raw random power spec.
%              frequencyAxis: [1×65536 double]  % frequency axis
%               prePowerSpec: [65536×1 double]  % raw background noise
%                                calculated using the preceding silence
%                         fs: 44100
%                    tspName: "setCAPRICEP100msRc.mat" % namr of TSP
%                  tResponse: 100
%                nRepetition: 30
%                    pinkLPC: [1×51 double]
%                 outChannel: 'L-ch'
%                numChannels: 1
%    initialAnalysisPosition: 107738            % for alignment later
%           rawShortResponse: [4410×4 double]   % set pf impulse responses
%       averageShortResponse: [4410×1 double]   % averaged impulse resp.
%              timeAxisShort: [4410×1 double]   % time axis for short ones
%               timeAxisLong: [17640×1 double]  % time axis for long ones
%           headMarginSample: 441               % preceding samples from 0
%        averageLongResponse: [17640×1 double]  % long impulse response
%           orthogonalSignal: [241937×4 double] % orthogonalized signals
%                  sumSignal: [241937×1 double] % for alignment
%                       lAeq: -30.6956
%                elapsedTime: 0.1631

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


startTic = tic;
% ----
y = analysisStr.yRecorded;
pinkLPC = analysisStr.pinkLPC;
fs = analysisStr.fs;
tResponse = analysisStr.tResponse;
switch fs
    case 44100
        switch tResponse
            case 100
                tspName = "setCAPRICEP100msRc.mat";
            case 200
                tspName = "setCAPRICEP200msRc.mat";
            case 400
                tspName = "setCAPRICEP400msRc.mat";
            case 800
                tspName = "setCAPRICEP800msRc.mat";
            otherwise
                output = [];
                disp("Available response time is 100, 200, and 400 ms, 800 ms, for the time being");
                return;
        end
        tmp = load(tspName);
    case 176400
        switch tResponse
            case 100
                tspName = "setCAPRICEP400msRc.mat";
            case 200
                tspName = "setCAPRICEP800msRc.mat";
            case 400
                tspName = "setCAPRICEP1600msRc.mat";
            otherwise
                output = [];
                disp("Available response time is 100, 200, and 400 ms with 176400 Hz for the time being");
                return;
        end
        tmp = load(tspName);
    case 192000
        switch tResponse
            case 100
                tspName = "setCAPRICEP400msRc.mat";
            case 200
                tspName = "setCAPRICEP800msRc.mat";
            case 400
                tspName = "setCAPRICEP1600msRc.mat";
            otherwise
                output = [];
                disp("Available response time is 100, 200, and 400 ms with 192000 Hz for the time being");
                return;
        end
        tmp = load(tspName);
    otherwise
        output = [];
        disp("Availabel sampling frequency is 44100 Hz, for the time being.")
        return;
end
xTSPSel = tmp.xTSPSel;
B4 = [1 1 1 1 1 1 1 1;1 -1 1 -1 1 -1 1 -1;1 1 -1 -1 1 1 -1 -1;1 1 1 1 -1 -1 -1 -1];
nto = round(tResponse / 1000 * fs);
outChannel = analysisStr.outChannel;
%numChannels = analysisStr.numChannels;
inChannel = size(y,2);
if isfield(analysisStr, 'selectedChannels')
    numChannels = length(analysisStr.selectedChannels);
else
    numChannels = inChannel;
end
%% 
%------ equalization of pink noise shaping
%
yr = filter(pinkLPC, 1, y);
fftl = size(xTSPSel, 1);
nRepetition = analysisStr.nRepetition;
nItem = 0;
simSilence = zeros(2 * length(yr), inChannel);
while nItem < length(yr) + 2 * fs
    endp = round(fs/2+rand(1,1)*fs/3);
%    tmpSig = [];
    tmpSig = [yr(200:endp, :);yr(endp-1:-1:201, :)];
    simSilence(nItem + (1:length(tmpSig)), :) = tmpSig;
    nItem = nItem + length(tmpSig);
end
simSilence = simSilence(fs + (1:length(yr)), :);
simSilence = simSilence - mean(simSilence);
%%
%------ pulse recovery (compression)
%
selectInvIndex = fftl / 2 + (3*nto:-1:-3*nto);
xTSPSelIv = xTSPSel(selectInvIndex, :);
compSignal = zeros(length(yr), inChannel, 4);
baseIdx = (1:length(yr));

compSignal(:, :, 1) = fftfilt(xTSPSelIv(:, 1), yr);
compSignal(:, :, 2) = fftfilt(xTSPSelIv(:, 2), yr);
compSignal(:, :, 3) = fftfilt(xTSPSelIv(:, 3), yr);
compSignal(:, :, 4) = fftfilt(xTSPSelIv(:, 4), yr);
compSilence = fftfilt(xTSPSelIv(:, 1), simSilence);
%
%------- orthogonalization
%
orthogonalSignal = zeros(length(yr), inChannel, 4);
orthogonalSilence = zeros(length(yr), inChannel);
for ii = 1:8
    tmpIdx = min(length(yr), baseIdx + (ii - 1) * nto);
    for jj = 1:4
        orthogonalSignal(:, :, jj) = orthogonalSignal(:, :, jj) ...
            + B4(jj, ii) * compSignal(tmpIdx, :, jj);
    end
    orthogonalSilence = orthogonalSilence ...
        + compSilence(tmpIdx, :);
end
%dmy = 1;
orthogonalSignal = orthogonalSignal / 8;
orthogonalSilence = orthogonalSilence / sqrt(8);
%
%% ----- analysis location alignment
%
tmpIdx = (1:length(yr))';
sumSignal = fftfilt(hanning(45), abs(sum(orthogonalSignal(:, 1, 1:3),3)));
maxSignal = max(abs(sumSignal));
sigPeaksRaw = tmpIdx(abs(sumSignal) > 0.95 * maxSignal & sumSignal > sumSignal([1 1:end-1]) & sumSignal >= sumSignal([2:end end])); 
sigPeaks = sigPeaksRaw;
peakId = 0;
ii = 1;
idxBase = 1:length(sigPeaksRaw);
while ii <= length(sigPeaksRaw)
    if sum(abs(sigPeaksRaw(ii) - sigPeaksRaw) < 441) > 1
        peakId = peakId + 1;
        sigPeaksTmp = sigPeaksRaw(abs(sigPeaksRaw(ii) - sigPeaksRaw) < 441);
        idxTmp = idxBase(abs(sigPeaksRaw(ii) - sigPeaksRaw) < 441);
        [~, idxPeak] = max(sumSignal(abs(sigPeaksRaw(ii) - sigPeaksRaw) < 441));
        sigPeaks(peakId) = sigPeaksTmp(idxPeak);
        ii = max(idxTmp) + 1;
    elseif sum(abs(sigPeaksRaw(ii) - sigPeaksRaw) < 441) == 1
        peakId = peakId + 1;
        sigPeaks(peakId) = sigPeaksRaw(ii);
        ii = ii + 1;
    end
end
sigPeaks = sigPeaks(1:peakId) - 22;

%safeSigPeaks = sigPeaks(2:end-1); % needs some safety checker hear!
safeSigPeaksLong = sigPeaks(2:2:end-1);
%%   initial point refinmement

headMargin = round(fs * 0.02); % 20 ms preceding margin
averageShortResponse = zeros(nto, inChannel);
for ii = 1:length(safeSigPeaksLong)
    selIdx = (1:nto) - headMargin + safeSigPeaksLong(1) + (ii - 1) * 8 * nto;
    averageShortResponse = averageShortResponse + sum(orthogonalSignal(selIdx, :, 1:3), 3);
end
tmp = averageShortResponse / length(safeSigPeaksLong) / 3;
cumPower = cumsum(tmp.^2)/sum(tmp.^2);
tmpIdx = 1:length(tmp);
headLocation = min(tmpIdx(cumPower > 0.01));
safeSigPeaksLong = safeSigPeaksLong - (headMargin - headLocation);
%
%% ----- individual and average short responses
%
headMargin = round(fs * 0.01); % 10 ms preceding margin
averageShortResponse = zeros(nto, inChannel);
rawShortResponse = zeros(nto, inChannel, 4);
averageSilence = zeros(nto, inChannel);
for ii = 1:length(safeSigPeaksLong)
    selIdx = (1:nto) - headMargin + safeSigPeaksLong(1) + (ii - 1) * 8 * nto;
    rawShortResponse = rawShortResponse + orthogonalSignal(selIdx, :, :);
    averageShortResponse = averageShortResponse + sum(orthogonalSignal(selIdx, :, 1:3), 3);
    averageSilence = averageSilence + orthogonalSilence(selIdx, :);
end
averageShortResponse = averageShortResponse / length(safeSigPeaksLong) / 3;
rawShortResponse = rawShortResponse / length(safeSigPeaksLong);
averageSilence = averageSilence / sqrt(length(safeSigPeaksLong));
%
%% ----- averaged long response (four times)
%
respExtendedRaw = (orthogonalSignal(:, :, 1) + orthogonalSignal(:, :, 2)) / 4 + orthogonalSignal(:, :, 3) / 2;
averageLongResponse = zeros(4 * nto, inChannel);
averageLongRTVcomp = zeros(4 * nto, inChannel);
for ii = 1:length(safeSigPeaksLong)
    selIdx = (1:4 * nto) - headMargin + safeSigPeaksLong(1) +  (ii - 1) * 8 * nto;
    averageLongResponse = averageLongResponse + respExtendedRaw(selIdx, :);
    averageLongRTVcomp = averageLongRTVcomp + orthogonalSignal(selIdx, :, 4);
end
averageLongResponse = averageLongResponse / length(safeSigPeaksLong);
averageLongRTVcomp = averageLongRTVcomp / length(safeSigPeaksLong);

%% ----- spectrum analysis
tmpAverageResponse = zeros(nto, inChannel, 3);
tmpAverageResponse(:, :, 1) = averageShortResponse;
tmpAverageResponse(:, :, 2) = averageShortResponse;
tmpAverageResponse(:, :, 3) = averageShortResponse;
deviationResponse = rawShortResponse(1:nto, :, 1:3) - tmpAverageResponse;
%    - [averageShortResponse averageShortResponse averageShortResponse];
%for ii = 1:3
%    for jj = 1:3
%        deviationResponse(:, ii) = deviationResponse(:, ii) + averageLongResponse(jj * nto + (1:nto))/3;
%    end
%end
for ii = 1:3
    deviationResponse(:, :, ii) = deviationResponse(:, :, ii) - mean(deviationResponse(:, :, ii));
end
deviationSpec = fft(deviationResponse, fftl);
randomSpec = fft(rawShortResponse(:,:, 4), fftl);
output.shortSpec = fft(averageShortResponse, fftl);
output.longSpec = fft(averageLongResponse, fftl);
output.deviationPowerSpec = sum(abs(deviationSpec) .^2, 3)/2 * length(safeSigPeaksLong);
output.randomPowerSpec = abs(randomSpec) .^2 * 8 * length(safeSigPeaksLong) / 3;
output.frequencyAxis = (0:fftl-1)/fftl * fs;
output.prePowerSpec = abs(fft(averageSilence, fftl)) .^2 / 3;
%% calibration information recovery
if analysisStr.caliblationConst == 0
weightFilt = weightingFilter('A-weighting' ,fs);
yAweight = weightFilt(analysisStr.yRecorded);
caliblationConst = analysisStr.lAeq - 20*log10(std(yAweight(round(length(analysisStr.yRecorded) / 2) + (-nto:nto), :)));
else
    caliblationConst = analysisStr.caliblationConst;
end

%%
output.fs = fs;
output.tspName = tspName;
output.tResponse = tResponse;
output.nRepetition = nRepetition;
output.pinkLPC = pinkLPC;
output.outChannel = outChannel;
output.numChannels = numChannels;
output.initialAnalysisPosition = safeSigPeaksLong(1);
output.rawShortResponse = rawShortResponse;
output.averageShortResponse = averageShortResponse;
output.timeAxisShort = ((1:nto)' - headMargin) / fs;
output.timeAxisLong = ((1:4*nto)' - headMargin) / fs;
output.headMarginSample = headMargin;
output.averageLongResponse = averageLongResponse;
output.averageLongRTVcomp = averageLongRTVcomp;
output.orthogonalSignal = orthogonalSignal;
output.sumSignal = sumSignal;
output.selectedChannels = analysisStr.selectedChannels;
output.caliblationConst = caliblationConst;
output.lAeq = analysisStr.lAeq;
output.elapsedTime = toc(startTic);
end

