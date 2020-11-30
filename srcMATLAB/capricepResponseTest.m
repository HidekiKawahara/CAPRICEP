function output = capricepResponseTest(fs, tResponse, nRepetition, outChannel, inChannel, varargin)
% Generate test signal and execute acoustic measurement using the test
% signal
% 
% output = 
%   capricepResponseTest(fs, tResponse, nRepetition, outChannel, inChannel)
% output =
%   capricepResponseTest(fs, tResponse, nRepetition, outChannel, inChannel
%      testMode, option)
%
% Argument
%   fs             : sampling frequency (Hz)
%   tResponse      : impulse response length (ms): acceptable lengths are
%                      100, 200, 400, 800
%   nRepetition    : number of repetitions of unit-CAPRICEPs
%                    30 or more is recommended
%   outChannel     : output channel ID, 1 or 2. 1 is for L-channel
%   inChannel      : nunber of input channels
%   testMode       : test mode, 'acoustic_system' (default) or 'simulator'
%   option         : structure with the following fields for simulator
%      calibrationConst   : a constant to convert to sound pressure level
%      selectedChannels   : a set of channel IDs to record
%                  : following fields are internally set
%        option.fs = fs;
%        option.nto = nto;
%        option.fftl = fftl;
%        option.nRepetition = nRepetition;
%
% Output
%   output  : structure variable with the following fields
%           : Below shows an example output
%          simulationOut: [1×1 struct] % only appears 'simulator mode'
%                     fs: 44100
%                tspName: "setCAPRICEP100msRc.mat" % file name of the unit
%              tResponse: 100
%            nRepetition: 30
%       calibrationConst: 0   % simulator did not used this
%                   fLow: 40  % low frequency limit for pink noise design
%                pinkLPC: [1×51 double] % paramet for IIR pink noise 
%              xTestPink: [197837×1 double] % used output test signal
%             outChannel: 'L-ch'
%            numChannels: 1    % number of channels assigned to A/D
%     testSignalDuration: 4.4861
%              yRecorded: [241937×1 double] % recorded signal
%             yRecovered: [241937×1 double] % pulse commpressed signal
%       orthogonalSignal: [241937×4 double] % orthogonalized signal
%    orthogonalSignalRef: [241937×4 double] % for debug. tentative
%            elapsedTime: 0.3602  % total elapsed time (s)
%                   lAeq: -30.6956  % calibrated sound pressure level in
%                                     A-weighted dB rel. 20 micro Pa
%     elapsedTimeRecover: 0.1838  % elapsed time for post processing (s)

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


switch nargin
    case 5
        testMode = 'acoustic_system';
        selectedChannels = 1:inChannel;
        calibrationConst = 0;
    case 7
        testMode = varargin{1};
        option = varargin{2};
        calibrationConst = option.calibrationConst;
        if isfield(option, 'selectedChannels')
            selectedChannels = option.selectedChannels;
            calibrationConst = option.calibrationConst(selectedChannels);
        else
            selectedChannels = 1:inChannel;
        end
end
startTic = tic;
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
                disp("Available response time is 100, 200, 400 ms, and 800 ms for the time being");
                return;
        end
        tmp = load(tspName);
    otherwise
        output = [];
        disp("Availabel sampling frequency is 44100 Hz, for the time being.")
        return;
end
%---- test sequence ----
xTSPSel = tmp.xTSPSel;
B4 = [1 1 1 1 1 1 1 1;1 -1 1 -1 1 -1 1 -1;1 1 -1 -1 1 1 -1 -1;1 1 1 1 -1 -1 -1 -1];
fftl = size(xTSPSel, 1);
nto = round(tResponse / 1000 * fs);
xTest = zeros(fftl + nRepetition * nto + 1, 1);
tTestSignal = length(xTest) / fs;
for ii = 1:nRepetition
    tmpIdx = (ii - 1) * nto + (1:fftl)';
    iiCycle = rem(ii - 1, 8) + 1;
    xTest(tmpIdx) = xTest(tmpIdx) ...
        + B4(1, iiCycle) * xTSPSel(:, 1) ...
        + B4(2, iiCycle) * xTSPSel(:, 2) ...
        + B4(3, iiCycle) * xTSPSel(:, 3);
end
%---- pink noise shaping ----
fLow = 40;
fx = (0:fftl - 1)' / fftl * fs;
fx(fx > fs/2) = fx(fx > fs/2) - fs;
g = 1.0 ./ abs(fx / fLow);
g(g > 1) = 1;
rawAutoCorr = real(ifft(g));
rawAutoCorr = rawAutoCorr / rawAutoCorr(1);
np = 50;
[pinkLPC, ~] = levinson(rawAutoCorr, np);
xTestPinkOrg = filter(1, pinkLPC, xTest);
xTestPink = 0.8 * xTestPinkOrg / max(abs(xTestPinkOrg));
%---- test target system ----
switch testMode
    case 'acoustic_system'
        switch outChannel
            case 'L-ch'
                xOut = [xTestPink xTestPink * 0];
            case 'R-ch'
                xOut = [xTestPink * 0 xTestPink];
            otherwise
                disp("use L-ch or R-ch");
                output = [];
                return;
        end
        nBits = 24;
        playerObj = audioplayer(xOut, fs, nBits);
        if audiodevinfo(1,-1,fs,nBits,inChannel)
            recorderObj = audiorecorder(fs, nBits, inChannel);
        else
            disp("This device does not support;  fs:" + num2str(fs) ...
                + " (Hs)  nBits:" + num2str(nBits) + "  nChannel:" ...
                + num2str(inChannel));
            output = [];
            return;
        end
        numChannels = get(recorderObj, 'NumChannels');
        record(recorderObj);
        pause(1)
        play(playerObj);
        pause(tTestSignal + 1);
        y = getaudiodata(recorderObj);
        stop(recorderObj);
        switch inChannel
            case 1
                y = y(:, 1);
            otherwise
                if length(selectedChannels) == 1
                    inChannel = 1;
                    y = y(:, selectedChannels);
                else
                    y = y(:, selectedChannels);
                end
        end
    case 'simulator'
        numChannels = 1;
        option.fs = fs;
        option.nto = nto;
        option.fftl = fftl;
        option.nRepetition = nRepetition;
        simulationOut = acousticSimulator([zeros(fs, 1) ;xTestPinkOrg], option);
        y = simulationOut.y;
        output.simulationOut = simulationOut;
end
%------
weightFilt = weightingFilter('A-weighting' ,fs);
yAweight = weightFilt(y);
lAeq = 20*log10(std(yAweight(round(length(y) / 2) + (-nto:nto), :))) + calibrationConst;
middleTic = tic;
yr = filter(pinkLPC, 1, y);
%------
selectInvIndex = fftl / 2 + (3*nto:-1:-3*nto);
xTSPSelIv = xTSPSel(selectInvIndex, :);
compSignal = zeros(length(yr), inChannel, 4);
compSignalRef = zeros(length(yr), 4);
baseIdx = (1:length(yr));
testIdx = (1:length(xTest)) + fs;

compSignal(:, :, 1) = fftfilt(xTSPSelIv(:, 1), yr);
compSignal(:, :, 2) = fftfilt(xTSPSelIv(:, 2), yr);
compSignal(:, :, 3) = fftfilt(xTSPSelIv(:, 3), yr);
compSignal(:, :, 4) = fftfilt(xTSPSelIv(:, 4), yr);

compSignalRef(testIdx, 1) = fftfilt(xTSPSelIv(:, 1), xTest);
compSignalRef(testIdx, 2) = fftfilt(xTSPSelIv(:, 2), xTest);
compSignalRef(testIdx, 3) = fftfilt(xTSPSelIv(:, 3), xTest);
compSignalRef(testIdx, 4) = fftfilt(xTSPSelIv(:, 4), xTest);

orthogonalSignal = zeros(length(yr), inChannel, 4);
orthogonalSignalRef = zeros(length(yr), 4);
for ii = 1:8
    tmpIdx = min(length(yr), baseIdx + (ii - 1) * nto);
    for jj = 1:4
        orthogonalSignal(:, :, jj) = orthogonalSignal(:, :, jj) ...
            + B4(jj, ii) * compSignal(tmpIdx, :, jj);
        orthogonalSignalRef(:, jj) = orthogonalSignalRef(:, jj) ...
            + B4(jj, ii) * compSignalRef(tmpIdx, jj);
    end
end


output.fs = fs;
output.tspName = tspName;
output.tResponse = tResponse;
output.nRepetition = nRepetition;
output.calibrationConst = calibrationConst;
output.fLow = fLow;
output.pinkLPC = pinkLPC;
output.xTestPink = xTestPinkOrg;
output.outChannel = outChannel;
output.numChannels = numChannels;
output.selectedChannels = selectedChannels;
output.testSignalDuration = tTestSignal;
output.yRecorded = y;
output.yRecovered = yr;
output.orthogonalSignal = orthogonalSignal;
output.orthogonalSignalRef = orthogonalSignalRef;
output.elapsedTime = toc(startTic);
output.lAeq = lAeq;
output.elapsedTimeRecover = toc(middleTic);
end

