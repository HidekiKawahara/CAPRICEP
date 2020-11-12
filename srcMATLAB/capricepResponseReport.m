function output = capricepResponseReport(specStr)
% Report generation from the analysis results
% output = capricepResponseReport(specStr)
% 
%  This version is tentative
%
% Argument
%  specStr    : structure variable with the following fields
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
%
% Output
%    output    : structure variable with the following fields
%           specStr                : copy of the analysis results
%           figureHandleRaw        : handle of the raw visualization
%           figureHandleReport     : handle of the report
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



fx = specStr.frequencyAxis;
fs = specStr.fs;
longSpec = specStr.longSpec;
longPowerSpec = abs(longSpec) .^ 2;
shortSpec = specStr.shortSpec;
shortPowerSpec = abs(shortSpec) .^ 2;
prePowerSpec = specStr.prePowerSpec;
totalLevel = 10*log10(mean(longPowerSpec));
fftl = length(fx);

if 1 == 2
    figureHandleRaw = figure;
    semilogx(fx, 10*log10(longPowerSpec)-totalLevel, 'linewidth', 2);grid on;
    hold all
    semilogx(fx, 10*log10(shortPowerSpec)-totalLevel, 'linewidth', 2);
    semilogx(fx, 10*log10(specStr.deviationPowerSpec)-totalLevel, 'linewidth', 2);
    semilogx(fx, 10*log10(specStr.randomPowerSpec)-totalLevel, 'linewidth', 2);
    semilogx(fx, 10*log10(prePowerSpec)-totalLevel, 'linewidth', 2);
    set(gca, 'fontsize', 16, 'linewidth', 2)
    legend('LTI-L', 'LTI-S', 'nonL-TI', 'RNTV', 'preBG', 'Location', "south", "numcolumns", 5);
    axis([10 fs/2 -80 15])
    xlabel('frequency (Hz)')
    ylabel('gain (rel. dB)')
    output.figureHandleRaw = figureHandleRaw;
end

%--
fxH = fx * 2^(1/6);
fxL = fx * 2^(-1/6);
fbw = fxH - fxL;
fbw(1) = fbw(2);
figureHandleReport = figure('unit', 'centimeters');
%set(figureHandleReport, 'position',[23.9889   19.6850   19.7556   14.8167]);
set(figureHandleReport, 'position', [33    7.0556   21   29.7])
subplot(311);
smLongPowerSpec = (interp1(fx, cumsum(longPowerSpec) * fx(2), fxH, 'linear', 'extrap') ...
    - interp1(fx, cumsum(longPowerSpec) * fx(2), fxL, 'linear', 'extrap')) ./ fbw;
smShortPowerSpec = (interp1(fx, cumsum(shortPowerSpec) * fx(2), fxH, 'linear', 'extrap') ...
    - interp1(fx, cumsum(shortPowerSpec) * fx(2), fxL, 'linear', 'extrap')) ./ fbw;
smDevPowerSpec = (interp1(fx, cumsum(specStr.deviationPowerSpec) * fx(2), fxH, 'linear', 'extrap') ...
    - interp1(fx, cumsum(specStr.deviationPowerSpec) * fx(2), fxL, 'linear', 'extrap')) ./ fbw;
smRandPowerSpec = (interp1(fx, cumsum(specStr.randomPowerSpec) * fx(2), fxH, 'linear', 'extrap') ...
    - interp1(fx, cumsum(specStr.randomPowerSpec) * fx(2), fxL, 'linear', 'extrap')) ./ fbw;
smPrePowerSpec = (interp1(fx, cumsum(prePowerSpec) * fx(2), fxH, 'linear', 'extrap') ...
    - interp1(fx, cumsum(prePowerSpec) * fx(2), fxL, 'linear', 'extrap')) ./ fbw;
semilogx(fx, 20*log10(abs(specStr.shortSpec))-totalLevel, 'linewidth', 1, 'color', 0.6*[1 1 1]);grid on;
hold all
%semilogx(fx, 10*log10(smLongPowerSpec)-totalLevel, 'linewidth', 2);grid on;
semilogx(fx, 10*log10(smShortPowerSpec)-totalLevel, 'k','linewidth', 2);
semilogx(fx, 10*log10(smDevPowerSpec)-totalLevel, 'linewidth', 2);
semilogx(fx, 10*log10(smRandPowerSpec)-totalLevel, 'linewidth', 2);
semilogx(fx, 10*log10(smPrePowerSpec)-totalLevel, 'linewidth', 2);
set(gca, 'fontsize', 14, 'linewidth', 2, 'fontname', 'Helvetica')
legend('LTI-L', 'LTI-S', 'nonL-TI', 'RNTV', 'preBG', 'Location', "south", "numcolumns", 5);
axis([10 fs/2 -70 20])
xlabel('frequency (Hz)')
ylabel('gain (rel. dB)')
title("Tr:" + num2str(specStr.tResponse) + " ms " + specStr.outChannel + " LAeq:" + num2str(specStr.lAeq) + " dB" );
subplot(529);
plot(specStr.timeAxisShort * 1000, ...
    specStr.averageShortResponse/max(abs(specStr.averageShortResponse)), 'k');grid on;
axis([-0.3 2 -1 1])
set(gca, 'fontsize', 14, 'linewidth', 2, 'fontname', 'Helvetica')
xlabel('time (ms)')
subplot(5,2,10);
pwResp = 10*log10(abs(fftfilt(hanning(21), specStr.averageShortResponse.^2)));
txTmp = specStr.timeAxisShort * 1000;
plot(txTmp, pwResp - max(pwResp), 'k');grid on
axis([txTmp(1) txTmp(end) -90 0]);
set(gca, 'fontsize', 14, 'linewidth', 2, 'fontname', 'Helvetica')
xlabel('time (ms)')
ylabel('level (dB)')

subplot(527)
resp = specStr.averageShortResponse;
fcv = fft(resp, fftl);
rgd = -angle(fcv([2:fftl fftl]) ./ fcv) /2/pi/fx(2);
semilogx(fx, (rgd - 0.01) * 1000 ,'k');grid on;
set(gca, 'fontsize', 14, 'linewidth', 2, 'fontname', 'Helvetica')
axis([50 1000 -15 30])
xlabel('frequency (Hz)')
ylabel('group delay (ms)')
tmpPos = get(gca, 'position');
set(gca, 'position' , tmpPos + [0 0.03 0 0])

subplot(513)
fcvt = fft(resp(882:end), fftl);
rgdt = -angle(fcvt([2:fftl fftl]) ./ fcvt) /2/pi/fx(2);
cumRgdt = cumsum(abs(fcvt([2:fftl fftl])) .* abs(fcvt) .* rgdt);
cumFcvt = cumsum(abs(fcvt([2:fftl fftl])) .* abs(fcvt));
avGdtH = interp1(fx, cumRgdt, fxH, 'linear', 'extrap');
avGdtL = interp1(fx, cumRgdt, fxL, 'linear', 'extrap');
avPwtH = interp1(fx, cumFcvt, fxH, 'linear', 'extrap');
avPwtL = interp1(fx, cumFcvt, fxL, 'linear', 'extrap');
avPw = avPwtH - avPwtL;
avPw(1) = avPw(2);
avGdt = (avGdtH - avGdtL) ./ avPw;
dispMax = (max(avGdt(fx > 100 & fx < 10000)) + 0.01) * 6.9077;

semilogx(fx, (avGdt+0.01) * 6.9077,'linewidth', 2);grid on;
axis([10 fs/2 0 dispMax * 1.1]);
set(gca, 'fontsize', 14, 'linewidth', 2, 'fontname', 'Helvetica')
xlabel('frequency (Hz)')
ylabel('rev. time (s)')

fxH = fx * 2^(1/2);
fxL = fx * 2^(-1/2);
cumRgdt = cumsum(abs(fcvt([2:fftl fftl])) .* abs(fcvt) .* rgdt);
cumFcvt = cumsum(abs(fcvt([2:fftl fftl])) .* abs(fcvt));
avGdtH = interp1(fx, cumRgdt, fxH, 'linear', 'extrap');
avGdtL = interp1(fx, cumRgdt, fxL, 'linear', 'extrap');
avPwtH = interp1(fx, cumFcvt, fxH, 'linear', 'extrap');
avPwtL = interp1(fx, cumFcvt, fxL, 'linear', 'extrap');
avPw = avPwtH - avPwtL;
avPw(1) = avPw(2);
avGdt = (avGdtH - avGdtL) ./ avPw;
hold all
semilogx(fx, (avGdt+0.01) * 6.9077,'linewidth', 2);grid on;
tmpPos = get(gca, 'position');
set(gca, 'position' , tmpPos + [0 0.05 0 0])

subplot(528)
plot([0 1],[0 1], '.');
axis off
grid off
set(gca, 'fontsize', 14, 'linewidth', 2, 'fontname', 'Helvetica')
text(0, 0.9, specStr.tspName);
if isfield(specStr, 'dataFile')
    text(0, 0.7, specStr.dataFile);
end

output.specStr = specStr;
output.figureHandleReport = figureHandleReport;

end