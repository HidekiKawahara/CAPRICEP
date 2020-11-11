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
figureHandleReport = figure;
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
semilogx(fx, 10*log10(smLongPowerSpec)-totalLevel, 'linewidth', 2);grid on;
hold all
semilogx(fx, 10*log10(smShortPowerSpec)-totalLevel, 'linewidth', 2);
semilogx(fx, 10*log10(smDevPowerSpec)-totalLevel, 'linewidth', 2);
semilogx(fx, 10*log10(smRandPowerSpec)-totalLevel, 'linewidth', 2);
semilogx(fx, 10*log10(smPrePowerSpec)-totalLevel, 'linewidth', 2);
set(gca, 'fontsize', 16, 'linewidth', 2)
legend('LTI-L', 'LTI-S', 'nonL-TI', 'RNTV', 'preBG', 'Location', "south", "numcolumns", 5);
axis([10 fs/2 -60 20])
xlabel('frequency (Hz)')
ylabel('gain (rel. dB)')
title("Tr:" + num2str(specStr.tResponse) + " ms " + specStr.outChannel + " LAeq:" + num2str(specStr.lAeq) + " dB" );

output.specStr = specStr;
output.figureHandleReport = figureHandleReport;

end