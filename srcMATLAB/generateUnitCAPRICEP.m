function output = generateUnitCAPRICEP(fs, averageFreqDistance, mag, ab, npFFT)
% generate unit extended TSP called CAPRICEP
% output = generateUnitCAPRICEP(fs, averageFreqDistance, mag, ab, npFFT)
%
% Argument
%   fs                   : sampling frequency (Hz)
%   averageFreqDistance  : average distance of phase manipulation (Hz)
%   mag                  : stretching manipulation width
%   ab                   : distribution shape controller
%   npFFT                : 2's exponent for FFT buffer size
%
% Output
%   output      : structure with the following fields
%    elapsedTime: elapsed time for procesing (s)
%             tx: time axis (s)
%           xTSP: generated exTSP signal
%           fftl: FFT siz
%             fx: freqiemcu axis (Hz)
%              h: frequency domain representation
%             gd: group delay (s)

%   Copyright 2020 Hideki Kawahara
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

startTic = tic;
fftl = 2 ^ npFFT;
fx = (0:fftl-1)'/fftl * fs;
nDivision = floor(fs / 2 / averageFreqDistance);
z = exp(-2 * 1i * pi * fx / fs);
idx = (1:fftl)';
idxF = circshift(idx, -1);
h = ones(fftl, 1);
randSeq1 = betaincinv(rand(nDivision + 1, 1), ab, ab);
rlocation = cumsum(randSeq1) / sum(randSeq1) * fs / 2;
rlocation = rlocation(1:end-1);
randSeq2 = rand(nDivision, 1);
for ii = 1:nDivision
    rpolarity = sign(randSeq2(ii) - 0.5);
    a = exp(-pi * (mag * averageFreqDistance) / fs) ...
        * exp(2 * 1i * pi * rlocation(ii) / fs);
    if rpolarity > 0
        h = h .* (z - conj(a)) .* (z - a) ./ (1 - a .* z) ./ (1 - conj(a) .* z);
    else
        h = h .* (1 - a .* z) .* (1 - conj(a) .* z) ./ (z - conj(a)) ./ (z - a);
    end
end
tx = ((1:fftl)' - fftl / 2) / fs;
gd1 = -angle(h(idxF) ./ h) / (fs / fftl * 2 * pi);
output.elapsedTime = toc(startTic);
output.tx = tx;
output.xTSP = fftshift(real(ifft(h)));
output.fftl = fftl;
output.fx = fx;
output.h = h;
output.gd = gd1;
end

