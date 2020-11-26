%%

close all
clear variables
%%
load setCAPRICEP025ms
%%
figure('Position',[61         559         500         440]);
semilogy(sort(corrList), (1:length(corrList))/length(corrList), 'k', 'linewidth', 3);
hold all
semilogy(sort(corrList), (length(corrList):-1:1)/length(corrList), 'k', 'linewidth', 1);
grid on;
set(gca, 'FontSize', 16, 'LineWidth', 2)
axis([0.06 0.17 10^(-5) 1])
ylabel('cumulative distribution')
xlabel('correlation')
print -depsc correlationDistributionN.eps
%%
load setCAPRICEP200ms
B4 = [1 1 1 1 1 1 1 1;1 -1 1 -1 1 -1 1 -1;1 1 -1 -1 1 1 -1 -1;1 1 1 1 -1 -1 -1 -1];
%% test signal
targetERD = 0.2;
samplr = abs(output.timeAxis) < 2 * targetERD;
lSample = length(output.timeAxis(samplr));
nSkip = round(targetERD * fs);
nRepetition = 44;
xTest = zeros(lSample + nSkip * nRepetition + 1, 1);
aSign = 1;
aaSign = [1 1 -1 -1];
for ii = 1:nRepetition
    xTest((ii - 1) * nSkip + (1:lSample)) = xTest((ii - 1) * nSkip + (1:lSample)) ...
        + xTSPs(samplr, 1) + aSign * xTSPs(samplr, 2) ...
        + aaSign(rem(ii - 1, 4) + 1) * xTSPs(samplr, 3);
    aSign = -aSign;
end
%%

orthSig = zeros(length(xTest) + lSample + nSkip * 8 + 1, 4);
B4 = [1 1 1 1 1 1 1 1;1 -1 1 -1 1 -1 1 -1;1 1 -1 -1 1 1 -1 -1;1 1 1 1 -1 -1 -1 -1];
for ii = 1:4
xx = xTSPs(end:-1:1, ii);
cmpTmp = fftfilt(xx(samplr), [xTest;zeros(lSample,1)]);
    for jj = 1:8
        tmpCf = B4(:, circshift([1 2 3 4 5 6 7 8], -(jj - 1)));
        orthSig((jj -1) * nSkip + (1:length(cmpTmp)), ii) = ...
            orthSig((jj -1) * nSkip + (1:length(cmpTmp)), ii) + tmpCf(ii, 1) * cmpTmp;
    end
end
txo = (1:size(orthSig, 1))' / fs;
figure;plot(txo, 20*log10(abs(orthSig(:,1))/8));grid on;
axis([txo(1) txo(end) -150 5]) 
%%
figure('Position',[681   466   560   513]);
subplot(511)
xx = xTSPs(end:-1:1, 1);
tx = (1:length(fftfilt(xx(samplr), [xTest;zeros(lSample,1)])))' / fs;
plot(tx,fftfilt(xx(samplr), [xTest;zeros(lSample,1)]), 'k');grid on;
set(gca, 'FontSize', 16, 'LineWidth', 2)
ylabel('q1')
axis([tx(1) tx(end) -1.1 1.1])
%figure('Position',[680   558   560   240]);
subplot(512)
xx = xTSPs(end:-1:1, 2);plot(tx, fftfilt(xx(samplr), [xTest;zeros(lSample,1)]), 'k');grid on;
set(gca, 'FontSize', 16, 'LineWidth', 2)
ylabel('q2')
axis([tx(1) tx(end) -1.1 1.1])
%figure('Position',[680   558   560   240]);
subplot(513)
xx = xTSPs(end:-1:1, 3);plot(tx, fftfilt(xx(samplr), [xTest;zeros(lSample,1)]), 'k');grid on;
set(gca, 'FontSize', 16, 'LineWidth', 2)
ylabel('q3')
axis([tx(1) tx(end) -1.1 1.1])
subplot(514)
ritr1 = fftfilt(xx(samplr), [xTest;zeros(lSample,1)]);
set(gca, 'FontSize', 16, 'LineWidth', 2)
plot(tx, 20*log10(abs(ritr1)), 'k');grid on;
ylabel('q1 (dB)')
set(gca, 'FontSize', 16, 'LineWidth', 2)
axis([tx(1) tx(end) -130 10])
subplot(515)
txo = (1:size(orthSig, 1))' / fs;
plot(txo - targetERD * 3, 20*log10(abs(orthSig(:,1))/8), 'k');grid on;
set(gca, 'FontSize', 16, 'LineWidth', 2)
axis([tx(1) tx(end) -130 10])
xlabel('time (s)')
ylabel('r1R (dB)')
print -deps orthSeqEx.eps

%%
fs = 44100;
npFFT = 15;
avFd = 40;
magList = 2 .^ (1/4);
abList = 7;
iteration = 1280 * 8;
npFFTSh = 9;
avFdSh = 1000;
magSh = 2;
abSh = 2;
selectrSig = (1:2^npFFT)' + 2^(npFFTSh - 1);
shapeFig = figure;
statFig = figure;
displayInterval = 100;
for ii = 1:length(magList)
    mag = magList(ii);
    for kk = 1:length(abList)
        rmsEnvTmp = zeros(2^npFFT, 1);
        wsDList = zeros(round(iteration / displayInterval) +1, 1) + NaN;
        erdList = zeros(round(iteration / displayInterval) +1, 1) + NaN;
        dispCount = 0;
        for jj = 1:iteration
            outp = generateUnitCAPRICEP(fs, avFd, mag, abList(kk), npFFT);
            outpSh = generateUnitCAPRICEP(fs, avFdSh, magSh, abSh, npFFTSh);
            xTSPfiltered = fftfilt(outpSh.xTSP, [outp.xTSP; zeros(2^npFFTSh, 1)]);
            xTSPfiltered = xTSPfiltered(selectrSig);
            rmsEnvTmp = rmsEnvTmp + xTSPfiltered .^ 2;
            if rem(jj, displayInterval) == 0
                dispCount = dispCount + 1;
                durRef = 0.564/avFd/sqrt(mag);
                totalEngy = sum(rmsEnvTmp);
                durList = (round(durRef * fs):round(3 * durRef * fs))';
                ctrIdx = outp.fftl / 2;
                tx = outp.tx;
                wsMeasure = zeros(length(durList), 1);
                for jk = 1:length(durList)
                    idx = ctrIdx + (-durList(jk):durList(jk));
                    wsMeasure(jk) = sum(abs(totalEngy / (length(idx)) - rmsEnvTmp(idx))) ...
                        + sum(rmsEnvTmp(tx < tx(idx(1)) | tx > tx(idx(end))));
                end
                [minV, minIdx]  = min(wsMeasure);
                wsDList(dispCount) = minV / jj;
                erdList(dispCount) =  (2*durList(minIdx)+1) / fs;
                tmpShape = zeros(outp.fftl, 1);
                tmpShape(outp.fftl / 2 + (-durList(minIdx):durList(minIdx))) = ...
                    totalEngy / (2*durList(minIdx)+1);
                figure(shapeFig)
                plot(outp.tx, tmpShape/jj, 'g', 'LineWidth', 3)
                hold all
                plot(outp.tx, rmsEnvTmp/jj, 'k');grid on;
                set(gca, 'xlim', durRef * [-4 4])
                title([num2str(jj)])
                hold off
                drawnow
                figure(statFig);
                subplot(211)
                plot((1:length(wsDList)) * displayInterval, wsDList);grid on;
                subplot(212)
                plot((1:length(wsDList)) * displayInterval, erdList);grid on;
                drawnow
            end
        end
    end
end
%%
figure('Position',[61         559         834         300]);
plot(outp.tx * avFd, tmpShape/jj, 'LineWidth', 5, 'color', 0.7 * [1 1 1]);grid on;
hold all
plot(outp.tx * avFd, rmsEnvTmp/jj, 'k', 'LineWidth', 2);grid on;
set(gca, 'FontSize', 20, 'LineWidth', 2, 'xlim', durList(minIdx) * [-2 2] * avFd / fs)
xlabel('time (normalized by 1/F_d)')
ylabel('sample variance')
print -depsc fixedShapingLinN.eps
%%
plot(tx(abs(tx) <  2 * durList(minIdx) / fs) * avFd, ...
    rmsEnvTmp(abs(tx) <  2 * durList(minIdx) / fs) / iteration);
grid on;
hold all
rectModel = tx * 0;
rectModel(abs(tx) <  durList(minIdx) / fs) = ...
    totalEngy / length(tx(abs(tx) <  durList(minIdx) / fs) ) / iteration;
plot(tx(abs(tx) <  2 * durList(minIdx) / fs) * avFd,...
    rectModel(abs(tx) <  2 * durList(minIdx) / fs), 'linewidth', 5, 'color', 0.7 * [1 1 1])
plot(tx(abs(tx) <  2 * durList(minIdx) / fs) * avFd, ...
    rmsEnvTmp(abs(tx) <  2 * durList(minIdx) / fs) / iteration, 'k');
set(gca, 'FontSize', 20, 'LineWidth', 2, 'xlim', 2 * durList(minIdx) / fs * [-1 1] * avFd)
xlabel('time (normalized by 1/F_d)')
ylabel('sample variance')
print -depsc rawShapingLinN.eps
%%
close all
clear variables
%%
fs = 44100;
npFFT = 15;
avFd = 40;
magList = 2 .^ (1/4);
abList = 7;
iteration = 1280 * 8;
npFFTSh = 9;
avFdSh = 1000;
magSh = 2;
abSh = 2;
outpSh = generateUnitCAPRICEP(fs, avFdSh, magSh, abSh, npFFTSh);
figure;plot(outpSh.xTSP);grid on;
%soundsc(outpSh.xTSP, fs)
%%
nXTSP = 10000;
xTSPset = zeros(length(outpSh.xTSP), nXTSP);
for ii = 1:nXTSP
    outpSh = generateUnitCAPRICEP(fs, avFdSh, magSh, abSh, npFFTSh);
    xTSPset(:, ii) = outpSh.xTSP;
end
%%
[x, fs] = audioread('PROLOGUEc19t27.wav');
snrL = zeros(nXTSP, 1);
snrR = zeros(nXTSP, 1);
dB1 = 20 * log10(std(x));
for ii = 1:nXTSP
    y1 = fftfilt(xTSPset(:, ii), [x;zeros(512,2)]);
    y1 = y1(256 + (1:length(x))');
    tmpDiff = y1 - x;
    dBdiff = 20 * log10(std(tmpDiff(512:end-512)));
    tmpSNR = dBdiff - dB1;
    snrL(ii) = tmpSNR(1);
    snrR(ii) = tmpSNR(2);
end
%%
figure('Position',[61         559        [ 400         300] * 1.3]);
plot(sort(-snrL(:)), (1:nXTSP)/nXTSP, 'k', 'LineWidth', 2);grid on;
hold all;plot(sort(-snrR(:)), (1:nXTSP)/nXTSP, 'k--', 'LineWidth', 2);grid on;
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('SNR (dB)')
ylabel('cumulative distribution');
legend('Lch', 'Rch','location', 'southeast')
print -depsc snrModn.eps

%%
[x, fs] = audioread('PROLOGUEc19t27.wav');
snrL = zeros(nXTSP, nXTSP);
snfR = zeros(nXTSP, nXTSP);
for ii = 1:nXTSP
    y1 = fftfilt(xTSPset(:, ii), x);
    dB1 = 20 * log10(std(y1));
    for jj = 1:nXTSP
        y2 = fftfilt(xTSPset(:, jj), x);
        dB2 = 20 * log10(std(y2));
        dBdiff = 20 * log10(std(y1 - y2));
        tmpSNR = dBdiff - (dB1 + dB2) /2;
        snrL(ii,jj) = tmpSNR(1);
        snrR(ii,jj) = tmpSNR(2);
    end
end
%%
figure('Position',[61         559        [ 400         300] * 1.3]);
plot(sort(-snrL(:)), (1:10000)/10000, 'k', 'LineWidth', 2);grid on;
hold all;plot(sort(-snrR(:)), (1:10000)/10000, 'k--', 'LineWidth', 2);grid on;
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('SNR (dB)')
ylabel('cumulative distribution');
legend('Lch', 'Rch','location', 'southeast')
print -depsc snrModn.eps

%%
fs = 44100;
npFFT = 15;
avFd = 40;
magList = 2 .^ (1/4);
abList = 7;
iteration = 1280 * 8;
npFFTSh = 9;
avFdSh = 1000;
magSh = 2;
abSh = 2;
outpSh = generateUnitCAPRICEP(fs, avFdSh, magSh, abSh, npFFTSh);
figure;plot(outpSh.xTSP);grid on;
%%
nXTSP = 100;
xTSPset = zeros(length(outpSh.xTSP), nXTSP);
for ii = 1:nXTSP
    outpSh = generateUnitCAPRICEP(fs, avFdSh, magSh, abSh, npFFTSh);
    xTSPset(:, ii) = outpSh.xTSP;
end
%%
[x, fs] = audioread('aaaaa.wav');
%%

xtrm = x(fs:fs*2);
figure;
semilogy(sort(xtrm), (1:length(xtrm))/length(xtrm),'k:', 'linewidth', 5);grid on;
y = fftfilt(xTSPset(:,1), x);
ytrm = y(fs:fs*2);
hold all
semilogy(sort(ytrm), (1:length(xtrm))/length(xtrm),'k');grid on;
semilogy(sort(xtrm), (length(xtrm):-1:1)/length(xtrm),'k:', 'linewidth', 5);grid on;
semilogy(sort(ytrm), (length(xtrm):-1:1)/length(xtrm),'k');grid on;
for ii = 6:17:100
    y = fftfilt(xTSPset(:,ii), x);
    ytrm = y(fs:fs*2);
    semilogy(sort(ytrm), (1:length(xtrm))/length(xtrm),'k');
    semilogy(sort(ytrm), (length(xtrm):-1:1)/length(xtrm),'k');
end
axis([-0.25 0.25 10^(-4) 1])
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('sampled value');
ylabel('cumulative distribution')
print -deps vowelAprob.eps

%%
[xx, fs] = audioread('nkym18t95eX.wav');
x = xx(:,1);
%
nXTSP = 1000;
xTSPset = zeros(length(outpSh.xTSP), nXTSP);
for ii = 1:nXTSP
    outpSh = generateUnitCAPRICEP(fs, avFdSh, magSh, abSh, npFFTSh);
    xTSPset(:, ii) = outpSh.xTSP;
end
%%

xtrm = x(10001:25000);
figure;
semilogy(sort(xtrm), (1:length(xtrm))/length(xtrm),'k', 'linewidth', 1);grid on;
y = fftfilt(xTSPset(:,1), x);
ytrm = y(10001:25000);
hold all
semilogy(sort(ytrm), (1:length(xtrm))/length(xtrm),'k');grid on;
semilogy(sort(xtrm), (length(xtrm):-1:1)/length(xtrm),'k', 'linewidth', 1);grid on;
semilogy(sort(ytrm), (length(xtrm):-1:1)/length(xtrm),'k');grid on;
for ii = 1:nXTSP
    y = fftfilt(xTSPset(:,ii), x);
    ytrm = y(10001:25000);
    semilogy(sort(ytrm), (1:length(xtrm))/length(xtrm),'color', 0.6*[1 1 1]);
    semilogy(sort(ytrm), (length(xtrm):-1:1)/length(xtrm),'color', 0.6*[1 1 1]);
end
semilogy(sort(xtrm), (1:length(xtrm))/length(xtrm),'k', 'linewidth', 3);grid on;
semilogy(sort(xtrm), (length(xtrm):-1:1)/length(xtrm),'k', 'linewidth', 3);grid on;
axis([-0.25 0.25 10^(-4) 1])
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('sampled value');
ylabel('cumulative distribution')
print -depsc vowelEprob.eps

