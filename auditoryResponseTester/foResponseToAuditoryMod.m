function output = foResponseToAuditoryMod(fullPathName, displayInd)
%% ausitory stimulation to speech fo response analysis
%   Argument
output = struct;
startTic = tic;
if displayInd == 2
    displayInd = 1;
    debug = 1;
else
    debug = 0;
end

[x,fs] = audioread(fullPathName);
FMcent = 25;
xOrg = x;
fileInfo = audioinfo(fullPathName);
tmp = sscanf(fileInfo.Comment, '%10f');
if length(tmp) > 1
    fo = tmp(1);
    foVoice = fo;
    calibrationConst = tmp(2);
    if length(tmp) > 2
        foVoice = tmp(3);
        if length(tmp) > 3
            FMcent = tmp(4);
        end
    end
else
    fo = tmp(1);
    calibrationConst = 0;
end

if displayInd
    disp(fileInfo)
    figure;
    subplot(211)
    plot((1:length(x))/fs, x(:,1));grid on;
    xlabel('time (s)')
    title(fullPathName);% ; ['L-ch signal fileInfo:' fileInfo.Comment]});
    subplot(212)
    plot((1:length(x))/fs, x(:,2));grid on;
    xlabel('time (s)')
    title('R-ch signal: ')
end
%x(:,1) = x(:,1) + 0.1 * x(:,1) .^2;

%% Test signal data read
testData = load('testSigBase400ms.mat');
xTSPSel = testData.xTSPSel;
c = testData.cOrdered;
fftl = testData.fftl;
fftlSeg = testData.fftlSeg;
output.fftl = fftl;
output.fftlSeg = fftlSeg;
w = sixtermCos(round(12*fs/fo));
tww = ((0:length(w)-1) - length(w) / 2)' / fs;
hw = w / sum(w) .* exp(1i * 2 * pi * fo * tww);
wV = sixtermCos(round(12*fs/foVoice));
twwV = ((0:length(wV)-1) - length(wV) / 2)' / fs;
hwV = wV / sum(wV) .* exp(1i * 2 * pi * foVoice * twwV);
if debug
    hw2 = w / sum(w) .* exp(1i * 2 * pi * fo * tww * 2);
end
%% Low cut filter design
fftlShort = 8192;
foLower = min(fo, foVoice) * 0.7;
filtGain = ones(fftlShort, 1);
fx = (0:fftlShort -1 ) /fftlShort * fs;
fxBi = fx;
fxBi(fx > fs / 2) = fxBi(fx > fs / 2) - fs;
filtGain(foLower > abs(fxBi)) = 0;
filtResp = fftshift(ifft(filtGain)) .* blackman(fftlShort);
wTmp = fftfilt(hw, [filtResp;zeros(fftlShort, 1)]);
hw = wTmp( round(6*fs/fo) + (1:fftlShort));
wTmp = fftfilt(hwV, [filtResp;zeros(fftlShort, 1)]);
hwV = wTmp( round(6*fs/foVoice) + (1:fftlShort));
%% Check test signal type

switch fileInfo.Artist
    case {'MFND','MFNDH'}
        if ~isnan(str2double(fileInfo.Title))
            combId = str2double(fileInfo.Title);
        else
            combId = 1;
        end
        outputg = generateTestFMSignal(testData, fo, fileInfo.Artist, combId, FMcent, x(:,1));
        segloc = (1:testData.fftlSeg * 8);
        checkLoc = ((1:testData.fftlSeg * 16));
        %slidingCorr = fftfilt(flip(x(segloc + testData.fftl*2, 2)), ...
        %    output.testSignal(segloc + testData.fftl*2));
        slidingCorr = fftfilt(flip(xOrg(segloc + testData.fftl*2, 2)), ...
            outputg.testSignal(checkLoc + testData.fftl*2));
        %xslidingCorr = fftfilt(flip(xOrg(segloc + testData.fftl*2, 2)), ...
        %    outputg.testSignal(checkLoc + testData.fftl*2));
        [maxCorr, idxx] = max(slidingCorr);
        if maxCorr < 0.9
            disp('test signal does not match with this test');
        else
            outputg = generateTestFMSignal(testData, fo, 'SINE', combId, FMcent, x(:,1));
            x(:,2) = outputg.testSignal(min(length(x), max(1, idxx + (1:length(x)) - testData.fftlSeg * 8)));
        end
    otherwise
        if ~isnan(str2double(fileInfo.Title))
            combId = str2double(fileInfo.Title);
        else
            combId = 1;
        end
end
output.fs = fs;
output.fo = fo;
output.signal = x;
output.fileInfo = fileInfo;
output.calibrationConst = calibrationConst;
output.combination = c(combId, :);

%% Fo extraction
xanl = fftfilt([hwV hw], [x; zeros(fftlShort, 2)]);
xanl = xanl(fftlShort/2 + (1:length(x)), :);
fxri = angle(xanl([2:end end],:) ./ xanl) / 2 / pi * fs;
if debug
    xanl2 = fftfilt(hw2, x);
    fxri2 = angle(xanl2([2:end end],:) ./ xanl2) / 2 / pi * fs / 2;
end

%% Display acquired fo trajectory
if displayInd
    figure;
    subplot(211)
    if debug
        plot((1:length(x))/fs, fxri2(:,1));grid on;
    else
        plot((1:length(x))/fs, fxri(:,1));grid on;
    end
    title('L-ch ');
    xlabel('time (s)')
    ylabel('fund. freq. (Hz)')
    axis([0 length(x)/fs foVoice/2^(1/3) foVoice*2^(1/3)])
    subplot(212)
    plot((1:length(x))/fs, fxri(:,2));grid on;
    title('R-ch ');
    xlabel('time (s)')
    ylabel('fund. freq. (Hz)')
    axis([0 length(x)/fs fo/2^(1/3) fo*2^(1/3)])
end
tmp1 =  max(foVoice/2^0.25,min(foVoice*2^0.25, fxri(:, 1)));
tmp2 = max(fo/2^0.25,min(fo*2^0.25, fxri(:, 2)));
rawFoi = [tmp1, tmp2];
%rawFoi = max(fo/2^0.25,min(fo*2^0.25, fxri));
if debug
    rawFoi2 = max(fo/2^0.25,min(fo*2^0.25, fxri2));
    rawFoi(:,1) = rawFoi2(:,1);
end

%% Pulse recovery
index1 = 1:10;
%minPeakIdx = 1;
otherElement = index1(~ismember(index1, c(combId,:)));
%element4 = otherElement(randi(7));
setId = [c(combId,:), otherElement];
output.rawFoi = rawFoi; % for debug
output.filterOutput = xanl;
output.rawFoValues = fxri;

avFoOut = averageReliableFo(output);
averageFoCent = 1200*log2(avFoOut.averageFo ./ [foVoice fo]);
averageFoCent(isnan(averageFoCent)) = 0;

rawFoCent = [1200 * log2(tmp1/foVoice), 1200 * log2(tmp2/fo)];
rawFoCent = rawFoCent - averageFoCent;
recoveredFoCentSp = fftfilt(xTSPSel(end:-1:1, setId), rawFoCent(:,1));
recoveredFoCentTx = fftfilt(xTSPSel(end:-1:1, setId), rawFoCent(:,2));

baseIdxSig = (1:length(x));
fixedIdx = circshift(baseIdxSig, -fftl/2);
if displayInd
    figure;
    subplot(211)
    plot((1:length(x))/fs, recoveredFoCentSp(fixedIdx, 1:3));grid on;
    %axis([0 length(x)/fs fo/2^(1/3) fo*2^(1/3)])
    subplot(212)
    plot((1:length(x))/fs, recoveredFoCentTx(fixedIdx, 1:3));grid on;
    %axis([0 length(x)/fs fo/2^(1/3) fo*2^(1/3)])
end

%% Orthogonalization
bMat = [1 1 1 1;1 -1 1 -1;1 1 -1 -1]';
bMatEx = [bMat ones(4, 7); bMat -ones(4, 7)];
orthFoCentSp = recoveredFoCentSp(fixedIdx, :);
orthFoCentTx = recoveredFoCentTx(fixedIdx, :);
for ii = 2:8
    rotatedIdx = circshift((1:length(x)), -(ii-1)*fftlSeg);
    orthFoCentSp = orthFoCentSp + recoveredFoCentSp(rotatedIdx, :) ...
        * diag([bMatEx(ii, 1), bMatEx(ii, 2), bMatEx(ii, 3), bMatEx(ii, 4:10)]);
    orthFoCentTx = orthFoCentTx + recoveredFoCentTx(rotatedIdx, :) ...
        * diag([bMatEx(ii, 1), bMatEx(ii, 2), bMatEx(ii, 3), bMatEx(ii, 4:10)]);
end
%%
orthIndex = circshift((1:length(x)), -4*fftlSeg);
orthFoCentTx = orthFoCentTx(orthIndex ,:);
orthFoCentSp = orthFoCentSp(orthIndex ,:);

%% Display raw result
txSig = (1:length(x))/fs;
if displayInd
    figure;plot(txSig, (orthFoCentTx(:,1)+ orthFoCentTx(:,2)+orthFoCentTx(:,3)*2)/4/8);grid on;
    hold all;plot(txSig, (orthFoCentSp(:,1)+ orthFoCentSp(:,2)+orthFoCentSp(:,3)*2)/4/8);grid on;
    hold all;plot(txSig, sum(orthFoCentSp(:,4:10),2)/8/sqrt(7*2.5))
    title([fileInfo.Artist ' ' fileInfo.Comment ])
end

%% Display synchronized averaged result
recovered4Signal = (orthFoCentSp(:,1)+ orthFoCentSp(:,2)+orthFoCentSp(:,3)*2)/4/8;
recovered4Pulse = (orthFoCentTx(:,1)+ orthFoCentTx(:,2)+orthFoCentTx(:,3)*2)/4/8;
backGound = sum(orthFoCentSp(:,4:10),2)/8/sqrt(7*2.5);%backGound = orthFoCentSp(:,4)/8;
output.recovered4Signal = recovered4Signal;
output.recovered4Pulse = recovered4Pulse;
output.backGound = backGound;

av4Pulse = recovered4Pulse;
av4Signal = recovered4Signal;
avBackGround = backGound;
nValid = 4;
for ii = 2:nValid
    rotatedIdx = circshift((1:length(x)), -(ii-1)*4*fftlSeg);
    av4Pulse = av4Pulse + recovered4Pulse(rotatedIdx);
    av4Signal = av4Signal + recovered4Signal(rotatedIdx);
    avBackGround = avBackGround + backGound(rotatedIdx);
end
av4Pulse = av4Pulse / nValid;
av4Signal = av4Signal / nValid;
avBackGround = avBackGround / sqrt(nValid);
output.av4Signal = av4Signal;
output.av4Pulse = av4Pulse;
output.avBackGround = avBackGround;


if displayInd
    disp(fileInfo.Filename)
    figure;plot(txSig, av4Pulse);grid on
    hold all;plot(txSig, av4Signal);grid on
    hold all;plot(txSig, avBackGround);grid on
    title([fileInfo.Artist ' ' fileInfo.Comment '  syncAv:' num2str(nValid)])
end
output.nValid = nValid;
output.originalData = xOrg;
output.filterResponse = [hwV hw];
output.inputFileName = fullPathName;
output.elapsedTime = toc(startTic);
end

