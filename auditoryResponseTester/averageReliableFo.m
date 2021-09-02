function avFoOut = averageReliableFo(output)

fs = output.fs;
spFiltered = output.filterOutput(:, 1);
spFilteredLevel = 20*log10(abs(spFiltered));
%nData = length(x);
nData = length(spFilteredLevel);
%figure;
%plot(sort(spFilteredLevel), (1:nData)/nData);grid on;
[sortedLevel, idxSrt] = sort(spFilteredLevel);
voicedLevel = mean(spFilteredLevel(spFilteredLevel > sortedLevel(round(0.8*nData))));
%%
voiced = spFilteredLevel > voicedLevel - 15;
dataIndex = (1:nData)';
onLocation = dataIndex(voiced > voiced([1 1:nData-1]));
offLocation = dataIndex(voiced < voiced([1 1:nData-1]));
if length(onLocation) > 1
    nSegments = length(onLocation);
    segments = zeros(nSegments, 1);
    for ii = 1:nSegments
        segments(ii) = offLocation(ii) - onLocation(ii);
    end
    [~, idxMax] = max(segments);
    onLocation = onLocation(idxMax);
    offLocation = offLocation(idxMax);
elseif onLocation > offLocation
    onLocation = round(fs*0.5);
    offLocation = nData - round(fs*0.5);
else
    onLocation = round(fs*0.5);
    offLocation = nData - round(fs*0.5);
end
indexMargin = round(0.3*fs); % 300ms for eliminating unstable region
spFoRaw = output.rawFoValues;
averageFo = exp(mean(log(spFoRaw(dataIndex > onLocation + indexMargin ...
    & dataIndex < offLocation - indexMargin))));
% check SD in terms of musical cent
spFoCent = 1200 * log2(spFoRaw);
stdFoCent = std(spFoCent(dataIndex > onLocation + indexMargin ...
    & dataIndex < offLocation - indexMargin, :));

%% Check signal RMS level

calibrationCf = output.calibrationConst;
tt = dataIndex / fs;
weightFiltA = weightingFilter('A-weighting',fs);
xc = weightFiltA(output.signal(:,1));
testSignalLevel = 20 * log10(std(xc((dataIndex > indexMargin ...
    & dataIndex < onLocation - indexMargin) ...
    | (dataIndex > offLocation + indexMargin ...
    & dataIndex < nData - indexMargin)))) + calibrationCf;
mixSignalLevel = 20 * log10(std(xc(dataIndex > onLocation + indexMargin ...
    & dataIndex < offLocation - indexMargin))) + calibrationCf;
spSignalLevel = 10 * log10(10^(mixSignalLevel/10) - 10^(testSignalLevel/10));
avFo = mean(output.rawFoi(onLocation:offLocation, :));
%%
avFoOut.averageFo = avFo;
avFoOut.stdFoCent = stdFoCent;
avFoOut.spSignalLevel = spSignalLevel;
avFoOut.outsideLevel = testSignalLevel;
return;

