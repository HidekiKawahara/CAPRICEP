%%
%-rw-r--r--  1 kawahara  staff   924852712  3 19 03:15 dataOutPEFr.mat
%-rw-r--r--  1 kawahara  staff   925329643  3 19 05:54 dataOutCREPEr.mat
%-rw-r--r--  1 kawahara  staff   945263067  3 19 11:06 dataOutNDF800.mat
%-rw-r--r--  1 kawahara  staff   926212731  3 19 13:34 dataOutPraatR.mat
%-rw-r--r--  1 kawahara  staff   925593954  3 23 19:58 dataOutREAPER.mat
%-rw-r--r--  1 kawahara  staff   924526149  3 23 21:34 dataOutNCF.mat
%-rw-r--r--  1 kawahara  staff   924830810  3 23 21:35 dataOutOpenSmile.mat
%-rw-r--r--  1 kawahara  staff   924976847  3 23 21:36 dataOutOpenSmiles.mat
%-rw-r--r--  1 kawahara  staff  1114743262  3 23 21:37 dataOutNINJALX2.mat
%-rw-r--r--  1 kawahara  staff   915864304  3 24 02:51 dataOutfxRAPT.mat
%-rw-r--r--  1 kawahara  staff   924736994  3 25 02:02 dataOutCEP.mat
%-rw-r--r--  1 kawahara  staff   924649004  3 25 02:16 dataOutLHS.mat
%-rw-r--r--  1 kawahara  staff   924841515  3 25 02:25 dataOutSRH.mat
%-rw-r--r--  1 kawahara  staff  1107817307  3 25 03:17 dataNINJAL.mat
%-rw-r--r--  1 kawahara  staff   928437226  3 25 15:28 dataXSX.mat
%-rw-r--r--  1 kawahara  staff   928358998  3 25 16:50 dataHarves.mat
%-rw-r--r--  1 kawahara  staff   924604754  3 26 02:00 dataSWIPEP.mat
%%
fmCent = 20;
fl = 80;
fh = 400;
nOct = 48;
dataOutYIN = makePitchTestSetSlim(fl, fh, nOct, fmCent, @pitchYIN_estimate);

%%
charCellArray = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U'};
dataStr = struct;
dataStr(1).name = "dataOutNCF";
dataStr(2).name = "dataOutCEP";
dataStr(3).name = "dataOutLHS";
dataStr(4).name = "dataOutSRH";
dataStr(5).name = "dataOutPEFr";
dataStr(6).name = "dataSWIPEP";
dataStr(7).name = "dataOutYINe";
dataStr(8).name = "dataOutOpenSmile";
dataStr(9).name = "dataOutOpenSmiles";
dataStr(10).name = "dataOutCREPEr";
dataStr(11).name = "dataOutPraatR";
dataStr(12).name = "dataOutREAPER";
dataStr(13).name = "dataOutfxRAPT";
dataStr(14).name = "dataXSX";
dataStr(15).name = "dataOutNDF800";
dataStr(16).name = "dataHarves";
dataStr(17).name = "dataNINJAL";
dataStr(18).name = "dataOutNINJALX2";
dataStr(19).name = "dataOutYANG";
dataStr(20).name = "dataOutcSRH";
nData = length(dataStr);
fqVar = zeros(nData, 1);
foVar = zeros(nData, 1);
figure;
for ii = 19:20%1:nData
    switch ii
        case 7
            frameStr = dataOutYIN.frameStr;
        case 8
            tmp = load(dataStr(ii).name);
            frameStr = tmp.dataOutOpSMILE.frameStr;
        case 9
            tmp = load(dataStr(ii).name);
            frameStr = tmp.dataOutOpSMILEs.frameStr;
        otherwise
            tmp = load(dataStr(ii).name);
            eval("frameStr = tmp." + dataStr(ii).name + ".frameStr;");
    end
    nFreq = length(frameStr);
    fx = frameStr(1).fxlResp;
    nBins = sum((fx > 1 & fx <= 10));
    respMap = zeros(nBins,nFreq);
    for jj = 1:nFreq
        ltiGain = 20*log10(abs(frameStr(jj).averageResponseL));
        respMap(:,jj) = ltiGain(fx > 1 & fx <= 10);
        respMap(isinf(respMap(:,jj)),jj) = 0;
    end
    ggg = diff(respMap);
    aaa = mean(respMap);
    fqVar(ii) = std(ggg(:)) / fx(2);
    foVar(ii) = std(diff(aaa)) * 4;
    loglog(fqVar(ii), foVar(ii), 'ko', 'LineWidth', 1,'MarkerSize',16);
    text(fqVar(ii), foVar(ii), charCellArray{ii}, 'fontsize', 14,'HorizontalAlignment','center','VerticalAlignment','middle');
    set(gca, 'FontSize',16,'LineWidth',2)
    xlabel('SD on modulation frequency change (dB/Hz)')
    ylabel('SD on pitch change (dB/semitone)')
    axis([0.03 4 0.01 20])
    grid on;
    hold all
    drawnow
end
hold off
%print -deps gainSDonFqandPitch.eps
%%
avBw = zeros(nData, 1);
avSNR = zeros(nData, 1);
figure;
for ii = 1:nData
    switch ii
        case 7
            frameStr = dataOutYIN.frameStr;
        case 8
            tmp = load(dataStr(ii).name);
            frameStr = tmp.dataOutOpSMILE.frameStr;
        case 9
            tmp = load(dataStr(ii).name);
            frameStr = tmp.dataOutOpSMILEs.frameStr;
        otherwise
            tmp = load(dataStr(ii).name);
            eval("frameStr = tmp." + dataStr(ii).name + ".frameStr;");
    end
    nFreq = length(frameStr);
    fx = frameStr(1).fxlResp;
    fxSeg = frameStr(1).fxSeg;
    nBins = sum((fx > 1 & fx <= 10));
    bwSet = zeros(nFreq,1);
    snrSet = zeros(nFreq,1);
    nCount = 0;
    for jj = 1:nFreq
        varNLTI = frameStr(jj).varNLTI;
        varNLTI(1) = varNLTI(2);
        varNLTIL = interp1(fxSeg,varNLTI,fx, "linear","extrap");
        totalErr = varNLTIL + frameStr(jj).varTVL;
        snr = 20*log10(abs(frameStr(jj).averageResponseL))-10*log10(totalErr);
        fxIntLim = min(fx(snr < 0 & fx >0));
        if ~isempty(fxIntLim)
            pwspec = abs(frameStr(jj).averageResponseL) .^2;
            bwSet(jj) = sqrt(sum(pwspec(fx<fxIntLim) .* fx(fx<fxIntLim).^2) ...
                / sum(pwspec(fx<fxIntLim)));
            snrSet(jj) = 10*log10(sum(pwspec(fx<bwSet(jj)))) ...
                -10*log10(sum(totalErr(fx<bwSet(jj))));
            if isnan(bwSet(jj)) || isnan(snrSet(jj))
                bwSet(jj) = 0;
                snrSet(jj) = 0;
            else
                nCount = nCount + 1;
            end
        end
    end
    avBw(ii) = sum(bwSet)/nCount;
    avSNR(ii) = sum(snrSet)/nCount;
    plot(avBw(ii), avSNR(ii), 'ko', 'LineWidth', 1,'MarkerSize',16);
    text(avBw(ii), avSNR(ii), charCellArray{ii}, 'fontsize', 14,'HorizontalAlignment','center','VerticalAlignment','middle');
    set(gca, 'FontSize',16,'LineWidth',2)
    xlabel('modulation frequency bandwidth (Hz)')
    ylabel('signal to noise ratio (dB)')
    %axis([0.03 4 0.01 20])
    grid on;
    hold all
    drawnow
end
hold off
print -deps linaverageBwSNR.eps
%%
%%
avBw = zeros(nData, 1);
avSNR = zeros(nData, 1);
figure;
set(gcf, 'Position', [3791        -130        1270        1019])
for ii = 1:nData
    switch ii
        case 7
            frameStr = dataOutYIN.frameStr;
        case 8
            tmp = load(dataStr(ii).name);
            frameStr = tmp.dataOutOpSMILE.frameStr;
        case 9
            tmp = load(dataStr(ii).name);
            frameStr = tmp.dataOutOpSMILEs.frameStr;
        otherwise
            tmp = load(dataStr(ii).name);
            eval("frameStr = tmp." + dataStr(ii).name + ".frameStr;");
    end
    nFreq = length(frameStr);
    fx = frameStr(1).fxlResp;
    fxSeg = frameStr(1).fxSeg;
    nBins = sum((fx > 1 & fx <= 10));
    bwSet = ones(nFreq,1);
    snrSet = zeros(nFreq,1);
    nCount = 0;
    for jj = 1:nFreq
        varNLTI = frameStr(jj).varNLTI;
        varNLTI(1) = varNLTI(2);
        varNLTIL = interp1(fxSeg,varNLTI,fx, "linear","extrap");
        totalErr = varNLTIL + frameStr(jj).varTVL;
        snr = 20*log10(abs(frameStr(jj).averageResponseL))-10*log10(totalErr);
        fxIntLim = min(fx(snr < 0 & fx >0));
        if ~isempty(fxIntLim)
            pwspec = abs(frameStr(jj).averageResponseL) .^2;
            bwSet(jj) = sqrt(sum(pwspec(fx<fxIntLim) .* fx(fx<fxIntLim).^2) ...
                / sum(pwspec(fx<fxIntLim)));
            snrSet(jj) = 10*log10(sum(pwspec(fx<bwSet(jj)))) ...
                -10*log10(sum(totalErr(fx<bwSet(jj))));
            if isnan(bwSet(jj)) || isnan(snrSet(jj)) || bwSet(jj) < 1
                bwSet(jj) = 1;
                snrSet(jj) = 0;
            else
                nCount = nCount + 1;
            end
        end
    end
    avBw(ii) = exp(sum(log(bwSet))/nCount);
    avSNR(ii) = sum(snrSet)/nCount;
    semilogx(avBw(ii), avSNR(ii), 'ko', 'LineWidth', 1,'MarkerSize',32);
    text(avBw(ii), avSNR(ii), charCellArray{ii}, 'fontsize', 28,'HorizontalAlignment','center','VerticalAlignment','middle');
    set(gca, 'FontSize',28,'LineWidth',2)
    xlabel('modulation frequency bandwidth (Hz)')
    ylabel('signal to noise ratio (dB)')
    %axis([0.03 4 0.01 20])
    grid on;
    hold all
    drawnow
end
hold off
set(gca,'XMinorGrid','off')
set(gca, 'xtick',[2 3 4 5 6 7 8 9 10 20 30],'xticklabel',{'2','3','4','5','6','7','8','9','10','20','30'});
axis([3 30 0 60])
%
nameList = {'A NCF(MATLAB)';'B CEP(MATLAB)';'C LHS(MATLAB)';'D SRH(MATLAB)';'E PEF(MATLAB)';'F SWIPEP';'G YIN estimate'; ...
    'H ARC(openSMILE)';'I SHS(openSMILE)';'J CREPE';'K Praat';'L REAPER';'M RAPT(VOICEBOX)'; ...
    'N XSX (TANDEM-STRAIGHT)';'O NDF (legacy-STRAIGHT)';'P Harvest (WORLD)';'Q NINJAL';'R NINJALX2';'S YANG';'T SRH(COVAREP)'};
%%
legend(nameList,'Location','northoutside','Orientation','horizontal', ...
    'NumColumns',4, 'FontSize',20);
print -deps logaverageBwSNRMag2.eps
%%
nData = length(dataStr);
fqVar = zeros(nData, 1);
foVar = zeros(nData, 1);
figure;
set(gcf, 'Position', [3791        -130        1270        1019])
for ii = 1:nData
    switch ii
        case 7
            frameStr = dataOutYIN.frameStr;
        case 8
            tmp = load(dataStr(ii).name);
            frameStr = tmp.dataOutOpSMILE.frameStr;
        case 9
            tmp = load(dataStr(ii).name);
            frameStr = tmp.dataOutOpSMILEs.frameStr;
        otherwise
            tmp = load(dataStr(ii).name);
            eval("frameStr = tmp." + dataStr(ii).name + ".frameStr;");
    end
    nFreq = length(frameStr);
    fx = frameStr(1).fxlResp;
    nBins = sum((fx > 1 & fx <= 10));
    respMap = zeros(nBins,nFreq);
    for jj = 1:nFreq
        ltiGain = 20*log10(abs(frameStr(jj).averageResponseL));
        respMap(:,jj) = ltiGain(fx > 1 & fx <= 10);
        respMap(isinf(respMap(:,jj)),jj) = 0;
    end
    ggg = diff(respMap);
    aaa = mean(respMap);
    fqVar(ii) = std(ggg(:)) / fx(2);
    foVar(ii) = std(diff(aaa)) * 4;
    loglog(fqVar(ii), foVar(ii), 'ko', 'LineWidth', 1,'MarkerSize',32);
    text(fqVar(ii), foVar(ii), charCellArray{ii}, 'fontsize', 28,'HorizontalAlignment','center','VerticalAlignment','middle');
    set(gca, 'FontSize',28,'LineWidth',2)
    xlabel('SD on modulation frequency change (dB/Hz)')
    ylabel('SD on pitch change (dB/semitone)')
    axis([0.03 4 0.01 20])
    grid on;
    hold all
    drawnow
end
legend(nameList,'Location','northoutside','Orientation','horizontal', ...
    'NumColumns',4, 'FontSize',20);
hold off
print -deps gainSDonFqandPitchMag2.eps
