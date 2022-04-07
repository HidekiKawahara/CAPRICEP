%%

load testSigBase400ms.mat
txAll = (1:fftl)'/fs - fftl/2/fs;
%%
lloc = (1:fftlSeg:fftl)'/fs - fftl/2/fs;
figure;
set(gcf, 'Position',[3433         557         857         371])
tiledlayout(2,1,"TileSpacing","tight")
nexttile
plot(txAll, xTSPSel(:,1),'k',"LineWidth",2);grid on;
for ii = 1:length(lloc)
    xline(lloc(ii),'-.','LineWidth',2)
end
set(gca,'LineWidth',2,'FontSize',16)
axis([-3 3 -0.04 0.04])
nexttile
plot(txAll, 20*log10(abs(xTSPSel(:,1))),'k',"LineWidth",2);grid on;
set(gca,'LineWidth',2,'FontSize',16)
axis([-3 3 -300 0])
xlabel('time (s)')
ylabel('level (dB)')
for ii = 1:length(lloc)
    xline(lloc(ii),'-.','LineWidth',2)
end
print -deps Figure8.eps

%%

fs = output.testSignalData.fs;
fftl = output.testSignalData.fftl;
fftlSeg = output.testSignalData.fftlSeg;
xTSPSel = output.testSigOut.xTSPSel;
combination = output.testSigOut.combination;
orthSeq = output.testSigOut.orthSeq;
xTR = xTSPSel(:,combination);
xMix = sum(orthSeq,2);
tt = (1:length(xMix))/fs;
figure;
plot(tt, xMix);grid on;
%%
q = fftfilt(xTR(end:-1:1,:),[xMix;zeros(fftl,1)]);
q = q(fftl/2+(1:length(tt)),:);
ddata = length(tt);
baseid = (1:ddata);
idx1 = max(1, min(ddata, baseid+fftlSeg+fftlSeg/2));
idx2 = max(1, min(ddata, baseid+0*fftlSeg+fftlSeg/2));
idx3 = max(1, min(ddata, baseid-fftlSeg+fftlSeg/2));
idx4 = max(1, min(ddata, baseid-2*fftlSeg+fftlSeg/2));
c1 = q(idx1,1) + q(idx2,1) + q(idx3,1) + q(idx4,1);
c2 = q(idx1,2) - q(idx2,2) + q(idx3,2) - q(idx4,2);
c3 = q(idx1,3) - q(idx2,3) - q(idx3,3) + q(idx4,3);
%%
figure;
set(gcf, 'Position',[3433          53         887         647])
tiledlayout(4,1,"TileSpacing","tight")
nexttile;
plot(tt, xMix,'k','LineWidth',2);grid on;
set(gca,'LineWidth',2,'FontSize',16,'XLim',[0 20],'YLim',[-0.05 0.05])
nexttile;
plot(tt, 20*log10(abs(q(:,1))),'k','LineWidth',2);grid on;
set(gca,'LineWidth',2,'FontSize',16,'XLim',[0 20],'YLim',[-300 0.05])
ylabel('level (rel. dB)')
nexttile;
plot(tt,  20*log10(abs(c1/4)),'k','LineWidth',2);grid on;
set(gca,'LineWidth',2,'FontSize',16,'XLim',[0 20],'YLim',[-300 0.05])
ylabel('level (rel. dB)')
nexttile;
plot(tt, 20*log10(abs(c2+c1+c3*2)/16),'k','LineWidth',2);grid on;
set(gca,'LineWidth',2,'FontSize',16,'XLim',[0 20],'YLim',[-300 0.05])
ylabel('level (rel. dB)')
xlabel('time (s)')
print -depsc Figure9.eps
