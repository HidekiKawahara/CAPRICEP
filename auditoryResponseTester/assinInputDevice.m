%% test audiorecorder
nAudioIn = audiodevinfo(1);
names = cell(nAudioIn,1);
idList = zeros(nAudioIn,1);
for ii = 1:nAudioIn
    names(ii) = {aaa.input(ii).Name};
    idList(ii) = aaa.input(ii).ID;
end
%%
[idx, tf] = listdlg('ListString', names, 'ListSize', [300, 200] ...
    ,'Name', 'Select input device');
%%
fs = 44100;
recorder = audiorecorder(fs,24,2,idList(idx));
%%
record(recorder,5);
pause(6);
y = getaudiodata(recorder);
stop(recorder);
