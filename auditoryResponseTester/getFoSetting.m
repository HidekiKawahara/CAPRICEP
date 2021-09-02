function output = getFoSetting(settingFilePath)
fid = fopen(settingFilePath,'r');
tLine = fgetl(fid);
dataType = 'NONE';
testFoValues = zeros(100, 1);
voicFoValues = zeros(100, 1);
nTest = 0;
nVoice = 0;
FMcent = 25; % defaultValue
inputStr = cell(1);
ii = 0;
while tLine ~= -1
    ii = ii + 1;
    inputStr(ii) = {tLine};
    switch tLine(1)
        case '#'
            dataType = tLine(2:end);
        otherwise
            switch dataType
                case 'TEST_fo'
                    nTest = nTest + 1;
                    testFoValues(nTest) = str2double(tLine);
                case 'VOICE_fo'
                    nVoice = nVoice + 1;
                    voicFoValues(nVoice) = str2double(tLine);
                case 'FMcent'
                    FMcent = str2double(tLine);
            end
    end
    tLine = fgetl(fid);
end
fclose(fid);
output.testFoValues = testFoValues(1:nTest);
output.voicFoValues = voicFoValues(1:nVoice);
output.FMcent = FMcent;
output.inputStr = inputStr;
output.settingFilePath = settingFilePath;
end