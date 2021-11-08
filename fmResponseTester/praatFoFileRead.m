function output = praatFoFileRead(fullpathname)
%  read and convert the pitch file (long format) of Praat
%  This is a quick and dirty implementation
%     output = praatFoFileRead(fullpathname)
%
%  Argument
%     fullpathname  : full pathname of the Praat pitch file
%  Output
%     output  : structure with the following fields
%         n      : number of fo analysis frames
%         fo     : fundamental frequency (Hz)
%         dtx    : interval between frames (s)
%         elapsedTime  : elapsed time of read and conversion (s)

% LICENSE: refer to LICENSE in this folder

startTic = tic;
fid = fopen(fullpathname);
tline = fgetl(fid);
n = 0;
m = 0;
while ~feof(fid)
    if ~isempty(tline)
        aa = textscan(tline, '%s');
        tmp = aa{1};
        switch tmp{1}
            case 'frames'
                fremeID = sscanf(tmp{2}, '[%f]:');
                if ~isempty(fremeID)
                    n = n + 1;
                    fgetl(fid);
                    fgetl(fid);
                    fgetl(fid);
                    fgetl(fid);
                    tline = fgetl(fid);
                    aa = textscan(tline, '%s');
                    tmp = aa{1};
                    fo(n) = str2double(tmp{3});
                end
            case 'nx'
                disp(tmp{3})
                nFrame = str2double(tmp{3});
                fo = zeros(nFrame, 1);
            case 'dx'
                dtx = str2double(tmp{3});
        end
    end
    tline = fgetl(fid);
end
fclose(fid);
output.n = n;
output.fo = fo;
output.dtx = dtx;
output.elapsedTime = toc(startTic);
end
