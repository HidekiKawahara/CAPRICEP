%% chromatic scale table

ao = 27.5;
noteHead = {'A','A','B','C','C','D','D','E','F','F','G','G'};
noteMod = {' ','#',' ',' ','#',' ','#',' ',' ','#',' ', '#'};
fo = ao;
foList = zeros(6*12, 1);
noteNameList = zeros(6*12,3);
noteId = 0;
for ii = 0:5
    for jj = 1:12
        if jj >= 4
            tmp = 1;
        else
            tmp = 0;
        end
        noteId = noteId + 1;
        foList(noteId) = fo;
        fo = fo * 2^(1/12);
        noteNameList(noteId, :) = [noteHead{jj} num2str(ii + tmp) noteMod{jj}];
        disp([num2str(foList(noteId), '%6.3f') '   ' noteNameList(noteId, :)]);
    end
end