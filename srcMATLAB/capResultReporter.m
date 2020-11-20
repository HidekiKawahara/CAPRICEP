function output = capResultReporter(matFileName, outFormat)
%  Visualize and report the CAPRICEP measurement result
%  output = capResultReporter(matFileName, outFormat)
%
%  Argument
%    matFileName      : pathname of the result file
%    outFormat        : output report format:
%                       'EPSF', 'PNG', 'JPEG', 'PDF'
% Output
%    output    : structure variable with the following fields
%           specStr          : copy of the analysis results
%           reportStr        : copy of the report results
%
%
% Licence 
% Copyright 2020 Hideki Kawahara
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

tmp = load(matFileName);
[~,name,ext] = fileparts(matFileName);
specStr = capricepResponseAnalysis(tmp.analysisStrCore);
specStr.dataFile = [name '.' ext];
reportStr = capricepResponseReport(specStr);
output.specStr = specStr;
output.reportStr = reportStr;
switch outFormat
    case 'EPSF'
        print(reportStr.figureHandleReport,'-depsc', ['capRep' name '.eps']);
    case 'PNG'
        print(reportStr.figureHandleReport,'-dpng', '-r200', ['capRep' name '.png']);
    case 'JPEG'
        print(reportStr.figureHandleReport,'-djpeg', '-r200', ['capRep' name '.jpg']);
    case 'PDF'
        print(reportStr.figureHandleReport,'-dpdf', '-bestfit', ['capRep' name '.pdf']);
end
end