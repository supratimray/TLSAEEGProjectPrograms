% This program displays the data saved in the analyzedData folder

% subjectNameList is a cell array, with each cell containing a list of
% subjects. For example, subjectNameList{1} could contain all MCI/AD subjects
% while subjectNameList{2} could contain all their controls.

function displayAnalyzedDataConn(folderSourceString,subjectNameLists,methodOptions,strList,projectName,refType,protocolType,stRange,freqRanges,freqRangeNames,removeMicroSaccadesFlag,useMedianFlag,spatialFrequenciesToRemove,useCleanData)

if ~exist('stRange','var');         stRange = [0.25 0.75];              end
if ~exist('freqRanges','var')
    freqRanges{1} = [8 12]; freqRangeNames{1} = 'Alpha'; % alpha
    freqRanges{2} = [20 34]; freqRangeNames{2} = 'Slow gamma'; % slow gamma
    freqRanges{3} = [36 66]; freqRangeNames{3} = 'Fast gamma'; % fast gamma
end
if ~exist('removeMicroSaccadesFlag','var'); removeMicroSaccadesFlag=0;  end
if ~exist('useMedianFlag','var');   useMedianFlag = 1;                  end
if ~exist('spatialFrequenciesToRemove','var'); spatialFrequenciesToRemove=[];  end
if ~exist('useCleanData','var');    useCleanData=0;                     end

numGroups = length(subjectNameLists);
numFreqRanges = length(freqRanges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataForDisplayAllGroups = cell(numGroups,numFreqRanges);

for i=1:numFreqRanges % Analysis done separately for each frequency
 
    % First, get subjectList for all groups
    if strcmp(methodOptions.comparisonType,'MidVsOld')
        
        if strcmp(methodOptions.controlType,'unmatched')
            subjectNameListTMP{1} = subjectNameLists{1};
            subjectNameListTMP{2} = subjectNameLists{2};
        elseif strcmp(methodOptions.controlType,'matched')
            subjectNameListTMP = getPowerMatchedSubjectLists(folderSourceString,subjectNameLists,projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.sideToShow,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
        end
        
        dataForDisplayAllGroups{1,i} = combineAnalyzedDataConn(folderSourceString,subjectNameListTMP{1},projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
        dataForDisplayAllGroups{2,i} = combineAnalyzedDataConn(folderSourceString,subjectNameListTMP{2},projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
        
    elseif strcmp(methodOptions.comparisonType,'CaseVsControl')
    
        numCases = length(subjectNameLists{2});
        
        dataForDisplayAllGroupsTMP = cell(numGroups,numCases);
        
        for iCase=1:numCases
            dataForDisplayAllGroupsTMP{2,iCase} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{2}{iCase},projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
                
            if strcmp(methodOptions.controlType,'unmatched') % For each case, simply average all controls - note that this is like the matched case in the ADGammaProject    
                dataForDisplayAllGroupsTMP{1,iCase} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1}{iCase},projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
                
            elseif strcmp(methodOptions.controlType,'matched')
                tmp = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1}{iCase},projectName,refType,protocolType,stRange,freqRanges(i),[],removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
                
                % Match power in tmp and dataForDisplayAllGroupsTMP{2,iCase}
            end
        end
    end
end

% Display Data
chanlocs = getMontageDetails(refType);
end
function chanlocs = getMontageDetails(refType)

capLayout = 'actiCap64';
clear cL bL chanlocs iElec electrodeList noseDir
switch refType
    case 'unipolar'
        cL = load([capLayout '.mat']);
        chanlocs = cL.chanlocs;
    case 'bipolar'
        cL = load(['bipolarChanlocs' capLayout '.mat']);
        chanlocs = cL.eloc;
end
end
function matchedSubjectNameLists = getPowerMatchedSubjectLists(folderSourceString,subjectNameLists,projectName,refType,protocolType,stRange,freqRange,sideToShow,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData)

powerMatchingBinWidth = 0.5; % dB

dataForDisplay{1} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1},projectName,refType,protocolType,stRange,freqRange,[],removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
dataForDisplay{2} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{2},projectName,refType,protocolType,stRange,freqRange,[],removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);

tmp = cell2mat(dataForDisplay{1}.diffPowerAllSubjects'); dataToMatch{1} = tmp(sideToShow,:);
tmp = cell2mat(dataForDisplay{2}.diffPowerAllSubjects'); dataToMatch{2} = tmp(sideToShow,:);

minVal = min(min(dataToMatch{1}),min(dataToMatch{2}));
maxVal = max(max(dataToMatch{1}),max(dataToMatch{2}));

powerBins = minVal:powerMatchingBinWidth:maxVal;
numPowerBins = length(powerBins);
d = powerBins(2)-powerBins(1);
powerBins = [powerBins powerBins(numPowerBins)+d];

matchedSubjectNameLists{1} = []; 
matchedSubjectNameLists{2} = [];
for i=1:numPowerBins
    pos1 = intersect(find(dataToMatch{1}>=powerBins(i)),find(dataToMatch{1}<powerBins(i+1)));
    pos2 = intersect(find(dataToMatch{2}>=powerBins(i)),find(dataToMatch{2}<powerBins(i+1)));
    
    [equalPos1,equalPos2] = getEqualNumOfIndices(pos1,pos2);
    matchedSubjectNameLists{1} = cat(2,matchedSubjectNameLists{1},subjectNameLists{1}(equalPos1));
    matchedSubjectNameLists{2} = cat(2,matchedSubjectNameLists{2},subjectNameLists{2}(equalPos2));
end
end
function [x2,y2] = getEqualNumOfIndices(x1,y1)

N1 = length(x1);
N2 = length(y1);

if (N1==0) || (N2==0) % one of the two is an empty array
    x2=[]; y2=[];
    
elseif N1==N2
    x2=x1; y2=y1;
    
elseif N1<N2
    x2 = x1;
    randVals = randperm(N2);
    y2 = y1(sort(randVals(1:N1)));
    
else %N1>N2
    y2 = y1;
    randVals = randperm(N1);
    x2 = x1(sort(randVals(1:N2)));
end
end