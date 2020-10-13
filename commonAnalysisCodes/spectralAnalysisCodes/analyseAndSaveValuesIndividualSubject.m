% analyseAndSaveValuesIndividualSubject
% This program calculates and saves PSDs and TF plots that are needed to
% generate all plots 
% Call this function from runAnalyseAndSaveValuesIndividualSubject.m

function analyseAndSaveValuesIndividualSubject(folderSourceString,subjectName,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag)

if ~exist('stRange','var');         stRange = [0.25 0.75];              end
if ~exist('removeMicroSaccadesFlag','var'); removeMicroSaccadesFlag=0;  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(protocolType,'TFCP')
    condVals=16; % SSEVEPFreq/2
else
    condVals=[];
end

dataFolder = fullfile(folderSourceString,'decimatedData',projectName,protocolType);
analyzedDataFolder = fullfile(pwd,'analyzedData',projectName,protocolType); % analysedFolder now in local project directory
makeDirectory(analyzedDataFolder);

% Analysis file
if removeMicroSaccadesFlag
    analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType ...
        '_stRange_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2)) '_NoMS.mat']);
    % numMSInRangePerProtocol ; % TODO
else
    analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType ...
        '_stRange_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2)) '.mat']);
    numMSInRangePerProtocol = [];
end

[expDates,protocolNames,capLayout,usableDataFlag] = getProtocolDetailsForAnalysis(projectName,subjectName,protocolType);

if usableDataFlag && ~isempty(expDates)
    clear fileLists
    for iProt = 1:length(expDates)
        fileLists{iProt} = [subjectName '-' expDates{iProt} '-' protocolNames{iProt} '.mat']; %#ok<AGROW>
    end
    
    % 1. Calculate PSDs and TF plots for selected electrodes
    electrodeList = getElectrodeList(capLayout{1},refType); % Selected electrodes
    [allProtocolsBLData,stPowerVsFreq,blPowerVsFreq,freqVals,tfPower,timeValsTF,freqValsTF,erpData,timeVals,numGoodTrials,numAnalysedElecs]=...
        getDataSingleSubject(dataFolder,fileLists,capLayout,electrodeList,stRange,1,numMSInRangePerProtocol,condVals);
    
    % 2. Calculate PSDs for topoplots - TF data is not stored
    electrodeList = getElectrodeList(capLayout{1},refType,1); % all electrodes
    [allProtocolsBLDataTopo,stPowerVsFreqTopo,blPowerVsFreqTopo,freqValsTopo]=...
        getDataSingleSubject(dataFolder,fileLists,capLayout,electrodeList,stRange,0,numMSInRangePerProtocol,condVals);
    
    save(analysisDetailsFile,'allProtocolsBLData','stPowerVsFreq','blPowerVsFreq',...
        'freqVals','tfPower','timeValsTF','freqValsTF','erpData','timeVals',...
        'numGoodTrials','numAnalysedElecs','allProtocolsBLDataTopo','stPowerVsFreqTopo','blPowerVsFreqTopo','freqValsTopo');
end
end