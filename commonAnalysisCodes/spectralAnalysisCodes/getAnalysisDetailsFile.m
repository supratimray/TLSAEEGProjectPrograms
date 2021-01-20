function [analysisDetailsFile,numMSInRangePerProtocol] = getAnalysisDetailsFile(analyzedDataFolder,subjectName,refType,protocolType,stRange,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData,temporalFrequencyToUse)

if removeMicroSaccadesFlag
    analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType ...
        '_stRange_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2)) '_NoMS']);
    % numMSInRangePerProtocol ; % TODO
else
    analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType ...
        '_stRange_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2))]);
    numMSInRangePerProtocol = [];
end

if ~isempty(spatialFrequenciesToRemove) && ~strcmp(protocolType,'TFCP')
    analysisDetailsFile = [analysisDetailsFile '_RemoveSF'];
    for i=1:length(spatialFrequenciesToRemove)
        analysisDetailsFile = cat(2,analysisDetailsFile,num2str(spatialFrequenciesToRemove(i)));
    end
end

if useCleanData
    analysisDetailsFile = cat(2,analysisDetailsFile,'_CleanData');
end

if strcmp(protocolType,'TFCP') && (temporalFrequencyToUse==0)
    analysisDetailsFile = cat(2,analysisDetailsFile,'_Static');
end

analysisDetailsFile = [analysisDetailsFile '.mat'];
end