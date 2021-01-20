% This program combines analyzedData across subjects

function dataForDisplay = combineAnalyzedData(folderSourceString,subjectNameList,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag,gamma1Range,gamma2Range,alphaRange,spatialFrequenciesToRemove,useCleanData,temporalFrequencyToUse)

if ~exist('stRange','var');         stRange = [0.25 0.75];              end
if ~exist('removeMicroSaccadesFlag','var'); removeMicroSaccadesFlag=0;  end
if ~exist('gamma1Range','var');     gamma1Range = [20 34];              end
if ~exist('gamma2Range','var');     gamma2Range = [36 66];              end
if ~exist('alphaRange','var');      alphaRange = [8 12];                end
if ~exist('spatialFrequenciesToRemove','var'); spatialFrequenciesToRemove=[];  end
if ~exist('useCleanData','var');    useCleanData=0;                     end
if ~exist('temporalFrequencyToUse','var'); temporalFrequencyToUse=[];   end

SSVEPFreqHz = 2*temporalFrequencyToUse;

analyzedDataFolder = fullfile(folderSourceString,'analyzedData',projectName,protocolType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numSubjects = length(subjectNameList);

erpDataAllSubjects = [];
logBLPowerVsFreqAllSubjects = [];
logSTPowerVsFreqAllSubjects = [];
powerDBAllSubjects = [];
dTFPowerDBAllSubjects = [];
powerDBTopoAllSubjects = [];

for iSub = 1:numSubjects
    subjectName = subjectNameList{iSub};
    
    % Analysis Input file
    analysisDetailsInputFile = getAnalysisDetailsFile(analyzedDataFolder,subjectName,refType,protocolType,stRange,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData,temporalFrequencyToUse);
    
    if ~isfile(analysisDetailsInputFile)
        disp(['fileName for subject ' subjectName ' does not exist']);
        usableDataFlag = 0;
    else
        analyzedData = load(analysisDetailsInputFile);
        usableDataFlag = any(analyzedData.allProtocolsBLData.goodProtFlag);
    end
    
    if ~usableDataFlag
        disp(['No useful protocols for subject ' subjectName]);
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~exist('timeVals','var')
            timeVals = analyzedData.timeVals;
            dataForDisplay.timeVals = timeVals;
        elseif ~isequal(timeVals,analyzedData.timeVals)
            error(['timeVals do not match for ' subjectName]);
        end
        erpDataAllSubjects = cat(1,erpDataAllSubjects,squeeze(mean(analyzedData.erpData,1)));
        
        %%%%%%%%%%%%%%%%%%%%%%% PSD and Band Power %%%%%%%%%%%%%%%%%%%%%%%%
        if ~exist('freqVals','var')
            freqVals = analyzedData.freqVals;
            dataForDisplay.freqVals = freqVals;
            
            if strcmpi(protocolType,'SF_ORI') || (strcmpi(protocolType,'TFCP') && (temporalFrequencyToUse==0))
                % Get alpha and gamma Pos
                posList = cell(1,3);
                rangeNames{1} = 'SG';
                rangeNames{2} = 'FG';
                rangeNames{3} = 'Alpha';
                dataForDisplay.rangeNames = rangeNames;
                
                badFreqPos = getBadFreqPos(freqVals);
                posList{1} = setdiff(intersect(find(freqVals>=gamma1Range(1)),find(freqVals<=gamma1Range(2))),badFreqPos);
                posList{2} = setdiff(intersect(find(freqVals>=gamma2Range(1)),find(freqVals<=gamma2Range(2))),badFreqPos);
                posList{3} = setdiff(intersect(find(freqVals>=alphaRange(1)),find(freqVals<=alphaRange(2))),badFreqPos);
            elseif strcmp(protocolType,'TFCP')
                rangeNames{1} = 'SSVEP';
                dataForDisplay.rangeNames = rangeNames;
                posList{1} = find(freqVals==SSVEPFreqHz);
            end
            
        elseif ~isequal(freqVals,analyzedData.freqVals)
            error(['freqVals do not match for ' subjectName]);
        end
        
        blPSD = squeeze(mean(analyzedData.blPowerVsFreq,1));
        stPSD = squeeze(mean(analyzedData.stPowerVsFreq,1));
        logBLPowerVsFreqAllSubjects = cat(1,logBLPowerVsFreqAllSubjects,log10(blPSD));
        logSTPowerVsFreqAllSubjects = cat(1,logSTPowerVsFreqAllSubjects,log10(stPSD));
        
        numPosList = length(posList);
        powerDB = zeros(1,numPosList);
        for j=1:numPosList
            powerDB(j) = 10*(log10(mean(stPSD(posList{j}))) - log10(mean(blPSD(posList{j}))));
        end
        powerDBAllSubjects = cat(1,powerDBAllSubjects,powerDB);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Time-Frequency plots %%%%%%%%%%%%%%%%%%
        if ~exist('timeValsTF','var')
            timeValsTF = analyzedData.timeValsTF;
            dataForDisplay.timeValsTF = timeValsTF;
            
            blRange = [-diff(stRange) 0];
            blPosTF = timeValsTF>=blRange(1) & timeValsTF<=blRange(2);
        elseif ~isequal(timeValsTF,analyzedData.timeValsTF)
            error(['timeValsTF do not match for ' subjectName]);
        end
        
        if ~exist('freqValsTF','var')
            freqValsTF = analyzedData.freqValsTF;
            dataForDisplay.freqValsTF = freqValsTF;
        elseif ~isequal(freqValsTF,analyzedData.freqValsTF)
            error(['freqValsTF do not match for ' subjectName]);
        end
        
        tfPower = squeeze(mean(analyzedData.tfPower,1));
        tfPowerBL = repmat(squeeze(mean(tfPower(blPosTF,:),1)),length(timeValsTF),1);
        dTFPower = 10.*log10(tfPower./tfPowerBL);
        dTFPowerDBAllSubjects = cat(3,dTFPowerDBAllSubjects,dTFPower);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isequal(freqVals,analyzedData.freqValsTopo)
            error(['freqValsTopo and frqVals dont match for subject ' subjectName '. Use different posLists for topoplot...']);
        else
            powerDBTopo = getDiffPowerSingleSubject(analyzedData.stPowerVsFreqTopo,analyzedData.blPowerVsFreqTopo,posList);
            if size(powerDBTopo,1)==50 % actiCap31Posterior
                % disp(['actiCap31Posterior was used for subject ' subjectName '. Readjusting topoplot...']);
                if ~exist('bipInds31To64','var')
                    bipInds31To64 = load('bipInds31To64.mat'); % Available in the montages folder
                    bipInds31To64 = bipInds31To64.bipInds31To64;
                end
                powerDBTopoLong = nan(112,numPosList);
                powerDBTopoLong(bipInds31To64,:) = powerDBTopo;
                powerDBTopoAllSubjects = cat(3,powerDBTopoAllSubjects,powerDBTopoLong);
            else
                powerDBTopoAllSubjects = cat(3,powerDBTopoAllSubjects,powerDBTopo);
            end
        end
    end
end

dataForDisplay.erpData = erpDataAllSubjects;
dataForDisplay.logBLPowerVsFreqAllSubjects = logBLPowerVsFreqAllSubjects;
dataForDisplay.logSTPowerVsFreqAllSubjects = logSTPowerVsFreqAllSubjects;
dataForDisplay.powerDBAllSubjects = powerDBAllSubjects;
dataForDisplay.dTFPowerDBAllSubjects = dTFPowerDBAllSubjects;
dataForDisplay.powerDBTopoAllSubjects = powerDBTopoAllSubjects;
end
