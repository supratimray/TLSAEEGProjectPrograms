% This program combines analyzedData across subjects. If connMethod is
% empty, it only returns the power values

function dataForDisplay = combineAnalyzedDataConn(folderSourceString,subjectNameList,projectName,refType,protocolType,stRange,freqRanges,connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData)

if ~exist('stRange','var');         stRange = [0.25 0.75];              end
if ~exist('removeMicroSaccadesFlag','var'); removeMicroSaccadesFlag=0;  end
if ~exist('spatialFrequenciesToRemove','var'); spatialFrequenciesToRemove=[];  end
if ~exist('useCleanData','var');    useCleanData=0;                     end

analyzedDataFolder = fullfile(folderSourceString,'analyzedData',projectName,protocolType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numSubjects = length(subjectNameList);

logBLPowerVsFreqAllSubjects = [];
logSTPowerVsFreqAllSubjects = [];
diffPowerAllSubjects = {};
stPowerAllSubjects = {};
blPowerAllSubjects = {};
diffPowerTopoAllSubjects = {};
stPowerTopoAllSubjects = {};
connFreqBandsAllSubjects = {};
connAllSubjects = {}; 
expDatesList = {};

for iSub = 1:numSubjects
    subjectName = subjectNameList{iSub};
    
    % Analysis Input file
    analysisDetailsFile = getAnalysisDetailsFile(analyzedDataFolder,subjectName,refType,protocolType,stRange,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData,[]);
    analysisDetailsFileFreq = [analysisDetailsFile(1:end-4) '_Freq.mat']; % for saving PSDs etc
    
    if ~isfile(analysisDetailsFileFreq)
        disp(['fileName for subject ' subjectName ' does not exist']);
        usableDataFlag = 0;
    else
        analyzedDataFreq = load(analysisDetailsFileFreq);
        usableDataFlag = any(analyzedDataFreq.goodProtFlag);
    end
    
    if ~usableDataFlag
        disp(['No useful protocols for subject ' subjectName]);
    else
        %%%%%%%%%%%%%%%%%%%%%%% PSD and Band Power %%%%%%%%%%%%%%%%%%%%%%%%
        freqVals = analyzedDataFreq.freqPost.freq;
        dataForDisplay.freqVals = freqVals;
        
        % Get PosList
        badFreqPos = getBadFreqPos(freqVals,4);
        numFreqRanges = length(freqRanges);
        posList = cell(1,numFreqRanges);
        for iFR=1:numFreqRanges
            posList{iFR} = setdiff(intersect(find(freqVals>=freqRanges{iFR}(1)),find(freqVals<=freqRanges{iFR}(2))),badFreqPos);
        end
        
        % averaging across electrode Groups
        [expDates,~,capLayout,~] = getProtocolDetailsForAnalysis('ADGammaProject',subjectName,protocolType);
        elecGroupsCell = getElectrodeList(capLayout{1},refType,0,1); % has electrode divisions
        numSides = length(elecGroupsCell)-1;
        
        % converting cell to double array
        elecGroups = cell(1,numSides);
        blPSD = cell(1,numSides);
        stPSD = cell(1,numSides);
        for iG = 1:numSides
            elecGroups{iG} = cell2mat(elecGroupsCell{iG+1}); % Adding 1 because the first one contains the union of all three sides
            blPSD{iG} = log10(squeeze(nanmean(analyzedDataFreq.freqPre.powspctrm(elecGroups{iG},:),1))); % HERE CHSNGE
            stPSD{iG} = log10(squeeze(nanmean(analyzedDataFreq.freqPost.powspctrm(elecGroups{iG},:),1)));
        end
        
        logBLPowerVsFreqAllSubjects = cat(1,logBLPowerVsFreqAllSubjects,blPSD);
        logSTPowerVsFreqAllSubjects = cat(1,logSTPowerVsFreqAllSubjects,stPSD);
        
        numAllElecs = 64;
        diffPower = zeros(numSides,numFreqRanges);
        stPower = zeros(numSides,numFreqRanges);
        blPower = zeros(numSides,numFreqRanges);
        for iSide = 1:numSides
            for j=1:numFreqRanges
                stPow = nanmean(nanmean(analyzedDataFreq.freqPost.powspctrm(elecGroups{iSide},posList{j}),2),1);
                blPow = nanmean(nanmean(analyzedDataFreq.freqPre.powspctrm(elecGroups{iSide},posList{j}),2),1);
                diffPower(iSide,j) = 10*(log10(stPow) - log10(blPow));
                stPower(iSide,j) = (log10(stPow));
                blPower(iSide,j) = (log10(blPow));
            end
        end
        
        stPowerTopo = zeros(numAllElecs,numFreqRanges);
        diffPowerTopo = zeros(numAllElecs,numFreqRanges);
        for j=1:numFreqRanges
            stPowAll = nanmean(analyzedDataFreq.freqPost.powspctrm(:,posList{j}),2);
            blPowAll = nanmean(analyzedDataFreq.freqPre.powspctrm(:,posList{j}),2);
            stPowerTopo(:,j) = (log10(stPowAll));
            diffPowerTopo(:,j) = 10*(log10(stPowAll) - log10(blPowAll));
        end
        
        diffPowerAllSubjects = cat(1,diffPowerAllSubjects,diffPower);
        stPowerAllSubjects = cat(1,stPowerAllSubjects,stPower);
        blPowerAllSubjects = cat(1,blPowerAllSubjects,blPower);
        stPowerTopoAllSubjects = cat(1,stPowerTopoAllSubjects,stPowerTopo);
        diffPowerTopoAllSubjects = cat(1,diffPowerTopoAllSubjects,diffPowerTopo);
        expDatesList = cat(1,expDatesList,expDates{1});
        
        % Compute connectivity measures also if needed
        if ~isempty(connMethod)
            analysisDetailsFileConn = [analysisDetailsFile(1:end-4) '_' connMethod '.mat']; % for saving connectivity '_bl' before '.mat'
            analyzedDataConn = load(analysisDetailsFileConn);
            
            % removing '1' from reference related (x,x) combination
            for elec = 1:numAllElecs
                analyzedDataConn.conn(elec,elec,:) = nan;
            end
            
            connFreqBands = cell(1,numSides);
            connVals = cell(1,numSides);          
            for iSide = 1:numSides
                for j = 1:numFreqRanges
                    connFreqBands{iSide}(j,:,:) = squeeze(nanmean(analyzedDataConn.conn(elecGroups{iSide},:,posList{j}),3));
                end
                connVals{iSide} = analyzedDataConn.conn(elecGroups{iSide},:,:);
            end
            connFreqBandsAllSubjects = cat(1,connFreqBandsAllSubjects,connFreqBands);
            connAllSubjects = cat(1,connAllSubjects,connVals);
        end
    end
end

dataForDisplay.logBLPowerVsFreqAllSubjects = logBLPowerVsFreqAllSubjects;
dataForDisplay.logSTPowerVsFreqAllSubjects = logSTPowerVsFreqAllSubjects;
dataForDisplay.diffPowerAllSubjects = diffPowerAllSubjects;
dataForDisplay.stPowerAllSubjects = stPowerAllSubjects;
dataForDisplay.blPowerAllSubjects = blPowerAllSubjects;
dataForDisplay.diffPowerTopoAllSubjects = diffPowerTopoAllSubjects;
dataForDisplay.stPowerTopoAllSubjects = stPowerTopoAllSubjects;
dataForDisplay.connAllSubjects = connAllSubjects;
dataForDisplay.connFreqBandsAllSubjects = connFreqBandsAllSubjects;
dataForDisplay.expDate = expDatesList;
dataForDisplay.freqVals = freqVals; % adding "freqVals" for further use, directly from "dataForDisplay"
end