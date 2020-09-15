% Needs to be fixed

function getMicrosaccadeData

% Extract microsaccades. % MD: 02-07-2019; 
% MD: modified 01-July-2020: removed option to save eye-data here, as it is saved in decimatedData file already.
if mainFlags.eyeDataFlag
try    
    clear saveFolder
    if strcmpi(protocolType,'CON'); protType = 'Contrast'; else protType = protocolType; end  %#ok<NASGU>
    saveFolder = eval(giveRawAnalysisSaveString); 
    
    analysisDetails = load(analysisDetailsFile,'gammaPowerValsAreaWise','allProtocolsBLData');
    allProtocolsBLData = analysisDetails.allProtocolsBLData;  
    goodProtVals = getDetails(analysisDetails,'allProtocolsBLData',subjectName);
    goodProtFlag = goodProtVals{2};
    
    if sum(goodProtFlag)>0
        clear fileLists
        for iProt = 1:length(expDates)
            fileLists{iProt} = [subjectName '-' expDates{iProt} '-' protocolNames{iProt} '.mat'];
        end
        fileLists(~goodProtFlag)=[];
        
        clear eyeDataDetails; eyeDataDetails = load(fullfile(cleanDataFolder,fileLists{1}),'eyeRangeMS');
        eyeRangeMS = eyeDataDetails.eyeRangeMS;
        if ~isempty(eyeRangeMS)
                        
            % Get eye data
            [eyeDataDegX,eyeDataDegY,~,~,trialNums,protNumForTrial,FsEyes] = ...
                getEyeDataIndividualSubject(fileLists,cleanDataFolder);
            
            % Save microsaccades
            if ~exist(fullfile(saveFolder,'eyeData',[subjectName '.fig']),'file') && MSFlag
                clear figH trialsForMSAnalysis
                figH = displayMicrosaccades_v2(fullfile(saveFolder,'eyeData'),eyeDataDegX,eyeDataDegY,...
                    trialNums,protNumForTrial,subjectName,protocolType,FsEyes,eyeRangeMS,1); drawnow;
                savefig(figH,fullfile(saveFolder,'eyeData',[subjectName '.fig']));
                close(figH);
            end
    
            % The results from displayMicrosaccades_v2 are stored in
            % separate folders for each subject which could be used as temporary locations.
            % Collating them into a single file here for storage and
            % distribution in a common folder.
            % Save data for SF_ORI_New also in SF_ORI folder for ease of
            % storage and access 
            microsaccades = load(fullfile(saveFolder,'eyeData','microsaccades.mat')); %#ok<NASGU>
            microsaccadesStats = load(fullfile(saveFolder,'eyeData','microsaccadesStats.mat')); %#ok<NASGU>
            save(fullfile(analysedEyeDataFolder,'SF_ORI',[subjectName '_AnalysedEyeData.mat']),'microsaccades','microsaccadesStats');
        else
            disp('Eye data absent')
        end
    else
        disp('No good protocols for MS calculation')
    end
    
catch err2
    disp(['............' err2.message]); 
    disp('Could not calculate MS data');
end

end

%%%%%%%%%%%%%%%%%%%%%%
% Calculate gamma power for trials containing no microsaccades
%%{
if mainFlags.noMSTrialsFlag
try
    AnalysisTypesSF_Ori_BipolarNoMS = {'gammaPowerValsAreaWiseNoMS' 'alphaPowerValsAreaWiseNoMS' 'allProtocolsBLData'};% 'alphaPowerValsAreaWiseNoMS' 'tfPower_sfAreaWiseNoMS' 'baseline_SFAreaWiseNoMS' 'logSTPowerVsFreq_SFAreaWiseNoMS'};
    clear analysisDetails; analysisDetails = load(analysisDetailsFile,AnalysisTypesSF_Ori_BipolarNoMS{:});
    analysisFlagSF_Ori_BipolarNoMS = checkForAnalysis(analysisDetails,subjectName,AnalysisTypesSF_Ori_BipolarNoMS);

%     analysisFlagSF_Ori_BipolarNoMS = 1;
    if analysisFlagSF_Ori_BipolarNoMS
        
        % temporary folder
        clear protType saveFolder
        if strcmpi(protocolType,'CON'); protType = 'Contrast'; else protType = protocolType; end %#ok<NASGU>
        saveFolder = eval(giveRawAnalysisSaveString); rawAnalysisFile = ['rawAnalysisNoMS_' refType '.mat'];            
        
        clear stPowerVsFreqSF_Ori blPowerVsFreqSF_Ori freqValsSF_Ori tfPowerSF_Ori timeValsTF freqValsTF gamma1Pos gamma2Pos numGoodTrialsSF_Ori erpData timeVals 
        try
            load(fullfile(saveFolder,rawAnalysisFile));
        catch
            allProtocolsBLData = analysisDetails.allProtocolsBLData;  
            goodProtVals = getDetails(analysisDetails,'allProtocolsBLData',subjectName);
            goodProtFlag = goodProtVals{2};
            clear fileLists
            for iProt = 1:length(expDates)
                fileLists{iProt} = [subjectName '-' expDates{iProt} '-' protocolNames{iProt} '.mat'];
            end
            fileLists(~goodProtFlag)=[];
            
            numMSInRangePerProtocol = load(fullfile(saveFolder,'eyeData','microsaccades.mat'),'numMSInRangePerProtocol');
            numMSInRangePerProtocol = numMSInRangePerProtocol.numMSInRangePerProtocol;
            [~,stPowerVsFreqSF_Ori,blPowerVsFreqSF_Ori,freqValsSF_Ori,tfPowerSF_Ori,timeValsTF,freqValsTF,erpData,timeVals,numGoodTrialsSF_Ori,numAnalysedElecs]=...
                getDataSingleSubject(cleanDataFolder,fileLists,capLayout,electrodeList,stRange,1,numMSInRangePerProtocol); 
            
            rawAnalysis = getStructFromVars(stPowerVsFreqSF_Ori,blPowerVsFreqSF_Ori,freqValsSF_Ori,tfPowerSF_Ori,timeValsTF,freqValsTF,erpData,timeVals,numGoodTrialsSF_Ori,numAnalysedElecs); %#ok<NASGU>            
            makeDirectory(saveFolder);            
            save(fullfile(saveFolder,rawAnalysisFile),'-struct','rawAnalysis');        
        end    
        
        if any(goodProtFlag)
                
            % Get alpha and gamma Pos
            clear badFreqPos alphaPos
            badFreqPos = getBadFreqPos(freqValsSF_Ori); 
            alphaPos = setdiff(intersect(find(freqValsSF_Ori>=alphaRange(1)),find(freqValsSF_Ori<=alphaRange(2))),badFreqPos);
            gamma1Pos = setdiff(intersect(find(freqValsSF_Ori>=gamma1Range(1)),find(freqValsSF_Ori<=gamma1Range(2))),badFreqPos);
            gamma2Pos = setdiff(intersect(find(freqValsSF_Ori>=gamma2Range(1)),find(freqValsSF_Ori<=gamma2Range(2))),badFreqPos);

            % Calculate power area wise for SF Protocol
            diffPower = getDiffPowerSingleSubject(stPowerVsFreqSF_Ori,blPowerVsFreqSF_Ori,gamma1Pos,gamma2Pos,alphaPos);

            % Trim TF plot to 0-100 Hz to save space: MD 09-10-2018
            analysisDetails = updateDetails(analysisDetails,'gammaPowerValsAreaWiseNoMS',subjectName,{diffPower(:,g1Num) diffPower(:,g2Num) numGoodTrialsSF_Ori numAnalysedElecs}); 
            analysisDetails = updateDetails(analysisDetails,'alphaPowerValsAreaWiseNoMS',subjectName,{diffPower(:,alphaNum) numGoodTrialsSF_Ori numAnalysedElecs}); 
        else
            warning('No useful protocols for the subject...');
            analysisDetails = updateDetails(analysisDetails,'gammaPowerValsAreaWiseNoMS',subjectName,{[]}); 
        end                                            
        save(fullfile(workbookFolder,[projectName 'AnalysisDetails_' refType '.mat']),'-struct','analysisDetails','-append');
    end
catch err1
    disp(['............' err1.message]); 
    disp('gamma power without MS could not be calculated');
end
end
%}


end