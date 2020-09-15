function [eyeDataDegX,eyeDataDegY,eyeDataArbUnitsP,timeValsEyeData,trialNums,protNumForTrial,FsEyes] = getEyeDataIndividualSubject(fileLists,cleanDataFolder)
    eyeDataDegX = []; eyeDataDegY = []; eyeDataArbUnitsP = []; trialNums = []; protNumForTrial = [];
    for iProt = 1:length(fileLists);
        clear eyeData; eyeData = load(fullfile(cleanDataFolder,fileLists{iProt}),'timeValsEye','eyeData');
        timeValsEyeData = eyeData.timeValsEye;
        eyeData = eyeData.eyeData;

        eyeDataDegX = cat(1,eyeDataDegX,eyeData.eyeDataDegX);
        eyeDataDegY = cat(1,eyeDataDegY,eyeData.eyeDataDegY);
        eyeDataArbUnitsP = cat(1,eyeDataArbUnitsP,eyeData.eyeDataArbUnitsP);
        trialNums = cat(1,trialNums,(1:size(eyeData.eyeDataDegX,1))');
        protNumForTrial = cat(1,protNumForTrial,repmat(iProt,size(eyeData.eyeDataDegX,1),1));
        FsEyes{iProt} = single(1/(timeValsEyeData(2)-timeValsEyeData(1))); %#ok<AGROW> % FsEye is constant across protocols.
    end
end