function [MSAllData,MSAmplitude,peakVelocity,peakVelocityAll,MSLengthAll,MSAmplitudeAll,trialVals,nMS,MSAllDataExtended] = getMSAmplitudeAndPeakVelocity(MS,MSLength,eyeSpeedMag,eyeDataDegX,eyeDataDegY,timeRange,timeValsEyePos,FsEye)
    numTrials = size(eyeSpeedMag,1);
    peakVelocityAll=[]; MSLengthAll=[]; MSAmplitudeAll=[];
    preMSOnsetPeriod = 0.030; % 30 ms baseline and 70 ms post-MS onset 
    postMSOnsetPeriod = 0.070; % 30 ms baseline and 70 ms post-MS onset 
    
    dbstop if error
    MSAllData = cell(1,numTrials);
    MSAllDataExtendedX = [];
    MSAllDataExtendedY = [];
    MSAmplitude = cell(1,numTrials);
    peakVelocity = cell(1,numTrials);
    trialVals = NaN(size(eyeSpeedMag));
    nMS = 0;
    for iTrial = 1:numTrials
        clear trialValsLogical MSForTrial MSToDiscard MSLengthForTrial
        trialValsLogical = false(1,size(eyeSpeedMag,2));
        MSForTrial = MS{iTrial};
        MSLengthForTrial = MSLength{iTrial};

        MSToDiscard = MSForTrial<timeRange(1)|MSForTrial>timeRange(2);
        MSForTrial(MSToDiscard)=[];
        MSLengthForTrial(MSToDiscard)=[];
        for iMS = 1:length(MSForTrial)
            clear MSTimeVals valsForTrial
            MSTimeVals = MSForTrial(iMS):1/FsEye:MSForTrial(iMS)+MSLengthForTrial(iMS)-1/FsEye;   
            valsForMS = ismember(int16(timeValsEyePos*FsEye),int16(MSTimeVals*FsEye));

            valsForMSIndices = find(valsForMS);
            xMS = eyeDataDegX(iTrial,valsForMSIndices(1):valsForMSIndices(end)+1);
            MSAllData{iTrial}{iMS}(1,:) = xMS-repmat(xMS(1),1,length(xMS));
            
            yMS = eyeDataDegY(iTrial,valsForMSIndices(1):valsForMSIndices(end)+1);
            MSAllData{iTrial}{iMS}(2,:) = yMS-repmat(yMS(1),1,length(yMS));
            
            if (MSForTrial(iMS)>(preMSOnsetPeriod+1/FsEye)) && (MSForTrial(iMS)<(timeValsEyePos(end)-postMSOnsetPeriod-1/FsEye))
                extendedMSTimeVals = MSForTrial(iMS)-preMSOnsetPeriod:1/FsEye:MSForTrial(iMS)+postMSOnsetPeriod-1/FsEye;
                extendedValsForTrial = ismember(int16(timeValsEyePos*FsEye),int16(extendedMSTimeVals*FsEye));            
                MSAllDataExtendedX = cat(1,MSAllDataExtendedX,eyeDataDegX(iTrial,extendedValsForTrial));
                MSAllDataExtendedY = cat(1,MSAllDataExtendedY,eyeDataDegY(iTrial,extendedValsForTrial));
            end
            
%             MSMagTotal = sqrt(diff([xMS(1) xMS(end)]).^2 + diff([yMS(1) yMS(end)]).^2);
%             MSMagPerMovement = max(sqrt(diff(xMS).^2 + diff(yMS).^2));
%             MSAmplitude{iTrial}(iMS) = max([MSMagTotal,MSMagPerMovement]);
            peakVelocity{iTrial}(iMS) = max(eyeSpeedMag(iTrial,valsForMS));
            
            clear Combi MSMagPerMovement
            Combi = nchoosek(1:length(xMS),2);
            for iMov = 1:size(Combi,1)
                MSMagPerMovement(iMov) = (sqrt(diff(xMS(Combi(iMov,:))).^2 + diff(yMS(Combi(iMov,:))).^2));
            end
            MSAmplitude{iTrial}(iMS) = max(MSMagPerMovement);

            peakVelocityAll = cat(2,peakVelocityAll,peakVelocity{iTrial}(iMS));
            MSLengthAll = cat(2,MSLengthAll,MSLengthForTrial(iMS));
            MSAmplitudeAll = cat(2,MSAmplitudeAll,MSAmplitude{iTrial}(iMS));

            trialValsLogical = trialValsLogical|valsForMS;
            nMS = nMS+1;

        end
        trialVals(iTrial,:) = double(trialValsLogical);
    end
    
    MSAllDataExtended(1,:,:) = MSAllDataExtendedX;
    MSAllDataExtended(2,:,:) = MSAllDataExtendedY;
end