pow_label = 'diff';
MiddleAgedColor = 'b'; ElderlyColor = 'c'; mciColor = 'r'; adColor = 'y'; 
Bands = {'Alpha','Slow gamma','Fast gamma'};
ElecGLabels = {'Left', 'Right', 'Back'};

data_mci = load(['CaseVsControl_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{methodOptions.sideToShow} '_MCI']);
% data_ad = load(['CaseVsControl_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{methodOptions.sideToShow} '_AD']);
data_healthy = load(['MidVsOld_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{methodOptions.sideToShow}]);

if(strcmp(pow_label,'diff'))
    t_xpos = [-2 1 1.5];
else % 'st' condition settings
    t_xpos = [10 -1 -5];
end

%%
figR = figure('numbertitle','off','name','Scatter-plot');
figR.PaperType = 'a4';
figR.PaperUnits = 'centimeters';
figR.PaperSize = [32 20];
figR.PaperOrientation = 'Landscape';
figR.PaperPosition = [0 0 figR.PaperSize];
figR.Color = [1 1 1]; % White background

for freqBand = 1:3
    subplot(3,1,freqBand);
    title(Bands{freqBand});
    middleAgedSubjs = find(data_healthy.subj_label{freqBand}==0); elderlySubjs = find(data_healthy.subj_label{freqBand}==1);
    scatter(data_healthy.pow_elec{freqBand}(middleAgedSubjs),data_healthy.slopes_elec{freqBand}(middleAgedSubjs),30,'filled','MarkerFaceColor',MiddleAgedColor,'MarkerFaceAlpha',1); hold on;
    scatter(data_healthy.pow_elec{freqBand}(elderlySubjs),data_healthy.slopes_elec{freqBand}(elderlySubjs),30,'filled','MarkerFaceColor',ElderlyColor,'MarkerFaceAlpha',0.75);
%     scatter(data_mci.pow_elec{freqBand},data_mci.slopes_elec{freqBand}{1},30,'filled','MarkerFaceColor',mciColor,'MarkerFaceAlpha',0.75);
%     scatter(data_ad.pow_elec{freqBand},data_ad.slopes_elec{freqBand},30,'filled','MarkerFaceColor',adColor,'MarkerFaceAlpha',1);

    if(strcmp(methodOptions.pow_label,'diff'))
        xlabel('Diff. power in dB');
    else
        xlabel('Absolute power in log10()');
    end
    ylabel([Bands{freqBand} '-Conn']);
    ylim([0 0.5]);
    if(freqBand ~=1)
        xlim([-2 4]);
    else
        xlim([-4 2]);
    end
    line([0 0],[0 0.5],'LineWidth',2,'Color','k');
    set(gca,'FontSize',10);
    if(freqBand == 3)
        legend('MiddleAged','Elderly','MCI');
    end
end
print(figR,'-painters',fullfile(pwd,['SuppFig_SlopeVsPow_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{methodOptions.sideToShow}]),'-dtiff','-r300');

%% regression1: need conn, age, change in pow
ageListRaw = data_healthy.methodOptions.MidVsOld.ageList';
ageGroup1Pos = (ageListRaw<65); ageGroup2Pos = (ageListRaw>=65);
age1 = ageListRaw(ageGroup1Pos);
age2 = ageListRaw(ageGroup2Pos);
ageList = [age1; age2];

t_xpos = [-3.5 -1.5 -1.5]; p_xpos = [0.6 2.5 2.5];

for freqBand = 1:3
    hold on;
    subplot(3,1,freqBand);
    mdl = fitlm([ageList data_healthy.pow_elec{freqBand}],data_healthy.slopes_elec{freqBand});
    bCoeffs = table2array(mdl.Coefficients(:,1));
    text(t_xpos(freqBand),0.45,['slope = ' num2str(bCoeffs(1)) '  + ' num2str(bCoeffs(2)) ' *age  + ' num2str(bCoeffs(3)) ' *pow'],'FontSize',6,'BackgroundColor','w');
    text(t_xpos(freqBand),0.4,['pval = (' num2str(table2array(mdl.Coefficients(2,4)),'%.2e') ',' num2str(table2array(mdl.Coefficients(3,4)),'%.3e') ')'],'FontSize',8);

    if(freqBand==1)
        selSubj = data_healthy.pow_elec{freqBand}<0;
    else
        selSubj = data_healthy.pow_elec{freqBand}>0;
    end
    
    mdl = fitlm([ageList(selSubj) data_healthy.pow_elec{freqBand}(selSubj)],data_healthy.slopes_elec{freqBand}(selSubj));
    bCoeffs = table2array(mdl.Coefficients(:,1));
    text(p_xpos(freqBand),0.55,['slope = ' num2str(bCoeffs(1)) '  + ' num2str(bCoeffs(2)) ' *age  + ' num2str(bCoeffs(3)) ' *pow'],'FontSize',6,'BackgroundColor','w');
    text(p_xpos(freqBand),0.5,['pval = (' num2str(table2array(mdl.Coefficients(2,4)),'%.2e') ',' num2str(table2array(mdl.Coefficients(3,4)),'%.3e') ')'],'FontSize',8);
    
end

print(figR,'-painters',fullfile(pwd,['SuppFig_SlopeVsPow_withRegress_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{methodOptions.sideToShow}]),'-dtiff','-r300');
