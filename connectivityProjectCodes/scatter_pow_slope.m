pow_label = 'diff';
data_mci = load(['CaseVsControl_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{methodOptions.sideToShow} '_MCI']);
data_ad = load(['CaseVsControl_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{methodOptions.sideToShow} '_AD']);
data_healthy = load(['MidVsOld_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{methodOptions.sideToShow}]);
healthyColor = 'g'; mciColor = 'r'; adColor = 'y'; 
Bands = {'Alpha','Slow gamma','Fast gamma'};
ElecGLabels = {'Left', 'Right', 'Back'};

if(strcmp(pow_label,'diff'))
    t_xpos = [-2 1 1.5];
else % 'st' condition settings
    t_xpos = [10 -1 -5];
end

hf = figure; % slopes with band power
for freqBand = 1:3
    subplot(3,1,freqBand);
    title(Bands{freqBand});
    scatter(data_healthy.pow_elec{freqBand},data_healthy.slopes_elec{freqBand},30,'filled','MarkerFaceColor',healthyColor,'MarkerFaceAlpha',1); hold on;
    scatter(data_mci.pow_elec{freqBand},data_mci.slopes_elec{freqBand},30,'filled','MarkerFaceColor',mciColor,'MarkerFaceAlpha',0.75);
    scatter(data_ad.pow_elec{freqBand},data_ad.slopes_elec{freqBand},30,'filled','MarkerFaceColor',adColor,'MarkerFaceAlpha',1);

    if(strcmp(methodOptions.pow_label,'diff'))
        xlabel('Diff. power in dB');
    else
        xlabel('Absolute power in log10()');
    end
    ylabel('Mean connectivity');
    ylim([0 0.5]);
    set(gca,'FontSize',6);
end
    
print(hf,'-painters',fullfile(pwd,['SuppFig1_SlopeVsPow_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{methodOptions.sideToShow}]),'-dtiff','-r300');
    