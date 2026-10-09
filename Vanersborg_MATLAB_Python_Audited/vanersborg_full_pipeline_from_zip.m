function vanersborg_full_pipeline_from_zip(zip1,zip2,outDir)
%VANERSBORG_FULL_PIPELINE_FROM_ZIP Complete ZIP-to-results Vänersborg workflow.
%
% Starts from the two original public raw data archives:
%   Dataset_steel_bridge_in_Sweden_1.zip
%   Dataset_steel_bridge_in_Sweden_2.zip
%
% The script performs:
%   1) raw ZIP extraction and CSV inventory;
%   2) cleaning and bridge-angle alignment of all 64 opening records;
%   3) assembly of 12 usable strain channels on a 256-point opening/closing grid;
%   4) physical-sensor benchmark with Area + RMS/P2P + Gram + DCT8 features;
%   5) healthy POD + pivoted-QR physical reference selection;
%   6) healthy-only ridge reconstruction of the seven withheld strain channels;
%   7) blind virtual-sensor validation on independent baseline records;
%   8) matched comparison of all physical, reduced QR physical, and QR+virtual;
%   9) event results, repeated-window resampling and all 20 baseline 3/3 splits;
%  10) CSV tables, MAT cache, PNG/PDF figures and summary text.
%
% IMPORTANT INTERPRETATION
% The virtual channels are data-driven reconstructions learned from historical
% healthy full-array measurements. They are NOT FE-generated virtual sensors
% and are NOT additional independent measurements.
%
% PRIMARY PROTOCOL
% Usable gauges: SG1, SG2, SG3, SG7, SG8, SG9, SG10, SG11, SG13, SG14, SG15, SG16
% Raw angle:    ch_25 (sign reversed)
% Temperature:  ch_30
% Fit records:  first 3 early opening events
% Calibration:  next 3 early opening events
% Crack file:   2023_03_09T23_45_25.csv
% Damage test:  crack record + six later records
% QR budget:    smallest healthy POD rank retaining >=99.9%% energy
% Ridge alpha:  0.01
% Truth zone:   SG9--SG10
%
% EXAMPLE
%   vanersborg_full_pipeline_from_zip( ...
%       'Dataset_steel_bridge_in_Sweden_1.zip', ...
%       'Dataset_steel_bridge_in_Sweden_2.zip', ...
%       'vanersborg_full_results_matlab');
%
% MATLAB R2021b or later is recommended. The implementation uses an internal
% DCT-II matrix and does not require Signal Processing Toolbox.

if nargin < 1 || isempty(zip1)
    zip1 = fullfile(fileparts(mfilename('fullpath')),'Dataset_steel_bridge_in_Sweden_1.zip');
end
if nargin < 2 || isempty(zip2)
    zip2 = fullfile(fileparts(mfilename('fullpath')),'Dataset_steel_bridge_in_Sweden_2.zip');
end
if nargin < 3 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')),'vanersborg_full_results_matlab');
end
zip1=char(zip1); zip2=char(zip2); outDir=char(outDir);
assert(isfile(zip1),'Missing raw archive: %s',zip1);
assert(isfile(zip2),'Missing raw archive: %s',zip2);
if ~exist(outDir,'dir'); mkdir(outDir); end

% Input manifest.
d1=dir(zip1); d2=dir(zip2);
InputManifest=table(string({d1.name;d2.name}),[d1.bytes;d2.bytes],string({zip1;zip2}), ...
    'VariableNames',{'file','size_bytes','path'});
writetable(InputManifest,fullfile(outDir,'00_input_manifest.csv'));

% Raw extraction/alignment.
[X,files,gauges,temp] = buildAlignedDataset(zip1,zip2,outDir);
assert(isequal(size(X),[64 256 12]),'Unexpected aligned array shape.');

% Plot a representative aligned full-array baseline record.
f0=figure('Color','w','Position',[100 100 950 540]); hold on;
for j=1:numel(gauges)
    plot(1:size(X,2),squeeze(X(1,:,j)),'LineWidth',0.8,'DisplayName',char("SG"+string(gauges(j))));
end
grid on; xlabel('Aligned opening/closing sample'); ylabel('Strain');
title('Vänersborg aligned baseline strain field: first reference record');
legend('Location','eastoutside');
exportgraphics(f0,fullfile(outDir,'00_aligned_baseline_strain_field.png'),'Resolution',300);
exportgraphics(f0,fullfile(outDir,'00_aligned_baseline_strain_field.pdf'),'ContentType','vector'); close(f0);

crackFile = "2023_03_09T23_45_25.csv";
crackIdx = find(files==crackFile,1);
assert(~isempty(crackIdx),'Documented crack file not found.');
damageIdx = crackIdx:numel(files);
fitIdx = 1:3;
calIdx = 4:6;
truthZone = "Z05_SG9_SG10";

% -------------------------------------------------------------------------
% QR reference selection from healthy baseline snapshots
% -------------------------------------------------------------------------
[refIdx,cumEnergy] = selectQRReferences(X,gauges,fitIdx,0.999);
assert(~isempty(refIdx) && all(refIdx>=1 & refIdx<=numel(gauges)), ...
    'QR selection returned no valid physical gauge indices.');
refGauges = gauges(refIdx);
fprintf('QR references (%d): %s\n',numel(refGauges),char(join("SG"+string(refGauges),", ")));
withheldIdx = setdiff(1:numel(gauges),refIdx,'stable');
withheldGauges = gauges(withheldIdx);

% -------------------------------------------------------------------------
% Frozen virtual strain transfer
% -------------------------------------------------------------------------
[Xvirtual,VirtualModel] = fitVirtualTransfer(X,refIdx,withheldIdx,fitIdx,0.01);

% Blind validation on the three calibration records
valEvent = strings(0,1); valGauge = strings(0,1);
valNRMSE = zeros(0,1); valCorr = zeros(0,1); valAmp = zeros(0,1);
for e = calIdx
    for j = withheldIdx
        measured = squeeze(X(e,:,j)).';
        predicted = squeeze(Xvirtual(e,:,j)).';
        rmse = sqrt(mean((predicted-measured).^2));
        valEvent(end+1,1) = files(e); %#ok<AGROW>
        valGauge(end+1,1) = "SG"+string(gauges(j)); %#ok<AGROW>
        valNRMSE(end+1,1) = rmse/max(range1(measured),1e-12); %#ok<AGROW>
        C = corrcoef(measured,predicted);
        if numel(C)>=4, cc=C(1,2); else, cc=NaN; end
        valCorr(end+1,1) = cc; %#ok<AGROW>
        valAmp(end+1,1) = range1(predicted)/max(range1(measured),1e-12); %#ok<AGROW>
    end
end
VirtualValidation = table(valEvent,valGauge,valNRMSE,valCorr,valAmp, ...
    'VariableNames',{'event','gauge','NRMSE_peak_to_peak','waveform_correlation','amplitude_ratio'});
writetable(VirtualValidation,fullfile(outDir,'03_virtual_sensor_blind_validation.csv'));

% Save the frozen healthy virtual-transfer model.
save(fullfile(outDir,'08_frozen_virtual_transfer_model.mat'),'VirtualModel','refIdx','refGauges','withheldIdx','withheldGauges','-v7');

% Representative blind validation waveform (first calibration event / first withheld gauge).
e0=calIdx(1); j0=withheldIdx(1);
fvs=figure('Color','w','Position',[100 100 900 500]); hold on;
plot(squeeze(X(e0,:,j0)),'LineWidth',1.2,'DisplayName',char("Measured SG"+string(gauges(j0))));
plot(squeeze(Xvirtual(e0,:,j0)),'--','LineWidth',1.2,'DisplayName',char("Virtual SG"+string(gauges(j0))));
grid on; xlabel('Aligned opening/closing sample'); ylabel('Strain'); title('Representative blind virtual-strain validation'); legend('Location','best');
exportgraphics(fvs,fullfile(outDir,'05_virtual_sensor_waveform_example.png'),'Resolution',300);
exportgraphics(fvs,fullfile(outDir,'05_virtual_sensor_waveform_example.pdf'),'ContentType','vector'); close(fvs);

% -------------------------------------------------------------------------
% Three matched configurations
% -------------------------------------------------------------------------
C1 = evaluateConfiguration("Physical sensors only",X,gauges,fitIdx,calIdx,damageIdx,crackIdx,truthZone,"12 measured gauges");
C1.physicalSensorCount = 12;
Xred = X(:,:,refIdx);
C2 = evaluateConfiguration("Reduced QR physical sensors",Xred,refGauges,fitIdx,calIdx,damageIdx,crackIdx,truthZone,join("SG"+string(refGauges),", "));
C2.physicalSensorCount = numel(refGauges);
C3 = evaluateConfiguration("QR references + virtual sensors",Xvirtual,gauges,fitIdx,calIdx,damageIdx,crackIdx,truthZone, ...
    string(numel(refGauges))+" measured refs + "+string(numel(withheldGauges))+" virtual gauges");
C3.physicalSensorCount = numel(refGauges);

Summary = table( ...
    [C1.label;C2.label;C3.label], ...
    [C1.physicalSensorCount;C2.physicalSensorCount;C3.physicalSensorCount], ...
    [C1.representation;C2.representation;C3.representation], ...
    [C1.threshold;C2.threshold;C3.threshold], ...
    [C1.detectedOf7;C2.detectedOf7;C3.detectedOf7], ...
    [C1.correctOf7;C2.correctOf7;C3.correctOf7], ...
    [C1.jointOf7;C2.jointOf7;C3.jointOf7], ...
    [C1.medianRatio;C2.medianRatio;C3.medianRatio], ...
    [C1.crackDetected;C2.crackDetected;C3.crackDetected], ...
    [C1.crackPredicted;C2.crackPredicted;C3.crackPredicted], ...
    [C1.preDocumentationAlerts;C2.preDocumentationAlerts;C3.preDocumentationAlerts], ...
    'VariableNames',{'configuration','physical_sensor_count','representation','threshold', ...
    'detected_of_7','SG9_SG10_top_zone_of_7','joint_detect_and_localise_of_7', ...
    'median_damage_score_over_threshold','documented_crack_detected', ...
    'documented_crack_predicted_zone','pre_documentation_alerts_of_51'});
writetable(Summary,fullfile(outDir,'01_three_configuration_summary.csv'));

% Crack/post event results
EventRows = table();
configs = {C1,C2,C3};
for cc = 1:3
    C = configs{cc};
    for e = damageIdx
        row = table(C.label,files(e),e==crackIdx,C.maxScore(e),C.threshold, ...
            C.maxScore(e)/C.threshold,C.detected(e),C.zoneNames(C.pred(e)), ...
            C.correct(e),C.detected(e)&&C.correct(e), ...
            'VariableNames',{'configuration','file','is_documented_crack_record','score','threshold', ...
            'score_over_threshold','detected','predicted_zone','correct_SG9_SG10_zone','joint_correct'});
        EventRows = [EventRows;row]; %#ok<AGROW>
    end
end
writetable(EventRows,fullfile(outDir,'02_crack_and_post_event_results.csv'));

% QR table / POD energy
QRTable = table((1:numel(refGauges)).',("SG"+string(refGauges)).', ...
    repmat("smallest POD rank reaching >=99.9% fit-record spatial energy",numel(refGauges),1), ...
    repmat(cumEnergy(numel(refGauges)),numel(refGauges),1), ...
    'VariableNames',{'QR_rank','physical_reference','selection_rule','cumulative_energy_at_selected_rank'});
writetable(QRTable,fullfile(outDir,'06_QR_reference_selection.csv'));
POD = table((1:numel(cumEnergy)).',cumEnergy(:),'VariableNames',{'POD_rank','cumulative_healthy_spatial_energy'});
writetable(POD,fullfile(outDir,'07_healthy_spatial_POD_energy.csv'));

% -------------------------------------------------------------------------
% Conditional resampling of the seven crack/post records
% -------------------------------------------------------------------------
rng(20261005,'twister');
Repeatability = table();
for cc=1:3
    C=configs{cc};
    for window=[1 3 5]
        calDraw = vanersborg_core('draws',calIdx,5000,window);
        calAccum = zeros(5000,size(C.scores,2));
        for k=1:window, calAccum=calAccum+C.scores(calDraw(:,k),:); end
        calMax=max(calAccum,[],2);
        thresholdAccum=empiricalQuantileHigher(calMax,0.95);

        dmgDraw = vanersborg_core('draws',damageIdx,1000,window);
        dmgAccum=zeros(1000,size(C.scores,2));
        for k=1:window, dmgAccum=dmgAccum+C.scores(dmgDraw(:,k),:); end
        [mx,pp]=max(dmgAccum,[],2);
        det=mx>thresholdAccum;
        correct=false(1000,1);
        for r=1:1000, correct(r)=C.zoneNames(pp(r))==truthZone; end
        row=table(C.label,window,1000,thresholdAccum,mean(det),mean(correct),mean(det&correct), ...
            'VariableNames',{'configuration','window_events','trials','accumulated_threshold', ...
            'detection_rate','SG9_SG10_selection_rate','joint_rate'});
        Repeatability=[Repeatability;row]; %#ok<AGROW>
    end
end
writetable(Repeatability,fullfile(outDir,'04_repeatability_bootstrap.csv'));

% -------------------------------------------------------------------------
% All 20 disjoint baseline 3/3 split sensitivity checks
% -------------------------------------------------------------------------
comb = nchoosek(1:6,3);
Split = table();
for rr=1:size(comb,1)
    tr=comb(rr,:); ca=setdiff(1:6,tr,'stable');
    A = evaluateSplit(X,gauges,tr,ca,damageIdx,truthZone);
    Split=[Split;splitRow("Physical sensors only",tr,"",A)]; %#ok<AGROW>

    % Matched five-reference QR budget for split sensitivity.
    ri=selectQRFixedCount(X,tr,5);
    rg=gauges(ri);
    A=evaluateSplit(X(:,:,ri),rg,tr,ca,damageIdx,truthZone);
    refText=join("SG"+string(rg),",");
    Split=[Split;splitRow("Reduced QR physical sensors",tr,refText,A)]; %#ok<AGROW>

    wi=setdiff(1:numel(gauges),ri,'stable');
    xv=fitVirtualTransfer(X,ri,wi,tr,0.01);
    A=evaluateSplit(xv,gauges,tr,ca,damageIdx,truthZone);
    Split=[Split;splitRow("QR references + virtual sensors",tr,refText,A)]; %#ok<AGROW>
end
writetable(Split,fullfile(outDir,'05_baseline_split_sensitivity.csv'));

% -------------------------------------------------------------------------
% Figures
% -------------------------------------------------------------------------
labels = categorical(["All physical","QR physical","QR + virtual"]);
labels = reordercats(labels,["All physical","QR physical","QR + virtual"]);
vals=[100*Summary.detected_of_7/7,100*Summary.SG9_SG10_top_zone_of_7/7,100*Summary.joint_detect_and_localise_of_7/7];
f=figure('Color','w','Position',[100 100 900 520]);
bar(labels,vals); ylim([0 105]); grid on;
ylabel('Fraction of seven crack/post-crack records (%)');
title('Vänersborg experimental localisation comparison');
legend({'Detected','SG9--SG10 top zone','Detected + correct'},'Location','southoutside','Orientation','horizontal');
exportgraphics(f,fullfile(outDir,'01_damage_localisation_comparison.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'01_damage_localisation_comparison.pdf'),'ContentType','vector'); close(f);

ug=unique(VirtualValidation.gauge,'stable'); med=zeros(numel(ug),1);
for i=1:numel(ug), med(i)=median(VirtualValidation.NRMSE_peak_to_peak(VirtualValidation.gauge==ug(i))); end
f=figure('Color','w','Position',[100 100 900 520]); bar(categorical(ug),100*med); grid on;
ylabel('Median NRMSE (% of measured peak-to-peak)'); title('Blind validation of virtual strain gauges');
exportgraphics(f,fullfile(outDir,'02_virtual_sensor_validation_NRMSE.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'02_virtual_sensor_validation_NRMSE.pdf'),'ContentType','vector'); close(f);

f=figure('Color','w','Position',[100 100 950 540]); hold on;
for cc=1:3
    C=configs{cc}; plot(1:numel(damageIdx),C.maxScore(damageIdx)/C.threshold,'-o','LineWidth',1.2,'DisplayName',C.label);
end
yline(1,'--'); grid on; xticks(1:7); xticklabels({'Crack','Post 1','Post 2','Post 3','Post 4','Post 5','Post 6'});
ylabel('Damage score / threshold'); title('Eventwise detection margin'); legend('Location','best');
exportgraphics(f,fullfile(outDir,'03_eventwise_detection_margin.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'03_eventwise_detection_margin.pdf'),'ContentType','vector'); close(f);

% Baseline split sensitivity: min/median/max joint counts.
cn=["Physical sensors only","QR references + virtual sensors","Reduced QR physical sensors"];
medv=zeros(3,1); minv=medv; maxv=medv;
for i=1:3
    q=Split.joint_of_7(Split.configuration==cn(i)); minv(i)=min(q); medv(i)=median(q); maxv(i)=max(q);
end
f=figure('Color','w','Position',[100 100 900 520]); errorbar(1:3,100*medv/7,100*(medv-minv)/7,100*(maxv-medv)/7,'o','LineWidth',1.2); grid on;
xlim([0.5 3.5]); ylim([-3 105]); xticks(1:3); xticklabels({'All physical','QR + virtual','QR physical'});
ylabel('Detected + correctly localised (%)'); title('Sensitivity to all 20 early-baseline 3/3 splits');
exportgraphics(f,fullfile(outDir,'04_baseline_split_sensitivity.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'04_baseline_split_sensitivity.pdf'),'ContentType','vector'); close(f);

% -------------------------------------------------------------------------
% Conceptual / manuscript-oriented result figures
% -------------------------------------------------------------------------
% 06 information-recovery pathway.
f=figure('Color','w','Position',[100 100 1150 520]); ax=axes(f); hold(ax,'on');
axis(ax,[0 12 0 5.2]); axis(ax,'off');
boxData={...
    0.4,1.55,2.2,2.0,"12 physical\nstrain gauges",sprintf('%d/%d joint\nlocalisation',C1.jointOf7,numel(damageIdx));...
    3.3,1.55,2.2,2.0,sprintf('QR selection\n%d references',numel(refIdx)),sprintf('%d/%d joint\nafter reduction',C2.jointOf7,numel(damageIdx));...
    6.2,1.55,2.2,2.0,"Healthy-trained\nvirtual transfer",sprintf('%d withheld channels\nreconstructed',numel(withheldIdx));...
    9.1,1.55,2.2,2.0,sprintf('%d measured +\n%d virtual',numel(refIdx),numel(withheldIdx)),sprintf('%d/%d joint\nlocalisation',C3.jointOf7,numel(damageIdx))};
for ii=1:size(boxData,1)
    x0=boxData{ii,1}; y0=boxData{ii,2}; ww=boxData{ii,3}; hh=boxData{ii,4};
    rectangle(ax,'Position',[x0 y0 ww hh],'Curvature',0.08,'LineWidth',1.1);
    text(ax,x0+ww/2,y0+1.25,sprintf(char(boxData{ii,5})),'HorizontalAlignment','center','FontWeight','bold','FontSize',12);
    text(ax,x0+ww/2,y0+0.52,sprintf(char(boxData{ii,6})),'HorizontalAlignment','center','FontSize',10);
end
for pair=[2.6 3.3;5.5 6.2;8.4 9.1].'
    x1=pair(1); x2=pair(2);
    quiver(ax,x1,2.55,x2-x1-0.05,0,0,'MaxHeadSize',0.35,'LineWidth',1.2);
end
text(ax,6.0,0.45,sprintf(['Virtual channels are historically calibrated reconstructions of removed physical gauges;\n' ...
    'they are not independent measurements and are not FE-generated sensors.']), ...
    'HorizontalAlignment','center','FontSize',10);
title(ax,'Vänersborg sensing experiment: information loss and recovery');
exportgraphics(f,fullfile(outDir,'06_information_recovery_pathway.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'06_information_recovery_pathway.pdf'),'ContentType','vector'); close(f);

% 07 documented-crack spatial fingerprint.
allZoneNames=["Z01_SG1_SG2","Z02_SG3","Z03_SG7","Z04_SG8", ...
              "Z05_SG9_SG10","Z06_SG11","Z07_SG13_SG14","Z08_SG15_SG16"];
M=nan(3,numel(allZoneNames));
for cc=1:3
    C=configs{cc};
    for jj=1:numel(C.zoneNames)
        col=find(allZoneNames==C.zoneNames(jj),1);
        M(cc,col)=C.scores(crackIdx,jj);
    end
end
Mnorm=M./max(M,[],2,'omitnan');
f=figure('Color','w','Position',[100 100 1100 480]);
imagesc(Mnorm); colorbar; caxis([0 1]);
set(gca,'YTick',1:3,'YTickLabel',{'All physical (12)','QR physical (5)','QR + virtual (5+7)'}, ...
    'XTick',1:8,'XTickLabel',{'SG1-2','SG3','SG7','SG8','SG9-10','SG11','SG13-14','SG15-16'});
xlabel('Sensor zone'); title('Spatial anomaly fingerprint at the documented crack record');
for ii=1:3
    for jj=1:8
        if isfinite(Mnorm(ii,jj)), txt=sprintf('%.2f',Mnorm(ii,jj)); else, txt='N/A'; end
        text(jj,ii,txt,'HorizontalAlignment','center','FontSize',9);
    end
end
rectangle('Position',[4.5 0.5 1 3],'EdgeColor','k','LineWidth',2.2);
text(5,0.25,'inspection truth','HorizontalAlignment','center','FontSize',9);
exportgraphics(f,fullfile(outDir,'07_crack_zone_fingerprint.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'07_crack_zone_fingerprint.pdf'),'ContentType','vector'); close(f);

% Crack-zone table used by the fingerprint plot.
CZ=table();
for cc=1:3
    C=configs{cc};
    for jj=1:8
        raw=M(cc,jj); normv=Mnorm(cc,jj);
        row=table(C.label,allZoneNames(jj),raw,normv,allZoneNames(jj)==truthZone, ...
            'VariableNames',{'configuration','zone','raw_score','normalized_score','inspection_truth'});
        CZ=[CZ;row]; %#ok<AGROW>
    end
end
writetable(CZ,fullfile(outDir,'10_documented_crack_zone_scores.csv'));

% 08 virtual-sensor validation map.
ug=unique(VirtualValidation.gauge,'stable'); vm=zeros(numel(ug),3);
for ii=1:numel(ug)
    q=VirtualValidation.gauge==ug(ii);
    vm(ii,1)=median(VirtualValidation.NRMSE_peak_to_peak(q),'omitnan');
    vm(ii,2)=median(VirtualValidation.waveform_correlation(q),'omitnan');
    vm(ii,3)=median(VirtualValidation.amplitude_ratio(q),'omitnan');
end
f=figure('Color','w','Position',[100 100 760 560]); scatter(100*vm(:,1),vm(:,2),65,'filled'); hold on;
for ii=1:numel(ug), text(100*vm(ii,1),vm(ii,2),"  "+ug(ii),'FontSize',9); end
grid on; xlabel('Median NRMSE (% of measured peak-to-peak strain)'); ylabel('Median waveform correlation');
title('Blind validation of reconstructed strain channels');
exportgraphics(f,fullfile(outDir,'08_virtual_sensor_validation_map.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'08_virtual_sensor_validation_map.pdf'),'ContentType','vector'); close(f);

% 09 all 20 baseline partitions rather than only min/median/max.
f=figure('Color','w','Position',[100 100 920 560]); hold on;
plotLabels=["Physical sensors only","Reduced QR physical sensors","QR references + virtual sensors"];
for ii=1:3
    vals=Split.joint_of_7(Split.configuration==plotLabels(ii));
    offs=linspace(-0.17,0.17,numel(vals)).';
    scatter(ii+offs,vals,34,'filled');
    medx=median(vals); plot([ii-0.22 ii+0.22],[medx medx],'LineWidth',2);
end
grid on; xlim([0.5 3.5]); ylim([-0.3 7.5]); yticks(0:7);
xticks(1:3); xticklabels({'All physical (12)','QR physical (5)','QR + virtual (5+7)'});
ylabel('Detected + correctly localised records (out of 7)');
title('Robustness across all 20 early-baseline 3/3 partitions');
exportgraphics(f,fullfile(outDir,'09_baseline_partition_robustness.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'09_baseline_partition_robustness.pdf'),'ContentType','vector'); close(f);

FigureManifest=table(["01_damage_localisation_comparison";"06_information_recovery_pathway"; ...
    "07_crack_zone_fingerprint";"08_virtual_sensor_validation_map";"09_baseline_partition_robustness"], ...
    ["Primary quantitative result";"Conceptual contribution";"Spatial localisation evidence"; ...
     "Virtual-sensor credibility";"Robustness / limitation"], ...
    'VariableNames',{'figure','role'});
writetable(FigureManifest,fullfile(outDir,'11_figure_manifest.csv'));

fprintf('\nSelected QR references: %s\n',char(join("SG"+string(refGauges),", ")));
fprintf('Withheld virtual gauges: %s\n',char(join("SG"+string(withheldGauges),", ")));
fprintf('99.9%% spatial-energy rank: %d (%.4f%%)\n',numel(refGauges),100*cumEnergy(numel(refGauges)));
disp(Summary(:,{'configuration','detected_of_7','SG9_SG10_top_zone_of_7','joint_detect_and_localise_of_7'}));
fprintf('Median virtual NRMSE: %.2f%%\n',100*median(VirtualValidation.NRMSE_peak_to_peak));
fprintf('Median virtual correlation: %.4f\n',median(VirtualValidation.waveform_correlation));
end

% =========================================================================
% Local functions
% =========================================================================
function [refIdx,cumEnergy]=selectQRReferences(X,gauges,fitIdx,targetEnergy)
[refIdx,cumEnergy]=vanersborg_core('select',X,gauges,fitIdx,targetEnergy);
end
function refIdx=selectQRFixedCount(X,fitIdx,k)
refIdx=vanersborg_core('select_fixed',X,fitIdx,k);
end
function [Xv,model]=fitVirtualTransfer(X,refIdx,withheldIdx,fitIdx,alpha)
[Xv,model]=vanersborg_core('virtual',X,refIdx,withheldIdx,fitIdx,alpha);
end

function C=evaluateConfiguration(label,Xin,availableGauges,fitIdx,calIdx,damageIdx,crackIdx,truthZone,representation)
C=vanersborg_core('evaluate',Xin,availableGauges,fitIdx,calIdx,damageIdx,crackIdx,char(truthZone));
C.label=string(label); C.representation=string(representation);
C.zoneNames=string(C.zoneNames); C.crackPredicted=string(C.crackPredicted);
end

function A=evaluateSplit(Xin,availableGauges,tr,ca,damageIdx,truthZone)
C=vanersborg_core('evaluate',Xin,availableGauges,tr,ca,damageIdx,damageIdx(1),char(truthZone));
A.detected=C.detectedOf7; A.correct=C.correctOf7; A.joint=C.jointOf7;
end

function row=splitRow(label,tr,references,A)
row=table(string(label),join(string(tr-1),","),string(references),A.detected,A.correct,A.joint, ...
    'VariableNames',{'configuration','fit_indices','references','detected_of_7','correct_zone_of_7','joint_of_7'});
end

function q=empiricalQuantileHigher(x,p)
x=sort(x(:)); idx=max(1,min(numel(x),ceil(p*numel(x)))); q=x(idx);
end

function r=range1(x)
r=max(x)-min(x);
end

% =========================================================================
% Raw data extraction and alignment
% =========================================================================
function [X,files,gauges,temp]=buildAlignedDataset(zip1,zip2,outDir)
gauges=[1 2 3 7 8 9 10 11 13 14 15 16];
strainCols="ch_"+string(gauges);
angleCol="ch_25"; tempCol="ch_30";

extractRoot=fullfile(outDir,'00_raw_extracted');
if isfolder(extractRoot), rmdir(extractRoot,'s'); end
mkdir(extractRoot);
unzip(zip1,extractRoot); unzip(zip2,extractRoot);
D=dir(fullfile(extractRoot,'**','*.csv'));
assert(numel(D)==64,'Expected 64 raw opening CSV files; found %d.',numel(D));
[~,ord]=sort(string({D.name})); D=D(ord);
files=string({D.name}).';
assert(any(files=="2023_03_09T23_45_25.csv"),'Documented crack file missing.');

X=zeros(64,256,numel(gauges)); temp=zeros(64,1);
rawSamples=zeros(64,1); i2v=zeros(64,1); i72v=zeros(64,1); j72v=zeros(64,1); j2v=zeros(64,1);
b0v=zeros(64,1); b1v=zeros(64,1); archive=strings(64,1); fullPath=strings(64,1);

for e=1:64
    path=fullfile(D(e).folder,D(e).name); fullPath(e)=string(path);
    if contains(path,'Dataset_steel_bridge_in_Sweden_1'), archive(e)="Dataset_steel_bridge_in_Sweden_1.zip";
    else, archive(e)="Dataset_steel_bridge_in_Sweden_2.zip"; end
    T=readtable(path,'VariableNamingRule','preserve');
    [Xe,te,audit]=alignRawEvent(T,strainCols,angleCol,tempCol,128);
    X(e,:,:)=Xe; temp(e)=te;
    rawSamples(e)=audit.rawSamples; i2v(e)=audit.i2; i72v(e)=audit.i72;
    j72v(e)=audit.j72; j2v(e)=audit.j2; b0v(e)=audit.b0; b1v(e)=audit.b1;
    if mod(e,8)==0 || e==64, fprintf('Aligned %d/64 raw events\n',e); end
end

AlignmentAudit=table(files,archive,fullPath,rawSamples,temp,i2v,i72v,j72v,j2v,b0v,b1v, ...
    'VariableNames',{'file','archive','path','raw_samples','temperature_C','opening_i2','opening_i72', ...
    'closing_j72','closing_j2','baseline_start','baseline_end_exclusive'});
writetable(AlignmentAudit,fullfile(outDir,'00_alignment_audit.csv'));

SensorMapping=table(strainCols(:),("SG"+string(gauges)).',repmat("usable strain",numel(gauges),1), ...
    'VariableNames',{'raw_channel','gauge','role'});
writetable(SensorMapping,fullfile(outDir,'00_sensor_mapping.csv'));

crackIdx=find(files=="2023_03_09T23_45_25.csv",1);
roles=[repmat("fit",3,1);repmat("calibration",3,1);repmat("evaluation",58,1)];
isCrackPost=(1:64).'>=crackIdx;
EventInventory=table(files,temp,roles,isCrackPost, ...
    'VariableNames',{'file','temperature_C','role_primary','documented_crack_or_post'});
writetable(EventInventory,fullfile(outDir,'00_event_inventory.csv'));

aligned=X; %#ok<NASGU>
save(fullfile(outDir,'00_aligned_opening_cache.mat'),'aligned','files','gauges','temp','-v7');
end

function [Xaligned,temp,audit]=alignRawEvent(T,strainCols,angleCol,tempCol,nq)
raw=zeros(height(T),numel(strainCols)+2);
for j=1:numel(strainCols), raw(:,j)=double(T{:,char(strainCols(j))}); end
raw(:,end-1)=double(T{:,char(angleCol)});
raw(:,end)=double(T{:,char(tempCol)});
[Xaligned,temp,audit]=vanersborg_core('align',raw,nq);
end
