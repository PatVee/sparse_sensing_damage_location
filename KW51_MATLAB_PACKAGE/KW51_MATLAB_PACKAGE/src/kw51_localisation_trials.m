function [results,thresholds]=kw51_localisation_trials(aligned,libraries,zones,losses,opt)
    rng(opt.RandomSeed,'twister');
    nz=height(zones); ns=numel(losses); nat=nz*ns;
    zone=repelem((0:nz-1).',ns); sev=repmat((1:ns).',nz,1);
    rows=cell(0,10); thresholdRows=cell(0,6);
    trialDir=fullfile(opt.OutputRoot,'trials'); kw51_ensureDir(trialDir);
    for fi=1:numel(libraries)
        lib=libraries{fi}; p=size(lib.healthy,2); atoms=cell(2,1);
        nCal=opt.NumberBootstrap;
        if ismember(lib.name,["rms_p2p_gram_dct8","strain_gram_dct16", ...
                "rms_p2p_gram_dct16","area_gram_dct8"])
            nCal=opt.NumberBootstrapExtended;
        end
        for g=1:2
            atoms{g}=reshape(permute(lib.Dw(g,:,:,:),[3 2 4 1]),nat,p);
        end
        [Hc,~,~]=kw51_scoreContributions(lib,aligned,aligned.calibrationIndices,atoms);
        [Ht,it,gt]=kw51_scoreContributions(lib,aligned,aligned.testIndices,atoms);
        crossA=2*(atoms{1}*atoms{1}.'); crossB=2*(atoms{2}*atoms{2}.');
        for n=opt.PassageWindows(:).'
            maxima=zeros(nCal,1);
            calibrationSampledEvents=zeros(nCal,n,'uint16');
            for b=1:nCal
                ix=randi(size(Hc,1),n,1); maxima(b)=max(sum(Hc(ix,:),1));
                calibrationSampledEvents(b,:)=uint16(aligned.calibrationIndices(ix).');
            end
            sorted=sort(maxima); threshold=sorted(ceil(.95*(numel(sorted)-1))+1);
            thresholdRows(end+1,:)={lib.name,n,threshold,mean(maxima), ...
                quantile(maxima,.99),nCal};
            predictedAtom=zeros(opt.NumberTrials,nat,'uint16');
            detected=false(opt.NumberTrials,nat); maxScore=zeros(opt.NumberTrials,nat);
            sampledEvents=zeros(opt.NumberTrials,n,'uint16');
            healthyDetected=false(opt.NumberTrials,1);
            for b=1:opt.NumberTrials
                ix=randi(size(Ht,1),n,1); H=sum(Ht(ix,:),1);
                wg=[sum(1./it(ix(gt(ix)==0)).^2),sum(1./it(ix(gt(ix)==1)).^2)];
                score=H+wg(1)*crossA+wg(2)*crossB;
                [mx,pr]=max(score,[],2);
                predictedAtom(b,:)=uint16(pr.'-1); detected(b,:)=mx.'>threshold;
                maxScore(b,:)=mx.';
                sampledEvents(b,:)=uint16(aligned.testIndices(ix).');
                healthyDetected(b)=max(H)>threshold;
            end
            truthZone=zone.'; truthSeverityIndex=sev.';
            correct=zone(double(predictedAtom)+1)==truthZone;
            severityCorrect=sev(double(predictedAtom)+1)==truthSeverityIndex;
            dc=sum(detected,1); lc=sum(correct,1); bc=sum(detected & correct,1); sc=sum(severityCorrect,1);
            for a=1:nat
                rows(end+1,:)={lib.name,n,zone(a),losses(sev(a)), ...
                    dc(a)/opt.NumberTrials,lc(a)/opt.NumberTrials,bc(a)/opt.NumberTrials, ...
                    sc(a)/opt.NumberTrials,threshold,mean(healthyDetected)};
            end
            feature=lib.name; n_passages=n; threshold_95=threshold; %#ok<NASGU>
            atomZone=zone; atomLoss=losses(sev); eventIDs=aligned.events; %#ok<NASGU>
            calibrationMaxScore=maxima; %#ok<NASGU>
            save(fullfile(trialDir,sprintf('%s_%dpass.mat',lib.name,n)), ...
                'predictedAtom','detected','maxScore','sampledEvents','truthZone', ...
                'truthSeverityIndex','feature','n_passages','threshold_95','healthyDetected', ...
                'atomZone','atomLoss','eventIDs','calibrationMaxScore','calibrationSampledEvents');
            if opt.Verbose, fprintf('  trials: %s, %d passages\n',lib.name,n); end
        end
    end
    results=cell2table(rows,'VariableNames',{'feature','n_passages','zone_index','loss', ...
        'detection','localization','both','severity_exact','threshold_95','healthy_alert_rate'});
    results=join(results,zones,'Keys','zone_index');
    thresholds=cell2table(thresholdRows,'VariableNames',{'feature','n_passages', ...
        'threshold_95','cal_max_mean','cal_max_q99','n_calibration_bootstraps'});
end
