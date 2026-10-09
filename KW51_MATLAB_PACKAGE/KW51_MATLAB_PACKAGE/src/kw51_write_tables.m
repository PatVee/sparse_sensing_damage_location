function report=kw51_write_tables(R,thresholds,zones,catalogue,dataDir,opt)
    writetable(R,fullfile(dataDir,'glrt_zone_results.csv'));
    writetable(thresholds,fullfile(dataDir,'glrt_thresholds.csv'));
    R=renamevars(R,{'detection','localization','both','severity_exact'}, ...
        {'detection_rate','localisation_rate','detected_and_correctly_localised_rate','severity_exact_rate'});
    R=join(R,catalogue,'Keys','feature');
    R.severity_pct=round(100*R.loss); R.total_trials=repmat(opt.NumberTrials,height(R),1);
    rates={'detection_rate','localisation_rate','detected_and_correctly_localised_rate'};
    stems={'detected','correctly_localised','detected_and_correctly_localised'};
    for k=1:3
        count=round(R.(rates{k})*opt.NumberTrials); [lo,hi]=kw51_wilson(count,opt.NumberTrials);
        R.([stems{k} '_trials'])=count;
        R.([stems{k} '_wilson95_low'])=lo; R.([stems{k} '_wilson95_high'])=hi;
    end
    R.detected_wrong_location_trials=R.detected_trials-R.detected_and_correctly_localised_trials;
    R.correctly_localised_without_detection_trials=R.correctly_localised_trials-R.detected_and_correctly_localised_trials;
    R.not_detected_trials=opt.NumberTrials-R.detected_trials;
    R.localisation_given_detection=R.detected_and_correctly_localised_trials./R.detected_trials;
    R.localisation_given_detection(R.detected_trials==0)=NaN;
    rankingRows=cell(0,9); names=catalogue.feature;
    for n=opt.PassageWindows(:).'
        for z=zones.zone_index(:).'
            for name=names(:).'
                q=sortrows(R(R.n_passages==n & R.zone_index==z & R.feature==name,:),'loss');
                rankingRows(end+1,:)={n,z,name,min(q.detected_and_correctly_localised_rate), ...
                    q.detected_and_correctly_localised_rate(1),mean(q.detected_and_correctly_localised_rate), ...
                    min(q.detection_rate),min(q.localisation_rate),q.feature_label(1)};
            end
        end
    end
    ranking=cell2table(rankingRows,'VariableNames',{'n_passages','zone_index','feature', ...
        'minimum_joint_rate_10_50','joint_rate_10pct','mean_joint_rate_10_50', ...
        'minimum_detection_rate_10_50','minimum_localisation_rate_10_50','feature_label'});
    ranking=sortrows(ranking,{'n_passages','zone_index','minimum_joint_rate_10_50', ...
        'joint_rate_10pct','mean_joint_rate_10_50','minimum_detection_rate_10_50', ...
        'minimum_localisation_rate_10_50','feature'}, ...
        {'ascend','ascend','descend','descend','descend','descend','descend','ascend'});
    ranking.robust_rank=repmat((1:numel(names)).',height(ranking)/numel(names),1);
    ranking.is_primary_robust_best=ranking.robust_rank==1;
    ranking.is_pareto_feature=false(height(ranking),1);
    stateBest=false(height(R),1);
    for n=opt.PassageWindows(:).'
        for z=zones.zone_index(:).'
            ri=find(ranking.n_passages==n & ranking.zone_index==z);
            values=zeros(numel(ri),5);
            for k=1:numel(ri)
                ix=find(R.n_passages==n & R.zone_index==z & R.feature==ranking.feature(ri(k)));
                [~,ord]=sort(R.loss(ix)); ix=ix(ord);
                values(k,:)=R.detected_and_correctly_localised_rate(ix).';
            end
            for k=1:numel(ri)
                dominated=any(all(values>=values(k,:)-1e-12,2) & any(values>values(k,:)+1e-12,2));
                ranking.is_pareto_feature(ri(k))=~dominated;
            end
            q=R.n_passages==n & R.zone_index==z;
            for severity=[10 20 30 40 50]
                ix=find(q & R.severity_pct==severity);
                stateBest(ix)=abs(R.detected_and_correctly_localised_rate(ix)- ...
                    max(R.detected_and_correctly_localised_rate(ix)))<=1e-12;
            end
        end
    end
    primary=ranking(ranking.robust_rank==1,:);
    primary.meets_retention_rule=primary.minimum_joint_rate_10_50>=opt.MinimumJointRate-1e-12;
    retainedIDs=unique(primary.zone_index(primary.meets_retention_rule));
    retained=zones(ismember(zones.zone_index,retainedIDs),:);
    retained.best_minimum_joint_rate=zeros(height(retained),1);
    retained.qualifying_passage_windows=strings(height(retained),1);
    retained.retention_threshold=repmat(opt.MinimumJointRate,height(retained),1);
    for n=opt.PassageWindows(:).'
        suffix=sprintf('%dpass',n);
        retained.(['best_feature_' suffix])=strings(height(retained),1);
        retained.(['best_feature_label_' suffix])=strings(height(retained),1);
        retained.(['minimum_joint_rate_' suffix])=nan(height(retained),1);
        retained.(['joint_rate_10pct_' suffix])=nan(height(retained),1);
        retained.(['mean_joint_rate_' suffix])=nan(height(retained),1);
    end
    for i=1:height(retained)
        q=primary(primary.zone_index==retained.zone_index(i),:);
        retained.best_minimum_joint_rate(i)=max(q.minimum_joint_rate_10_50);
        retained.qualifying_passage_windows(i)=strjoin(string(q.n_passages(q.meets_retention_rule)),';');
        for j=1:height(q)
            suffix=sprintf('%dpass',q.n_passages(j));
            retained.(['best_feature_' suffix])(i)=q.feature(j);
            retained.(['best_feature_label_' suffix])(i)=q.feature_label(j);
            retained.(['minimum_joint_rate_' suffix])(i)=q.minimum_joint_rate_10_50(j);
            retained.(['joint_rate_10pct_' suffix])(i)=q.joint_rate_10pct(j);
            retained.(['mean_joint_rate_' suffix])(i)=q.mean_joint_rate_10_50(j);
        end
    end
    retained=sortrows(retained,{'best_minimum_joint_rate','zone_index'},{'descend','ascend'});
    retained.retained_location_number=(1:height(retained)).';
    selected=innerjoin(R,primary(:,{'n_passages','zone_index','feature','robust_rank','minimum_joint_rate_10_50'}), ...
        'Keys',{'n_passages','zone_index','feature'});
    selected=sortrows(selected,{'zone_index','n_passages','severity_pct'});
    retainedTrials=selected(ismember(selected.zone_index,retainedIDs),:);
    writetable(R,fullfile(dataDir,'all_location_feature_trials.csv'));
    writetable(ranking,fullfile(dataDir,'feature_rankings.csv'));
    writetable(R(stateBest,:),fullfile(dataDir,'all_best_features_by_location_and_severity.csv'));
    writetable(ranking(ranking.is_pareto_feature,:),fullfile(dataDir,'pareto_features_by_location.csv'));
    writetable(primary,fullfile(dataDir,'robust_best_feature_for_every_location.csv'));
    writetable(selected,fullfile(dataDir,'best_feature_trials_for_every_location.csv'));
    writetable(retained,fullfile(dataDir,'retained_best_locations.csv'));
    writetable(retainedTrials,fullfile(dataDir,'retained_location_trials.csv'));
    perDir=fullfile(dataDir,'per_location'); kw51_ensureDir(perDir);
    for z=zones.zone_index(:).'
        writetable(R(R.zone_index==z,:),fullfile(perDir,sprintf('zone_%02d_all_features.csv',z)));
    end
    report=struct('retainedLocations',retained,'retainedTrials',retainedTrials);
end
