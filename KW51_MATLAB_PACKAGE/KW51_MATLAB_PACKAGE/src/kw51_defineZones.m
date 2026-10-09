function zones=kw51_defineZones(R)
    rows=cell(0,9);
    names=sort(unique(R.memberName(R.family=="diagonal")));
    for name=names(:).'
        rows(end+1,:)=kw51_zoneRow(name,"diagonal",R(R.memberName==name & R.family=="diagonal",:),size(rows,1));
    end
    sides=["A","C"]; ys=[0 11.75];
    for j=1:2
        q=R(R.family=="main_girder" & abs(R.y_m-ys(j))<1e-7,:);
        for c=12:12:108
            qq=q(q.x_m>=c-6 & q.x_m<c+6,:);
            if ~isempty(qq)
                rows(end+1,:)=kw51_zoneRow(sprintf('main_girder_%s_X%03d',sides(j),c),"main_girder",qq,size(rows,1));
            end
        end
    end
    tx=[30.625 34.558 46.358 50.291 58.158 81.758 85.691 97.391 101.225];
    names=sort(unique(R.memberName(R.family=="deck_transverse")));
    for name=names(:).'
        q=R(R.family=="deck_transverse" & R.memberName==name,:);
        if min(abs(tx-q.x_m(1)))<.01
            rows(end+1,:)=kw51_zoneRow(name,"deck_transverse",q,size(rows,1));
        end
    end
    ys=[1.75 4.75 10]; centres={[36 48 84 96],60,[36 48 84 96]};
    for j=1:3
        q=R(R.family=="deck_longitudinal" & abs(R.y_m-ys(j))<1e-7,:);
        for c=centres{j}
            qq=q(q.x_m>=c-6 & q.x_m<c+6,:);
            if ~isempty(qq)
                name=sprintf('longitudinal_Y%05.2f_X%03d',ys(j),c);
                rows(end+1,:)=kw51_zoneRow(name,"deck_longitudinal",qq,size(rows,1));
            end
        end
    end
    zones=cell2table(rows,'VariableNames',{'zone_index','zone_name','family', ...
        'x_m','y_m','z_m','physical_length_m','n_elements','element_ids'});
    zones.zone_name=string(zones.zone_name); zones.family=string(zones.family);
    zones.element_ids=string(zones.element_ids);
    assert(height(zones)==68,'Expected the original 68 physical zones.');
end
