function [libraries,catalogue,information,qrSets]=kw51_build_features(aligned,passage,zones,verbose)
    names=["strain_gram","rms_p2p_gram","area_rms_p2p","signed_dct16", ...
        "absolute_rms_p2p_gram","qr_hybrid_16","qr_hybrid_24", ...
        "qr_hybrid_32","qr_hybrid_48","rms_p2p_gram_dct8", ...
        "strain_gram_dct16","rms_p2p_gram_dct16","area_gram_dct8"];
    counts=[66 88 33 176 88 16 24 32 48 176 242 264 187];
    labels=["Normalised strain Gram","RMS/P2P + normalised Gram", ...
        "Signed area + RMS/P2P","Signed DCT16","Absolute RMS/P2P + Gram", ...
        "QR hybrid (16)","QR hybrid (24)","QR hybrid (32)","QR hybrid (48)", ...
        "RMS/P2P + Gram + DCT8","Gram + DCT16", ...
        "RMS/P2P + Gram + DCT16","Area + Gram + DCT8"];
    catalogue=table(names.',labels.',counts.','VariableNames',{'feature','feature_label','n_features'});
    catalogue.feature_components=["66 signed normalised upper-triangular Gram terms"; ...
        "11 centred log-RMS + 11 centred log-P2P + 66 Gram terms"; ...
        "11 normalised signed passage means + 22 RMS/P2P terms"; ...
        "16 orthonormal DCT coefficients for each of 11 gauges"; ...
        "11 absolute log-RMS + 11 absolute log-P2P + 66 absolute Gram terms"; ...
        "16 pivoted-QR terms from RMS/P2P+Gram+DCT8"; ...
        "24 pivoted-QR terms from RMS/P2P+Gram+DCT8"; ...
        "32 pivoted-QR terms from RMS/P2P+Gram+DCT8"; ...
        "48 pivoted-QR terms from RMS/P2P+Gram+DCT8"; ...
        "88 RMS/P2P+Gram terms + 88 signed DCT8 terms"; ...
        "66 Gram terms + 176 signed DCT16 terms"; ...
        "88 RMS/P2P+Gram terms + 176 signed DCT16 terms"; ...
        "33 area/RMS/P2P terms + 66 Gram terms + 88 signed DCT8 terms"];
    master=kw51_featureBatch(aligned.Y,'rms_p2p_gram_dct8',[]);
    sd=std(master(aligned.fitIndices,:),0,1);
    sd=max(sd,median(sd)*.05+1e-8);
    nz=height(zones); ns=numel(passage.losses);
    S10=zeros(2*nz,size(master,2));
    for g=1:2
        base=kw51_featureOne(squeeze(passage.base(g,:,:)),'rms_p2p_gram_dct8',[]);
        for z=1:nz
            h=double(squeeze(passage.response(g,z,1,:,:)));
            S10((g-1)*nz+z,:)=kw51_featureOne(h,'rms_p2p_gram_dct8',[])-base;
        end
    end
    A=S10./sd; A=A-mean(A,1);
    [~,~,P]=qr(A); [~,piv]=max(P,[],1);
    qrSets=struct;
    for k=[16 24 32 48], qrSets.(sprintf('qr_hybrid_%d',k))=piv(1:k); end
    libraries=cell(numel(names),1); infoRows=cell(0,7);
    for fi=1:numel(names)
        name=names(fi); indices=[];
        if startsWith(name,'qr_hybrid_'), indices=qrSets.(char(name)); end
        F=kw51_featureBatch(aligned.Y,name,indices);
        m=kw51_fitHealthy(F(aligned.fitIndices,:),aligned.groups(aligned.fitIndices), ...
            aligned.temperature(aligned.fitIndices));
        Dw=zeros(2,nz,ns,size(F,2));
        for g=1:2
            bf=kw51_featureOne(squeeze(passage.base(g,:,:)),name,indices);
            for z=1:nz
                for si=1:ns
                    h=double(squeeze(passage.response(g,z,si,:,:)));
                    delta=kw51_featureOne(h,name,indices)-bf;
                    Dw(g,z,si,:)=kw51_whiten(m,delta);
                end
            end
        end
        radial=zeros(2,nz); margin=zeros(2,nz); exact=false(2,nz);
        for g=1:2
            X=reshape(Dw(g,:,1,:),nz,[]); norms=sqrt(sum(X.^2,2));
            Xn=X./max(norms,1e-15); corr=Xn*Xn.';
            [~,pred]=max(corr,[],2); exact(g,:)=pred.'==(1:nz);
            corr(1:nz+1:end)=-Inf;
            margin(g,:)=1-max(corr,[],2).'; radial(g,:)=norms.'/sqrt(size(F,2));
        end
        for z=1:nz
            r=min(radial(:,z)); c=min(margin(:,z));
            infoRows(end+1,:)={name,z-1,size(F,2),r,c,all(exact(:,z)),r*max(c,0)};
        end
        libraries{fi}=struct('name',name,'healthy',F,'model',m,'Dw',Dw,'qr_indices',indices);
        if verbose, fprintf('  feature %s (%d terms)\n',name,size(F,2)); end
    end
    information=cell2table(infoRows,'VariableNames',{'feature','zone_index','n_features', ...
        'min_whitened_radial_10pct','min_cosine_margin_10pct', ...
        'noiseless_exact_both_classes_10pct','information_score_10pct'});
    information=join(information,zones,'Keys','zone_index');
end
