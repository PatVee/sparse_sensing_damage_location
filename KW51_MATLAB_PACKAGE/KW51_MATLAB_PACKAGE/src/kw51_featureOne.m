function f=kw51_featureOne(X,name,indices)
    X=double(X); Q=size(X,2);
    r=max(sqrt(mean(X.^2,2)),1e-12); p=max(max(X,[],2)-min(X,[],2),1e-12);
    lr=log(r); lr=lr-mean(lr); lp=log(p); lp=lp-mean(lp);
    Y=X/max(sqrt(mean(X.^2,'all')),1e-12);
    G=(Y*Y.')/Q; GA=(X/1e-5)*(X/1e-5).'/Q;
    gram=zeros(1,66); ag=gram; idx=0;
    for i=1:11
        for j=i:11, idx=idx+1; gram(idx)=G(i,j); ag(idx)=GA(i,j); end
    end
    coeff=dct(Y.').'/sqrt(Q);
    dc8=reshape(coeff(:,1:8).',1,[]); dc16=reshape(coeff(:,1:16).',1,[]);
    area=mean(X,2); area=area/max(norm(area),1e-12);
    rp=[lr.' lp.']; rg=[rp gram]; ar=[area.' rp];
    switch char(name)
        case 'strain_gram', f=gram;
        case 'rms_p2p_gram', f=rg;
        case 'area_rms_p2p', f=ar;
        case 'signed_dct16', f=dc16;
        case 'absolute_rms_p2p_gram', f=[log(r*1e6).' log(p*1e6).' ag];
        case 'rms_p2p_gram_dct8', f=[rg dc8];
        case 'strain_gram_dct16', f=[gram dc16];
        case 'rms_p2p_gram_dct16', f=[rg dc16];
        case 'area_gram_dct8', f=[ar gram dc8];
        otherwise
            assert(startsWith(name,'qr_hybrid_'),'Unknown feature family.');
            f=[rg dc8]; f=f(indices);
    end
end
