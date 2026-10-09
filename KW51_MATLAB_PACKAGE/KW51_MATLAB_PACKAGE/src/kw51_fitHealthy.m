function m=kw51_fitHealthy(F,g,t)
    m.tref=median(t,'omitnan'); m.tscale=max(std(t,1,'omitnan'),1);
    assert(isfinite(m.tref),'Training temperature data are missing.');
    D=kw51_design(m,g,t); m.coef=pinv(D)*F; m.xtxi=pinv(D.'*D);
    R=F-D*m.coef;
    m.scale=sqrt(sum((R-mean(R,1)).^2,1)/(size(F,1)-size(D,2)));
    m.scale=max(m.scale,max(median(m.scale)*.05,1e-4));
    Z=R./m.scale; C=Z.'*Z/max(size(F,1)-size(D,2),1);
    m.alpha=.05; cov=(1-m.alpha)*C+m.alpha*diag(diag(C)+1e-8)+eye(size(F,2))*1e-8;
    [V,E]=eig((cov+cov.')/2,'vector');
    m.W=(V./sqrt(max(E,1e-10)).')*V.';
end
