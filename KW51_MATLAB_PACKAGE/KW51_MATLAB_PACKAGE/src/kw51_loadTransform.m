function L=kw51_loadTransform(x,front,axles,loads)
    L=zeros(numel(x),numel(front)); step=x(2)-x(1);
    for a=1:numel(axles)
        u=(front-axles(a)-x(1))/step; ii=floor(u); frac=u-ii;
        for off=0:1
            idx=ii+off+1; ok=idx>=1 & idx<=numel(x);
            if off==0, weights=1-frac; else, weights=frac; end
            cols=find(ok); indices=sub2ind(size(L),idx(ok),cols);
            L(indices)=L(indices)+loads(a)*weights(ok);
        end
    end
end
