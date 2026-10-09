function F=kw51_featureBatch(X,name,indices)
    first=kw51_featureOne(squeeze(X(1,:,:)),name,indices);
    F=zeros(size(X,1),numel(first)); F(1,:)=first;
    for e=2:size(X,1), F(e,:)=kw51_featureOne(squeeze(X(e,:,:)),name,indices); end
end
