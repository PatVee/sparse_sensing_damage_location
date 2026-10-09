function [H,infl,g]=kw51_scoreContributions(lib,aligned,indices,atoms)
    m=lib.model; g=aligned.groups(indices); t=aligned.temperature(indices);
    D=kw51_design(m,g,t); raw=lib.healthy(indices,:)-D*m.coef;
    infl=sqrt(1+sum((D*m.xtxi).*D,2)); r=kw51_whiten(m,raw)./infl;
    H=zeros(numel(indices),size(atoms{1},1));
    for e=1:numel(indices)
        A=atoms{g(e)+1}/infl(e);
        H(e,:)=2*(r(e,:)*A.')-sum(A.^2,2).';
    end
end
