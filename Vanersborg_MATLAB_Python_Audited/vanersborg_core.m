function varargout=vanersborg_core(action,varargin)
% Shared numerical implementation used by MATLAB pipeline and Octave audit.
% Raw alignment input columns: usable strains, raw angle, temperature.
switch action
 case 'select', [varargout{1:nargout}]=selectQRReferences(varargin{:});
 case 'select_fixed', [varargout{1:nargout}]=selectQRFixedCount(varargin{:});
 case 'virtual', [varargout{1:nargout}]=fitVirtualTransfer(varargin{:});
 case 'draws', [varargout{1:nargout}]=drawEvents(varargin{:});
 case 'evaluate', [varargout{1:nargout}]=evaluateNumeric(varargin{:});
 case 'scores', [varargout{1:nargout}]=zoneScores(varargin{:});
 case 'features', [varargout{1:nargout}]=areaGramDCT8(varargin{:});
 case 'align', [varargout{1:nargout}]=alignNumeric(varargin{:});
 otherwise, error('Unknown numerical operation: %s',action);
end
end

function [refIdx,cumEnergy]=selectQRReferences(X,gauges,fitIdx,targetEnergy)
snap=[];
for e=fitIdx, snap=[snap;squeeze(X(e,:,:))]; end %#ok<AGROW>
mu=mean(snap,1); sigma=std(snap,0,1); sigma(sigma<=1e-12)=1;
Z=(snap-mu)./sigma;
[~,S,V]=svd(Z,'econ'); sv=diag(S); cumEnergy=cumsum(sv.^2)/sum(sv.^2);
assert(all(isfinite(cumEnergy)), 'Healthy POD energy is nonfinite. Check aligned baseline strains.');
r=find(cumEnergy>=targetEnergy,1,'first');
assert(~isempty(r) && r>=1, 'No POD rank meets the requested energy threshold.');
Phi=V(:,1:r);
[~,~,p]=qr(Phi','vector');
refIdx=p(1:r);
end

function refIdx=selectQRFixedCount(X,fitIdx,k)
snap=[]; for e=fitIdx, snap=[snap;squeeze(X(e,:,:))]; end %#ok<AGROW>
mu=mean(snap,1); sigma=std(snap,0,1); sigma(sigma<=1e-12)=1;
Z=(snap-mu)./sigma; [~,~,V]=svd(Z,'econ'); Phi=V(:,1:k);
[~,~,p]=qr(Phi','vector'); refIdx=p(1:k);
end

function [Xv,model]=fitVirtualTransfer(X,refIdx,withheldIdx,fitIdx,alpha)
xt=[]; yt=[];
for e=fitIdx
    xt=[xt;reshape(X(e,:,refIdx),size(X,2),numel(refIdx))]; %#ok<AGROW>
    yt=[yt;reshape(X(e,:,withheldIdx),size(X,2),numel(withheldIdx))]; %#ok<AGROW>
end
xm=mean(xt,1); xs=std(xt,0,1); xs(xs<=1e-12)=1;
ym=mean(yt,1); ys=std(yt,0,1); ys(ys<=1e-12)=1;
xz=(xt-xm)./xs; yz=(yt-ym)./ys;
B=(xz'*xz+alpha*eye(size(xz,2)))\(xz'*yz);
Xv=zeros(size(X));
for e=1:size(X,1)
    xr=reshape(X(e,:,refIdx),size(X,2),numel(refIdx)); yp=((xr-xm)./xs)*B.*ys+ym;
    Xv(e,:,refIdx)=xr; Xv(e,:,withheldIdx)=yp;
end
model=struct('xm',xm,'xs',xs,'ym',ym,'ys',ys,'B',B,'alpha',alpha);
end

function [zoneNames,scores]=zoneScores(Xin,availableGauges,trainIdx)
zoneNames={'Z01_SG1_SG2','Z02_SG3','Z03_SG7','Z04_SG8','Z05_SG9_SG10','Z06_SG11','Z07_SG13_SG14','Z08_SG15_SG16'};
zoneMembers={ [1 2],3,7,8,[9 10],11,[13 14],[15 16] };
assert(~isempty(Xin) && size(Xin,3)==numel(availableGauges), ...
    'Strain array must contain events x samples x available gauges.');
assert(~isempty(availableGauges) && all(ismember(availableGauges,[1 2 3 7 8 9 10 11 13 14 15 16])), ...
    'Available gauges must be physical SG numbers, not QR column positions.');
keepNames=cell(0,1); cols={};
for z=1:numel(zoneNames)
    members=intersect(zoneMembers{z},availableGauges,'stable');
    if isempty(members), continue; end
    idx=zeros(1,numel(members)); for j=1:numel(members), idx(j)=find(availableGauges==members(j),1); end
    F=zeros(size(Xin,1),numel(areaGramDCT8(reshape(Xin(1,:,idx),size(Xin,2),numel(idx)))));
    for e=1:size(Xin,1), F(e,:)=areaGramDCT8(reshape(Xin(e,:,idx),size(Xin,2),numel(idx))); end
    sd=std(F(trainIdx,:),0,1); pos=sd(sd>1e-12);
    if isempty(pos), floorv=1+1e-8; else, floorv=median(pos)*0.05+1e-8; end
    sd=max(sd,floorv); A=F(trainIdx,:)./sd; Y=F./sd;
    best=inf(size(Y,1),1);
    for a=1:size(A,1)
        err=mean((Y-A(a,:)).^2,2); best=min(best,err);
    end
    pairs=nchoosek(1:size(A,1),2);
    for p=1:size(pairs,1)
        a=A(pairs(p,1),:); b=A(pairs(p,2),:); v=b-a; den=max(v*v',1e-30);
        t=((Y-a)*v')/den; t=min(max(t,0),1); pred=a+t.*v;
        err=mean((Y-pred).^2,2); best=min(best,err);
    end
    keepNames(end+1,1)=zoneNames(z); %#ok<AGROW>
    cols{end+1}=best; %#ok<AGROW>
end
zoneNames=keepNames; scores=zeros(size(Xin,1),numel(cols));
for i=1:numel(cols), scores(:,i)=cols{i}; end
end

function f=areaGramDCT8(arr)
arr=double(arr); if isvector(arr), arr=arr(:); end
Q=size(arr,1); m=size(arr,2); eps0=1e-12;
meanv=mean(arr,1); area=meanv/max(norm(meanv),eps0);
rmsv=sqrt(mean(arr.^2,1)); p2p=max(arr,[],1)-min(arr,[],1);
lr=log(max(rmsv,eps0)); lr=lr-mean(lr); lp=log(max(p2p,eps0)); lp=lp-mean(lp);
scale=sqrt(max(mean(arr(:).^2),eps0)); Y=arr/scale; G=(Y'*Y)/Q;
mask=tril(true(m)); GT=G.'; gram=GT(mask).';
D=dctMatrix(8,Q); coeff=(D*Y)/sqrt(Q); dct8=reshape(coeff,1,[]);
f=[area,lr,lp,gram,dct8];
end

function D=dctMatrix(kmax,N)
% First kmax rows of orthonormal DCT-II matrix.
n=0:N-1; D=zeros(kmax,N);
for k=0:kmax-1
    if k==0, a=sqrt(1/N); else, a=sqrt(2/N); end
    D(k+1,:)=a*cos(pi/N*(n+0.5)*k);
end
end

function [Xaligned,temp,audit]=alignNumeric(raw,nq)
n=size(raw,2)-2; N=size(raw,1); X=zeros(N,n);
for j=1:n, X(:,j)=cleanInterp(raw(:,j)); end
angle=-cleanInterp(raw(:,n+1));
temp=median(cleanInterp(raw(:,n+2)));

i2=firstCross(angle,2,1,true); assert(~isempty(i2),'No opening 2-degree crossing.');
i72=firstCross(angle,72,i2,true); assert(~isempty(i72),'No opening 72-degree crossing.');
hi=find(angle>=72); assert(~isempty(hi),'No high-angle samples.'); j72=hi(end);
j2=firstCross(angle,2,j72,false);
if isempty(j2)
    low=find((1:numel(angle)).'>j72 & angle<=4,1,'first');
    assert(~isempty(low),'No closing low-angle segment.'); j2=low;
end

% 200 Hz source: 7.5 s to 0.5 s before opening starts.
b0=max(1,i2-1500); b1=max(2,i2-100);
base=median(X(b0:b1-1,:),1); X=X-base;
op=interpBranch(angle,X,i2,i72,true,nq);
cl=interpBranch(angle,X,j72,j2,false,nq);
Xaligned=[op;cl];
audit=struct('rawSamples',N,'i2',i2,'i72',i72,'j72',j72,'j2',j2,'b0',b0,'b1',b1);
end

function v=cleanInterp(v)
v=double(v(:)); bad=~isfinite(v) | v<=-9e5 | abs(v)>5000;
if all(bad), v=zeros(size(v)); return; end
if any(bad)
    good=find(~bad); badIdx=find(bad);
    if numel(good)==1
        v(bad)=v(good(1));
    else
        inside=badIdx(badIdx>=good(1) & badIdx<=good(end));
        if ~isempty(inside), v(inside)=interp1(good,v(good),inside,'linear'); end
        v(badIdx(badIdx<good(1)))=v(good(1));
        v(badIdx(badIdx>good(end)))=v(good(end));
    end
end
end

function idx=firstCross(a,level,startIdx,isUp)
a=a(:); startIdx=max(1,startIdx);
if isUp
 hit=a(startIdx:end-1)<level & a(startIdx+1:end)>=level;
else
 hit=a(startIdx:end-1)>level & a(startIdx+1:end)<=level;
end
k=find(hit,1,'first');
if isempty(k), idx=[]; else, idx=startIdx+k; end
end

function out=interpBranch(angle,X,i0,i1,ascending,nq)
a=angle(i0:i1); y=X(i0:i1,:);
if ascending
    target=linspace(2,72,nq).'; targetAsc=target;
else
    target=linspace(72,2,nq).'; a=flipud(a); y=flipud(y); targetAsc=flipud(target);
end
[a,ord]=sort(a); y=y(ord,:); % Match NumPy round-to-even at exact half-centidegree ties.
scaled=a*100; rounded=round(scaled); ties=(scaled-floor(scaled))==0.5;
rounded(ties)=2*round(scaled(ties)/2); ar=rounded/100;
[u,~,ic]=unique(ar,'sorted'); ym=zeros(numel(u),size(y,2));
counts=accumarray(ic,1,[numel(u),1]);
for j=1:size(y,2), ym(:,j)=accumarray(ic,y(:,j),[numel(u),1])./counts; end
out=zeros(nq,size(y,2));
% NumPy interp uses endpoint values outside the observed angle range.
% Clamping avoids NaNs when the first crossing sample lies above 2 degrees.
query=min(max(targetAsc,u(1)),u(end));
for j=1:size(y,2)
    if numel(u)==1, out(:,j)=ym(1,j);
    else, out(:,j)=interp1(u,ym(:,j),query,'linear'); end
end
if ~ascending, out=flipud(out); end
end

function C=evaluateNumeric(Xin,gauges,fitIdx,calIdx,damageIdx,crackIdx,truthZone)
[names,scores]=zoneScores(Xin,gauges,fitIdx);
nEvents=size(Xin,1);
assert(size(scores,1)==nEvents && size(scores,2)>0, ...
    'Expected one score row per event and at least one represented zone.');
assert(all(isfinite(scores(:))), 'Zone scores contain NaN or Inf.');
idx=[fitIdx(:);calIdx(:);damageIdx(:);crackIdx];
assert(all(idx>=1 & idx<=nEvents & idx==fix(idx)), ...
    'Requested event indices exceed the available event array.');
cal=scores(calIdx,:); threshold=max(cal(:));
[maxScore,pred]=max(scores,[],2);
detected=maxScore>threshold;
correct=reshape(strcmp(names(pred),truthZone),[],1);
C=struct(); C.zoneNames=names; C.scores=scores;
C.threshold=threshold; C.maxScore=maxScore; C.pred=pred;
C.detected=detected; C.correct=correct;
C.detectedOf7=sum(detected(damageIdx)); C.correctOf7=sum(correct(damageIdx));
C.jointOf7=sum(detected(damageIdx)&correct(damageIdx));
C.medianRatio=median(maxScore(damageIdx)/threshold);
C.crackDetected=detected(crackIdx); C.crackPredicted=names{pred(crackIdx)};
C.preDocumentationAlerts=sum(detected(7:crackIdx-1));
end

function draws=drawEvents(eventIdx,nTrials,window)
assert(~isempty(eventIdx),'Cannot resample an empty event list.');
draws=reshape(eventIdx(randi(numel(eventIdx),nTrials,window)),nTrials,window);
end
