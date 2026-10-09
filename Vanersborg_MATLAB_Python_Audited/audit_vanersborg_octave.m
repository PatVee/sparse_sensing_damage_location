function audit_vanersborg_octave(fixtureDir,outDir)
% Runs the exact numerical core used by the MATLAB pipeline.
if ~exist(outDir,'dir'), mkdir(outDir); end
R=load(fullfile(fixtureDir,'reference.mat'));
fid=fopen(fullfile(outDir,'octave_audit.txt'),'w');
fprintf(fid,'Runtime: %s\n',version);
X=zeros(size(R.aligned)); temps=zeros(64,1);
for e=1:64
 D=load(fullfile(fixtureDir,sprintf('raw_%02d.mat',e)));
 [v,temps(e)]=vanersborg_core('align',D.raw,128); X(e,:,:)=v;
 assert(all(isfinite(v(:))));
 if mod(e,8)==0, fprintf('Octave aligned %d/64 raw events\n',e); end
end
save('-mat7-binary',fullfile(outDir,'octave_aligned.mat'),'X','temps');
alignmentError=max(abs(X(:)-R.aligned(:)));
fprintf(fid,'Alignment max absolute difference: %.12g\n',alignmentError);
assert(alignmentError<1e-8,'Raw alignment differs from Python.');
assert(max(abs(temps(:)-R.temp(:)))<1e-10);
for channels=[1 2 5 12]
 f=vanersborg_core('features',R.(sprintf('featureInput%d',channels)));
 expected=R.(sprintf('featureExpected%d',channels));
 assert(max(abs(f(:)-expected(:)))<1e-12);
end
fprintf(fid,'Feature checks for 1, 2, 5, 12 channels: PASS\n');
[ri,energy]=vanersborg_core('select',X,R.gauges,1:3,0.999);
assert(isequal(ri(:),R.refIdx(:)));
assert(max(abs(energy(:)-R.energy(:)))<1e-12);
wi=setdiff(1:12,ri,'stable');
xv=vanersborg_core('virtual',X,ri,wi,1:3,0.01);
virtualError=max(abs(xv(:)-R.virtual(:)));
fprintf(fid,'Virtual-strain max absolute difference: %.12g\n',virtualError);
assert(virtualError<1e-8);
xx={X,X(:,:,ri),xv}; gg={R.gauges,R.gauges(ri),R.gauges};
counts=zeros(3,3);
for k=1:3
 C=vanersborg_core('evaluate',xx{k},gg{k},1:3,4:6,R.damageIdx,R.damageIdx(1),'Z05_SG9_SG10');
 assert(numel(C.detected)==64 && numel(C.correct)==64);
 expected=R.(sprintf('scores%d',k));
 relativeError=max(abs(C.scores(:)-expected(:))./max(1,abs(expected(:))));
 fprintf(fid,'Configuration %d score scaled error: %.12g\n',k,relativeError);
 assert(relativeError<1e-7);
 counts(k,:)=[C.detectedOf7 C.correctOf7 C.jointOf7];
end
assert(isequal(counts,[7 7 7;2 7 2;7 7 7]));
fprintf(fid,'Primary detection/correct/joint counts: 7/7/7; 2/7/2; 7/7/7 PASS\n');
comb=nchoosek(1:6,3); splitCounts=zeros(60,3); splitRefs=zeros(20,5);
for i=1:20
 tr=comb(i,:); ca=setdiff(1:6,tr);
 ri=vanersborg_core('select_fixed',X,tr,5); wi=setdiff(1:12,ri,'stable');
 splitRefs(i,:)=R.gauges(ri);
 xv=vanersborg_core('virtual',X,ri,wi,tr,0.01);
 xx={X,X(:,:,ri),xv}; gg={R.gauges,R.gauges(ri),R.gauges};
 for k=1:3
  C=vanersborg_core('evaluate',xx{k},gg{k},tr,ca,R.damageIdx,R.damageIdx(1),'Z05_SG9_SG10');
  splitCounts(3*(i-1)+k,:)=[C.detectedOf7 C.correctOf7 C.jointOf7];
 end
end
% Catch the empty-zone and bad-index conditions before counting.
failed=false; try, vanersborg_core('evaluate',X(:,:,[]),[],1:3,4:6,R.damageIdx,58,'Z05_SG9_SG10'); catch, failed=true; end
assert(failed);
failed=false; try, vanersborg_core('evaluate',X,R.gauges,1:3,4:6,65,58,'Z05_SG9_SG10'); catch, failed=true; end
assert(failed);
fprintf(fid,'All 20 baseline partitions executed; invalid-input regression checks: PASS\n');
rng(818);
for window=[1 3 5]
 draws=vanersborg_core('draws',4:6,5000,window);
 assert(isequal(size(draws),[5000 window]));
 assert(numel(unique(draws(:,1)))==3);
 draws=vanersborg_core('draws',R.damageIdx,1000,window);
 assert(isequal(size(draws),[1000 window]));
 assert(all(ismember(draws(:),R.damageIdx)));
end
fprintf(fid,'Bootstrap draw shapes for windows 1, 3, 5: PASS\n');
fclose(fid);
save('-mat7-binary',fullfile(outDir,'octave_audit_results.mat'),'counts','splitCounts','splitRefs','alignmentError','virtualError');
fprintf('Octave numerical audit completed.\n');
end
