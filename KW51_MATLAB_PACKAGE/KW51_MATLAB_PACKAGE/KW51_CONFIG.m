function cfg=KW51_CONFIG(packageRoot)
% Edit analysis settings here. No absolute model paths are required.
if nargin==0, packageRoot=fileparts(mfilename('fullpath')); end
cfg.archives={fullfile(packageRoot,'input','traindata_201901_part1.zip'), ...
    fullfile(packageRoot,'input','traindata_201901_part2.zip')};
cfg.outputRoot=fullfile(packageRoot,'outputs');
cfg.numberTrials=300;
cfg.numberBootstrap=2500;
cfg.numberBootstrapExtended=2000;
cfg.passageWindows=[16 32];
cfg.randomSeed=20260917;
cfg.minimumJointRate=0.50; % Across all five loss states, no top-N cap.
cfg.rebuild=false; % Valid caches are reused; changed input/source/settings invalidate them.
cfg.plotAllLocations=true; % All 68 zone probability figures, plus retained-zone spatial panels.
end
