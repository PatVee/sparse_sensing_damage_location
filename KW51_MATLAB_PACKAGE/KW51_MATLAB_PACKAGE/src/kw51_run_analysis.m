function outputs = kw51_run_analysis(part1Zip, part2Zip, varargin)
% Pipeline coordinator; called by RUN_KW51. See README.md for inputs and outputs.
    if nargin==0
        part1Zip=kw51_chooseInput('*.zip','Select traindata January part 1');
        part2Zip=kw51_chooseInput('*.zip','Select traindata January part 2');
    elseif nargin<2
        error('KW51:Inputs','Supply both ZIP paths, or use no arguments for file selection.');
    end
    thisDir = fileparts(mfilename('fullpath'));

    parser = inputParser;
    parser.FunctionName = mfilename;
    addRequired(parser, 'part1Zip', @(x) ischar(x) || isstring(x));
    addRequired(parser, 'part2Zip', @(x) ischar(x) || isstring(x));
    addParameter(parser, 'OutputRoot', fullfile(thisDir, '..', 'outputs'), ...
        @(x) ischar(x) || isstring(x));
    addParameter(parser, 'CandidateFile', ...
        fullfile(thisDir, '..', 'reference', 'candidate_regions.csv'), ...
        @(x) ischar(x) || isstring(x));
    addParameter(parser, 'FEReferenceFile', ...
        fullfile(thisDir, '..', 'reference', 'base_and_all_gauge_calibration.mat'), ...
        @(x) ischar(x) || isstring(x));
    addParameter(parser, 'RebuildRecords', false, @(x) islogical(x) && isscalar(x));
    addParameter(parser, 'RebuildFE', false, @(x) islogical(x) && isscalar(x));
    addParameter(parser, 'RebuildFeatures', false, @(x) islogical(x) && isscalar(x));
    addParameter(parser, 'RebuildTrials', false, @(x) islogical(x) && isscalar(x));
    addParameter(parser, 'NumberTrials', 300, @kw51_positiveInteger);
    addParameter(parser, 'NumberBootstrap', 2500, @kw51_positiveInteger);
    addParameter(parser, 'NumberBootstrapExtended', 2000, @kw51_positiveInteger);
    addParameter(parser, 'PassageWindows', [16 32], ...
        @(x) isnumeric(x) && isvector(x) && ~isempty(x) && all(isfinite(x) & x>0 & x==floor(x)));
    addParameter(parser, 'RandomSeed', 20260917, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0 && x==floor(x) && x<=2^32-1);
    addParameter(parser, 'MinimumJointRate', 0.50, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    addParameter(parser, 'Verbose', true, @(x) islogical(x) && isscalar(x));
    parse(parser, part1Zip, part2Zip, varargin{:});
    opt = parser.Results;

    part1Zip = char(part1Zip);
    part2Zip = char(part2Zip);
    opt.OutputRoot = char(opt.OutputRoot);
    opt.CandidateFile = char(opt.CandidateFile);
    opt.FEReferenceFile = char(opt.FEReferenceFile);
    opt.PassageWindows=unique(opt.PassageWindows(:).');
    kw51_mustExist(part1Zip, 'raw traindata archive part 1');
    kw51_mustExist(part2Zip, 'raw traindata archive part 2');
    kw51_mustExist(opt.CandidateFile, 'candidate-region CSV');
    kw51_mustExist(opt.FEReferenceFile, 'KW51 FE reference MAT file');
    kw51_requireSignalFunctions();

    dataDir = fullfile(opt.OutputRoot, 'tables');
    cacheDir = fullfile(opt.OutputRoot, 'cache');
    workDir = fullfile(opt.OutputRoot, 'work');
    kw51_ensureDir(opt.OutputRoot);
    kw51_ensureDir(dataDir);
    kw51_ensureDir(cacheDir);
    kw51_ensureDir(workDir);

    recordsCache = fullfile(cacheDir, 'aligned_records.mat');
    feCache = fullfile(cacheDir, 'physical_zone_passages.mat');
    featureCache = fullfile(cacheDir, 'feature_libraries.mat');
    trialCache = fullfile(cacheDir, 'glrt_results.mat');

    % Validate every cache against both source files and numerical settings.
    % Hashing happens only when YOU run this function, not while writing it.
    runSignature=struct('script',kw51_sha256File([mfilename('fullpath') '.m']), ...
        'raw1',kw51_sha256File(part1Zip),'raw2',kw51_sha256File(part2Zip), ...
        'fe',kw51_sha256File(opt.FEReferenceFile),'geometry',kw51_sha256File(opt.CandidateFile), ...
        'nTrials',opt.NumberTrials,'nBootstrap',opt.NumberBootstrap, ...
        'nBootstrapExtended',opt.NumberBootstrapExtended, ...
        'windows',opt.PassageWindows,'seed',opt.RandomSeed);

    sourceFiles=dir(fullfile(thisDir,'*.m'));
    runSignature.sourceNames=string({sourceFiles.name});
    runSignature.sourceHashes=strings(size(runSignature.sourceNames));
    for k=1:numel(sourceFiles)
        runSignature.sourceHashes(k)=kw51_sha256File(fullfile(thisDir,sourceFiles(k).name));
    end

    if opt.Verbose
        fprintf('\nKW51 full localisation pipeline\n');
        fprintf('Output: %s\n', opt.OutputRoot);
    end

    recordsRebuilt=opt.RebuildRecords || ~kw51_validCache(recordsCache,runSignature);
    if recordsRebuilt
        kw51_stageMessage(opt, 1, 'Extracting, auditing and aligning raw train records');
        [aligned, eventInventory, channelQuality, channelAvailability, alignedTable] = ...
            kw51_extract_traindata({part1Zip, part2Zip}, workDir, opt.Verbose);
        save(recordsCache, 'aligned', 'eventInventory', 'channelQuality', ...
            'channelAvailability', 'alignedTable', 'runSignature', '-v7.3');
    else
        kw51_stageMessage(opt, 1, 'Loading cached aligned train records');
        cached = load(recordsCache);
        aligned = cached.aligned;
        eventInventory = cached.eventInventory;
        channelQuality = cached.channelQuality;
        channelAvailability = cached.channelAvailability;
        alignedTable = cached.alignedTable;
    end
    writetable(eventInventory, fullfile(dataDir, 'event_inventory.csv'));
    writetable(channelQuality, fullfile(dataDir, 'channel_event_quality.csv'));
    writetable(channelAvailability, fullfile(dataDir, 'channel_availability.csv'));
    writetable(alignedTable, fullfile(dataDir, 'aligned_records.csv'));
    allocation=table(aligned.events,repmat("fit",62,1),'VariableNames',{'event','subset'});
    allocation.subset(aligned.calibrationIndices)="calibration";
    allocation.subset(aligned.testIndices)="test";
    writetable(allocation,fullfile(dataDir,'record_allocation.csv'));

    upstreamChanged = recordsRebuilt;
    feRebuilt=opt.RebuildFE || ~kw51_validCache(feCache,runSignature);
    if feRebuilt
        kw51_stageMessage(opt, 2, 'Building full-matrix physical-zone strain passages');
        [zones, passage, feAudit] = kw51_build_fe_library( ...
            opt.FEReferenceFile, opt.CandidateFile, opt.Verbose);
        save(feCache, 'zones', 'passage', 'feAudit', 'runSignature', '-v7.3');
    else
        kw51_stageMessage(opt, 2, 'Loading cached physical-zone strain passages');
        cached = load(feCache);
        zones = cached.zones;
        passage = cached.passage;
        feAudit = cached.feAudit;
    end
    writetable(zones, fullfile(dataDir, 'zone_definitions.csv'));
    copyfile(opt.CandidateFile, fullfile(dataDir, 'candidate_regions.csv'), 'f');
    writetable(feAudit, fullfile(dataDir, 'fe_execution_audit.csv'));

    upstreamChanged = upstreamChanged || feRebuilt;
    featuresRebuilt=opt.RebuildFeatures || upstreamChanged || ~kw51_validCache(featureCache,runSignature);
    if featuresRebuilt
        kw51_stageMessage(opt, 3, 'Deriving measured and FE feature libraries');
        [featureLibraries, featureCatalogue, informationTable, qrSets] = ...
            kw51_build_features(aligned, passage, zones, opt.Verbose);
        save(featureCache, 'featureLibraries', 'featureCatalogue', ...
            'informationTable', 'qrSets', 'runSignature', '-v7.3');
    else
        kw51_stageMessage(opt, 3, 'Loading cached feature libraries');
        cached = load(featureCache);
        featureLibraries = cached.featureLibraries;
        featureCatalogue = cached.featureCatalogue;
        informationTable = cached.informationTable;
        qrSets = cached.qrSets; %#ok<NASGU>
    end
    writetable(featureCatalogue, fullfile(dataDir, 'feature_catalogue.csv'));
    writetable(informationTable, fullfile(dataDir, 'zone_feature_information.csv'));
    pivotRows=cell(0,3);
    for name=string(fieldnames(qrSets)).'
        values=qrSets.(char(name));
        for k=1:numel(values), pivotRows(end+1,:)={name,k,values(k)}; end %#ok<AGROW>
    end
    writetable(cell2table(pivotRows,'VariableNames', ...
        {'feature','selection_rank','master_feature_index_1based'}), ...
        fullfile(dataDir,'qr_feature_selection.csv'));

    upstreamChanged = upstreamChanged || featuresRebuilt;
    if opt.RebuildTrials || upstreamChanged || ~kw51_validCache(trialCache,runSignature)
        kw51_stageMessage(opt, 4, 'Calibrating GLRT and running detection/localisation trials');
        [glrtResults, thresholdTable] = kw51_localisation_trials( ...
            aligned, featureLibraries, zones, passage.losses, opt);
        save(trialCache, 'glrtResults', 'thresholdTable', 'runSignature', '-v7.3');
    else
        kw51_stageMessage(opt, 4, 'Loading cached GLRT results');
        cached = load(trialCache);
        glrtResults = cached.glrtResults;
        thresholdTable = cached.thresholdTable;
    end

    kw51_stageMessage(opt, 5, 'Writing source GLRT CSVs and retained-location trial tables');
    report = kw51_write_tables(glrtResults, thresholdTable, zones, ...
        featureCatalogue, dataDir, opt);

    audit = struct;
    audit.created_utc = char(datetime('now', 'TimeZone', 'UTC', ...
        'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
    audit.part1_zip = kw51_absolutePath(part1Zip);
    audit.part2_zip = kw51_absolutePath(part2Zip);
    audit.part1_sha256 = runSignature.raw1;
    audit.part2_sha256 = runSignature.raw2;
    audit.fe_reference_file = kw51_absolutePath(opt.FEReferenceFile);
    audit.fe_reference_sha256 = runSignature.fe;
    audit.candidate_file = kw51_absolutePath(opt.CandidateFile);
    audit.n_unique_events = size(aligned.Y, 1);
    audit.n_active_strain_channels = size(aligned.Y, 2);
    audit.n_phase_samples = size(aligned.Y, 3);
    audit.n_physical_zones = height(zones);
    audit.losses = passage.losses;
    audit.n_feature_families = numel(featureLibraries);
    audit.n_trials_per_state = opt.NumberTrials;
    audit.n_calibration_bootstraps = opt.NumberBootstrap;
    audit.n_calibration_bootstraps_extended = opt.NumberBootstrapExtended;
    audit.passage_windows = opt.PassageWindows;
    audit.matlab_random_seed = opt.RandomSeed;
    audit.minimum_joint_rate = opt.MinimumJointRate;
    audit.n_retained_locations = height(report.retainedLocations);
    audit.sensor_policy='11 measured strain channels; sgBD1718C excluded; no acceleration features in this final study.';
    audit.response_class_note='RMS-ratio response surrogate, not measured track metadata.';
    audit.rng_note = ['MATLAB twister differs from the earlier Python PCG64 stream. ' ...
        'These CSVs contain fresh MATLAB results, not the previously reported counts.'];
    audit.validation_scope='Exploratory semi-synthetic FE damage plus resampled healthy residuals.';
    audit.independent_fit_cal_test_records=[32 14 16];
    audit.threshold_quantile=0.95;
    audit.covariance_diagonal_shrinkage=0.05;
    audit.selection_scope='Frozen 68-zone shortlist; no new exhaustive 1603-element search.';
    kw51_writeText(fullfile(dataDir, 'full_pipeline_run_audit.json'), jsonencode(audit, 'PrettyPrint',true));

    outputs = struct;
    outputs.outputRoot = opt.OutputRoot;
    outputs.dataDirectory = dataDir;
    outputs.alignedRecords = recordsCache;
    outputs.physicalZonePassages = feCache;
    outputs.featureLibraries = featureCache;
    outputs.glrtResults = trialCache;
    outputs.retainedLocations = report.retainedLocations;
    outputs.retainedTrials = report.retainedTrials;

    if opt.Verbose
        fprintf('\nComplete. CSV products are in:\n  %s\n', dataDir);
        fprintf('Retained %d locations with no top-N cap.\n', ...
            height(report.retainedLocations));
    end
end
