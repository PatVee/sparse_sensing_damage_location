function [aligned, eventInventory, quality, availability, alignedTable] = ...
        kw51_extract_traindata(zipPaths, workDir, verbose)
% Exact MATLAB analogue of python_reference/process_records.py followed by
% the alignment portion of build_features.py.

    labels = ["sgBD1011A","sgBD1415A","sgBD1718A","sgBD1718C", ...
        "sgBD2324AB","sgBD2324AT","sgBD2324C","sgBD2728A", ...
        "sgDI20ALB","sgDI20ALL","sgDI23ALB","sgDI23ALL"];
    active = [1 2 3 5 6 7 8 9 10 11 12];
    targetFs = 50;
    nPhase = 128;

    extractRoot = fullfile(workDir, 'extracted_archives');
    recordsDir = fullfile(workDir, 'processed_records');
    if isfolder(extractRoot)
        rmdir(extractRoot, 's');
    end
    if isfolder(recordsDir)
        rmdir(recordsDir, 's');
    end
    kw51_ensureDir(extractRoot);
    kw51_ensureDir(recordsDir);

    candidateFiles = struct('name', {}, 'folder', {}, 'date', {}, 'bytes', {}, ...
        'isdir', {}, 'datenum', {});
    for k = 1:numel(zipPaths)
        destination = fullfile(extractRoot, sprintf('part%d', k));
        kw51_ensureDir(destination);
        unzip(zipPaths{k}, destination);
        found = dir(fullfile(destination, '**', '*.mat'));
        candidateFiles = [candidateFiles; found(:)]; %#ok<AGROW>
    end
    if isempty(candidateFiles)
        error('KW51:NoRecords', 'No MAT records were found in the two archives.');
    end

    % Deduplicate using basename plus SHA-256, preserving chronological order.
    seen = containers.Map('KeyType', 'char', 'ValueType', 'char');
    uniqueFiles = strings(0,1);
    uniqueNames = strings(0,1);
    uniqueTimes = NaT(0,1);
    uniqueHashes = strings(0,1);
    for k = 1:numel(candidateFiles)
        name = candidateFiles(k).name;
        path = fullfile(candidateFiles(k).folder, name);
        digest = kw51_sha256File(path);
        if isKey(seen, name)
            if ~strcmp(seen(name), digest)
                error('KW51:ConflictingDuplicate', ...
                    'Duplicate filename has different bytes: %s', name);
            end
            continue;
        end
        seen(name) = digest;
        token = regexp(name, 'traindata_(\d{8}_\d{6})\.mat$', 'tokens', 'once');
        if isempty(token)
            continue;
        end
        uniqueFiles(end+1,1) = string(path); %#ok<AGROW>
        uniqueNames(end+1,1) = string(name); %#ok<AGROW>
        uniqueTimes(end+1,1) = datetime(token{1}, 'InputFormat', 'yyyyMMdd_HHmmss'); %#ok<AGROW>
        uniqueHashes(end+1,1) = string(digest); %#ok<AGROW>
    end
    [uniqueTimes, order] = sort(uniqueTimes);
    uniqueFiles = uniqueFiles(order);
    uniqueNames = uniqueNames(order);
    uniqueHashes = uniqueHashes(order);

    nEvents = numel(uniqueFiles);
    inventoryRows = cell(nEvents, 9);
    qEvent = strings(nEvents*numel(labels),1);
    qChannel = strings(nEvents*numel(labels),1);
    qFinite = nan(nEvents*numel(labels),1);
    qGap = nan(nEvents*numel(labels),1);
    qUsable = false(nEvents*numel(labels),1);
    qRange = nan(nEvents*numel(labels),1);
    qrow = 0;

    for e = 1:nEvents
        raw = load(uniqueFiles(e));
        if ~isfield(raw, 'predat_sg')
            error('KW51:MissingPredatSG', '%s has no predat_sg.', uniqueFiles(e));
        end
        sg = raw.predat_sg;
        Xraw = double(sg.tdata);
        rawLabels = kw51_normaliseLabels(sg.labels);
        fs = double(sg.fs);
        nOriginal = size(Xraw, 1);
        targetLength = ceil(nOriginal * targetFs / fs);
        strain = nan(numel(labels), targetLength);
        good = false(numel(labels),1);

        for j = 1:numel(labels)
            qrow = qrow + 1;
            hit = find(rawLabels == labels(j), 1);
            if isempty(hit)
                x = nan(nOriginal,1);
            else
                x = Xraw(:,hit);
            end
            finite = isfinite(x);
            finiteFraction = mean(finite);
            maxGap = kw51_maximumMissingGap(finite) / fs;
            values = x(finite);
            if isempty(values)
                amplitude = NaN;
            else
                amplitude = max(values) - min(values);
            end
            usable = finiteFraction >= 0.99 && maxGap <= 0.1 && ...
                ~isempty(values) && amplitude > 1e-8;
            if usable
                xp = kw51_constantEndInterpolation(x, finite);
                denominators=1:10000; numerators=round((targetFs/fs)*denominators);
                [~,best]=min(abs(numerators./denominators-targetFs/fs));
                p=numerators(best); q=denominators(best); common=gcd(p,q); p=p/common; q=q/common;
                xp = resample(xp, p, q);
                if numel(xp) < targetLength
                    xp(end+1:targetLength) = xp(end);
                else
                    xp = xp(1:targetLength);
                end
                strain(j,:) = xp(:).';
            end
            good(j) = usable;
            qEvent(qrow) = uniqueNames(e);
            qChannel(qrow) = labels(j);
            qFinite(qrow) = finiteFraction;
            qGap(qrow) = maxGap;
            qUsable(qrow) = usable;
            qRange(qrow) = amplitude * 1e6;
        end

        time = (0:targetLength-1) / targetFs;
        for j = find(good).'
            edge = min(100, targetLength);
            firstMedian = median(strain(j,1:edge), 'omitnan');
            finalMedian = median(strain(j,end-edge+1:end), 'omitnan');
            if targetLength > 1
                baseline = firstMedian + (finalMedian-firstMedian) * time/time(end);
            else
                baseline = firstMedian;
            end
            strain(j,:) = strain(j,:) - baseline;
        end

        temperature = NaN;
        if isfield(raw, 'predat_env') && isstruct(raw.predat_env)
            envLabels = kw51_normaliseLabels(raw.predat_env.labels);
            envData = double(raw.predat_env.data);
            hit = find(envLabels == "tBD31A", 1);
            if ~isempty(hit)
                temperature = envData(hit);
            end
        end
        fs50 = targetFs;
        event = erase(uniqueNames(e), '.mat');
        save(fullfile(recordsDir, char(event + ".mat")), 'strain', 'time', ...
            'labels', 'good', 'temperature', 'fs50');

        duration = (nOriginal-1)/fs;
        inventoryRows(e,:) = {uniqueNames(e), uniqueTimes(e), uniqueHashes(e), ...
            nOriginal, fs, duration, temperature, sum(good), string(uniqueFiles(e))};
        if verbose && mod(e,10) == 0
            fprintf('  read %d/%d unique records\n', e, nEvents);
        end
    end

    eventInventory = cell2table(inventoryRows, 'VariableNames', ...
        {'event','datetime','sha256','n_samples_original','fs_original', ...
        'duration_s','temperature_C','n_usable','source_path'});
    quality = table(qEvent, qChannel, qFinite, qGap, qUsable, qRange, ...
        'VariableNames', {'event','channel','finite_fraction','max_gap_s', ...
        'usable','range_microstrain'});

    availability = table(labels.', zeros(numel(labels),1), zeros(numel(labels),1), ...
        nan(numel(labels),1), nan(numel(labels),1), ...
        'VariableNames', {'channel','n_events','n_usable','min_finite','max_range_microstrain'});
    for j = 1:numel(labels)
        ix = quality.channel == labels(j);
        availability.n_events(j) = sum(ix);
        availability.n_usable(j) = sum(quality.usable(ix));
        availability.min_finite(j) = min(quality.finite_fraction(ix));
        availability.max_range_microstrain(j) = max(quality.range_microstrain(ix));
    end

    if nEvents ~= 62
        error('KW51:EventCount', ['This benchmark uses a fixed chronological ' ...
            '32/14/16 split and expects 62 unique events; found %d.'], nEvents);
    end

    [bLow, aLow] = butter(4, 2/(targetFs/2), 'low');
    Y = nan(nEvents, numel(active), nPhase);
    temperatures = nan(nEvents,1);
    groups = zeros(nEvents,1);
    directions = zeros(nEvents,1);
    ratios = nan(nEvents,1);
    windowStart = nan(nEvents,1);
    windowEnd = nan(nEvents,1);
    eventStems = erase(uniqueNames, '.mat');

    for e = 1:nEvents
        rec = load(fullfile(recordsDir, char(eventStems(e) + ".mat")));
        X = rec.strain(active,:);
        if any(~isfinite(X), 'all')
            error('KW51:IncompleteActiveChannels', ...
                'An active strain channel is incomplete in %s.', eventStems(e));
        end
        X = filtfilt(bLow, aLow, X.').';
        energy = sqrt(sum(X.^2, 1));
        ix = find(energy > 0.04*max(energy));
        if isempty(ix)
            error('KW51:NoPassageWindow', 'No active passage found in %s.', eventStems(e));
        end
        lo = max(1, ix(1)-15);
        hi = min(size(X,2), ix(end)+15);
        phaseGrid = linspace(lo, hi, nPhase);
        Ya = zeros(numel(active), nPhase);
        for j = 1:numel(active)
            Ya(j,:) = interp1(1:size(X,2), X(j,:), phaseGrid, 'linear');
        end
        axial = (Ya(8,:) + Ya(9,:))/2;
        firstLobe = find(abs(axial) > 0.3*max(abs(axial)), 1);
        direction = 1;
        if axial(firstLobe) <= 0
            direction = -1;
            Ya = Ya(:,end:-1:1);
        end
        ratio = rms(Ya(1,:)) / max(rms(Ya(6,:)), 1e-15);
        group = double(ratio >= 0.45);

        Y(e,:,:) = Ya;
        temperatures(e) = rec.temperature;
        groups(e) = group;
        directions(e) = direction;
        ratios(e) = ratio;
        windowStart(e) = (lo-1)/targetFs;
        windowEnd(e) = (hi-1)/targetFs;
    end

    aligned = struct;
    aligned.Y = Y;
    aligned.groups = groups;
    aligned.temperature = temperatures;
    aligned.events = eventStems;
    aligned.labels = labels(active);
    aligned.fitIndices = (1:32).';
    aligned.calibrationIndices = (33:46).';
    aligned.testIndices = (47:62).';
    aligned.fs = targetFs;
    aligned.nPhase = nPhase;

    alignedTable = table(eventStems, uniqueTimes, temperatures, windowStart, windowEnd, ...
        directions, groups, ratios, 'VariableNames', {'event','datetime','temperature_C', ...
        'window_start_s','window_end_s','passage_direction_surrogate', ...
        'response_class','class_ratio'});
end
