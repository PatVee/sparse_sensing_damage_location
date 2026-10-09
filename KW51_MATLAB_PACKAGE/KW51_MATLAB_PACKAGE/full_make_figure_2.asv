function report = full_make_figure(outputDir, dpi, varargin)
%FULL_MAKE_FIGURE Rebuild KW51 figures from an UNZIPPED output folder.
%
%   full_make_figure('C:\project\output')
%   report = full_make_figure('C:\project\output', 300);
%   full_make_figure('C:\project\output', 300, ...
%       'PerLocation', false, 'Spatial', false);   % six main figures only
%
% Input folder: output/tables/*.csv and output/trials/*.mat.
% The folder may be called output, outputs, or any explicitly supplied name.
% This function never extracts archives, reruns FE analysis, calibrates a
% classifier, re-ranks features, or substitutes selected-feature counts into
% the common-family plot. Common-family bars match all four keys explicitly:
% feature, zone_name, severity_pct and n_passages.
%
% Name-value options:
%   'FigureDir'    : destination (default output/figures/full_make_figure_matlab)
%   'PerLocation' : true -> all saved zone-specific probability figures
%   'Spatial'     : true -> selected AND common-family spatial predictions
%   'CommonFeature': 'area_gram_dct8' (one exact ID, shared by every zone)
%   'LossPct'     : 10 (saved stiffness-loss level for stacked bar figures)
%   'LocationOrder': optional string array of retained zone names
%   'Visible'     : 'off' (default) or 'on'; exported figures are then closed
%
% Default for the supplied run: 6 main + 68 per-zone + 20 spatial figures.
% Every figure is exported as PDF and PNG. Source-row CSVs, LaTeX tables and
% run_manifest.json are written beneath FigureDir. Original files are not
% overwritten. Requires MATLAB tables, string arrays and tiledlayout; no
% specialised toolbox is used. This file is self-contained.

if nargin < 1, outputDir = ''; end
if nargin < 2 || isempty(dpi), dpi = 300; end
validateattributes(dpi, {'numeric'}, {'scalar','real','finite','positive','integer'});
p = inputParser;
p.addParameter('FigureDir', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
p.addParameter('PerLocation', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Spatial', true, @(x) islogical(x) && isscalar(x));
p.addParameter('CommonFeature', 'area_gram_dct8', @(x) ischar(x) || (isstring(x) && isscalar(x)));
p.addParameter('LossPct', 10, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('LocationOrder', strings(0,1), @(x) isstring(x) || iscellstr(x));
p.addParameter('Visible', 'off', @(x) any(strcmpi(string(x), ["off","on"])));
p.parse(varargin{:}); opt = p.Results;
opt.CommonFeature = string(opt.CommonFeature);
root = resolve_output_dir(outputDir);
if strlength(string(opt.FigureDir)) == 0
    figureDir = fullfile(root, 'figures', 'full_make_figure_matlab');
else
    figureDir = char(opt.FigureDir);
end
fprintf('full_make_figure 2026-09-24-output-folder-v1\nInput folder: %s\nFigure folder: %s\n', root, figureDir);
D = load_inputs(root, opt);
if opt.Spatial
    [spatialItems, spatialAudit, trialFiles] = prepare_spatial(D);
    D.sourceFiles = [D.sourceFiles; trialFiles];
else
    spatialItems = struct([]); spatialAudit = table();
end
fprintf('Loaded %d result rows. Common feature: %s\n', height(D.all), char(D.commonFeature));
ensure_dir(figureDir);
outputs = cell(0,1);
% Common-family figure is drawn first. A missing family/row is an error.
[files, commonLow] = plot_low_severity_trials(D, figureDir, dpi, opt, 'common');
outputs = [outputs; files];
[files, selectedLow] = plot_low_severity_trials(D, figureDir, dpi, opt, 'selected');
outputs = [outputs; files];
outputs = [outputs; plot_overview(D, figureDir, dpi, opt, false)];
outputs = [outputs; plot_overview(D, figureDir, dpi, opt, true)];
outputs = [outputs; plot_joint_curves(D, figureDir, dpi, opt, 'selected')];
outputs = [outputs; plot_joint_curves(D, figureDir, dpi, opt, 'common')];
write_tables(D, figureDir, selectedLow, commonLow);
if opt.PerLocation
    outputs = [outputs; plot_per_location(D, figureDir, dpi, opt)];
end
if opt.Spatial
    outputs = [outputs; plot_spatial(D, spatialItems, figureDir, dpi, opt)];
    write_plot_data(spatialAudit, fullfile(figureDir, 'plot_data', 'mat_csv_count_checks.csv'));
end
report = struct();
report.version = '2026-09-24-output-folder-v1';
report.implementation = 'MATLAB';
report.generated_utc = char(datetime('now','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss'Z'"));
report.input_folder = root; report.figure_folder = figureDir; report.dpi = dpi;
report.common_feature = char(D.commonFeature); report.low_loss_pct = D.lossPct;
report.locations = cellstr(D.order); report.passage_windows = D.windows;
report.severities_pct = D.severities; report.all_feature_result_rows = height(D.all);
report.common_plot_rows = height(commonLow); report.raw_mat_csv_checks = height(spatialAudit);
report.figure_count = numel(outputs)/2; report.png_and_pdf_files = outputs;
report.source_files = D.sourceFiles;
report.interpretation = ['Saved exploratory semi-synthetic outcomes. A common family does not ' ...
    'retrospectively make feature/zone selection independently validated.'];
report.arch_guide_note = ['01b uses visual interpolation through diagonal upper ends, ' ...
    'not saved arch-element connectivity.'];
try
    json = jsonencode(report, 'PrettyPrint', true);
catch
    json = jsonencode(report);
end
write_lines(fullfile(figureDir, 'run_manifest.json'), {json});
fprintf('Completed: %d figures (PNG + PDF), source-row CSVs and LaTeX tables.\n', report.figure_count);
fprintf('Common-family source counts at %g%% loss:\n', D.lossPct);
disp(sortrows(commonLow(:, {'feature','zone_name','n_passages', ...
    'detected_and_correctly_localised_trials','neither_trials','display_location_number'}), ...
    {'n_passages','display_location_number'}));
end


function root = resolve_output_dir(requested)
if strlength(string(requested)) > 0
    root = char(requested);
    if isfile(root)
        error('KW51:OutputFolder', 'Expected an unzipped output folder, not a file: %s', root);
    end
    [parent, last] = fileparts(root);
    if strcmpi(last, 'tables') && isfile(fullfile(root, 'all_location_feature_trials.csv'))
        root = parent;
    end
    if ~isfile(fullfile(root, 'tables', 'all_location_feature_trials.csv'))
        error('KW51:OutputFolder', ['Cannot find %s. Pass the output folder itself, ' ...
            'which contains tables and trials.'], fullfile(root, 'tables', 'all_location_feature_trials.csv'));
    end
    return;
end
here = fileparts(mfilename('fullpath'));
bases = {pwd, here, fileparts(here)};
for i = 1:numel(bases)
    candidates = {fullfile(bases{i}, 'output'), fullfile(bases{i}, 'outputs'), bases{i}};
    for j = 1:numel(candidates)
        if isfile(fullfile(candidates{j}, 'tables', 'all_location_feature_trials.csv'))
            root = candidates{j}; return;
        end
    end
end
error('KW51:OutputFolder', 'No output folder found. Supply full_make_figure(''C:\\project\\output'').');
end


function D = load_inputs(root, opt)
keys = {'regions','zones','sensors','retained','selected','all','catalogue'};
files = {'candidate_regions.csv','zone_definitions.csv','strain_sensor_locations.csv', ...
    'retained_best_locations.csv','retained_location_trials.csv','all_location_feature_trials.csv', ...
    'feature_catalogue.csv'};
if opt.PerLocation
    keys{end+1} = 'best'; files{end+1} = 'best_feature_trials_for_every_location.csv';
end
D = struct(); D.sourceFiles = cell(0,1);
for i = 1:numel(keys)
    path = fullfile(root, 'tables', files{i});
    if ~isfile(path), error('KW51:MissingTable', 'Required result table is missing: %s', path); end
    T = readtable(path, 'TextType', 'string');
    for name = ["feature","feature_label","zone_name","family","gauge_name","element_ids","memberName"]
        if ismember(char(name), T.Properties.VariableNames)
            T.(char(name)) = strtrim(string(T.(char(name))));
        end
    end
    D.(keys{i}) = T; D.sourceFiles{end+1,1} = path;
end
validate_trial_table(D.all, 'all_location_feature_trials.csv');
validate_trial_table(D.selected, 'retained_location_trials.csv');
if opt.PerLocation, validate_trial_table(D.best, 'best_feature_trials_for_every_location.csv'); end
require_columns(D.regions, {'elementID','family','memberName','x1_m','y1_m','z1_m','x2_m','y2_m','z2_m'}, 'candidate_regions');
require_columns(D.zones, {'zone_index','zone_name','x_m','y_m','z_m','element_ids','physical_length_m','n_elements'}, 'zone_definitions');
require_columns(D.sensors, {'gauge_name','x_m','y_m','z_m'}, 'strain_sensor_locations');
require_columns(D.catalogue, {'feature','feature_label'}, 'feature_catalogue');
checks = {'zones','zone_name'; 'zones','zone_index'; 'regions','elementID'; ...
    'sensors','gauge_name'; 'retained','zone_name'};
for i = 1:size(checks,1)
    values = D.(checks{i,1}).(checks{i,2});
    if numel(unique(values)) ~= numel(values)
        error('KW51:DuplicateID', '%s: duplicate %s.', checks{i,1}, checks{i,2});
    end
end
if isempty(D.retained), error('KW51:NoRetainedZones', 'No retained zones are present.'); end
preferred = ["diagonal_A_11","diagonal_A_14","main_girder_A_X012", ...
    "main_girder_C_X084","main_girder_C_X108"];
savedNames = string(D.retained.zone_name(:)).';
if isempty(opt.LocationOrder)
    order = preferred(ismember(preferred, savedNames));
    order = [order, savedNames(~ismember(savedNames, order))];
else
    order = string(opt.LocationOrder(:)).';
    if isempty(order) || numel(unique(order)) ~= numel(order) || ~all(ismember(order, savedNames))
        error('KW51:LocationOrder', 'LocationOrder must contain unique names from retained_best_locations.csv.');
    end
end
if ~all(ismember(order, D.zones.zone_name)), error('KW51:Zones', 'A retained zone has no geometry.'); end
D.root = root; D.order = order; D.commonFeature = opt.CommonFeature;
D.lossPct = opt.LossPct; D.windows = sort(unique(double(D.selected.n_passages))).';
if ~isequal(D.windows, [16 32])
    error('KW51:Windows', 'These manuscript layouts require passage windows [16 32].');
end
D.severities = sort(unique(double(D.all.severity_pct))).';
if ~any(abs(D.severities-D.lossPct) < 1e-8)
    error('KW51:LossLevel', 'Loss %g%% is not saved.', D.lossPct);
end
if ~ismember(D.commonFeature, D.catalogue.feature)
    error('KW51:Feature', 'Unknown common family: %s', char(D.commonFeature));
end
D.labels = strings(size(order));
for i = 1:numel(order), D.labels(i) = zone_label(get_zone(D, order(i))); end
if ismember('retention_threshold', D.retained.Properties.VariableNames)
    threshold = unique(double(D.retained.retention_threshold));
else
    threshold = 0.5;
end
if numel(threshold) ~= 1, error('KW51:Threshold', 'Expected one retention threshold.'); end
D.retentionThreshold = threshold;
% Verify saved selections without choosing new features.
selectionKeys = {'selected'};
if opt.PerLocation, selectionKeys{end+1} = 'best'; end
countFields = [{'total_trials','detected_trials','correctly_localised_trials'}, count_columns()];
for k = 1:numel(selectionKeys)
    T = D.(selectionKeys{k}); pairs = unique(T(:, {'zone_name','n_passages'}), 'rows');
    for j = 1:height(pairs)
        zone = string(pairs.zone_name(j)); n = double(pairs.n_passages(j));
        feature = feature_for(T, zone, n);
        q = T(T.zone_name == zone & T.n_passages == n, :);
        if ~isequal(sort(double(q.severity_pct(:))), D.severities(:))
            error('KW51:MissingSeverity', '%s: missing or duplicate severities for %s, %d passages.', selectionKeys{k}, char(zone), n);
        end
        if strcmp(selectionKeys{k}, 'selected')
            r = D.retained(D.retained.zone_name == zone, :);
            field = sprintf('best_feature_%dpass', n);
            if height(r) ~= 1 || ~ismember(field, r.Properties.VariableNames) || string(r.(field)) ~= feature
                error('KW51:SelectionMismatch', 'Saved feature summary disagrees for %s, %d passages.', char(zone), n);
            end
        end
        for s = 1:height(q)
            original = match_row(D.all, feature, zone, n, double(q.severity_pct(s)));
            if ~isequal(double(q{s,countFields}), double(original{1,countFields}))
                error('KW51:CountMismatch', 'Selected/full-result counts disagree for %s, %d passages.', char(zone), n);
            end
        end
    end
end
for i = 1:numel(D.order)
    zone_geometry(D, D.order(i));
    for n = D.windows
        for loss = D.severities
            match_row(D.all, D.commonFeature, D.order(i), n, loss);
        end
    end
end
end


function require_columns(T, columns, context)
missing = columns(~ismember(columns, T.Properties.VariableNames));
if ~isempty(missing)
    error('KW51:MissingColumns', '%s: missing columns: %s', context, strjoin(missing, ', '));
end
end


function cols = count_columns()
cols = {'detected_and_correctly_localised_trials','detected_wrong_location_trials', ...
    'correctly_localised_without_detection_trials'};
end


function validate_trial_table(T, context)
keys = {'feature','zone_name','n_passages','severity_pct'};
cols = [{'total_trials','detected_trials','correctly_localised_trials'}, count_columns()];
require_columns(T, [keys, {'zone_index'}, cols], context);
if isempty(T), error('KW51:EmptyTable', '%s is empty.', context); end
if any(any(ismissing(T(:, keys))))
    error('KW51:MissingKey', '%s has missing key values.', context);
end
if height(unique(T(:, keys), 'rows')) ~= height(T)
    error('KW51:DuplicateRows', '%s has duplicate feature/location/window/loss rows.', context);
end
a = double(T{:,cols});
if any(~isfinite(a(:))) || any(a(:) < 0) || any(abs(a(:)-round(a(:))) > 1e-8)
    error('KW51:Counts', '%s: counts must be finite non-negative integers.', context);
end
N = double(T.total_trials); counts = double(T{:,count_columns()});
if any(N <= 0) || any(sum(counts,2) > N)
    error('KW51:Counts', '%s: invalid denominator or category counts exceed total_trials.', context);
end
if ~isequal(double(T.detected_trials), counts(:,1)+counts(:,2)) || ...
        ~isequal(double(T.correctly_localised_trials), counts(:,1)+counts(:,3))
    error('KW51:Counts', '%s: detection/localisation totals disagree with category counts.', context);
end
stems = {'detected','correctly_localised','detected_and_correctly_localised'};
rates = {'detection_rate','localisation_rate','detected_and_correctly_localised_rate'};
for k = 1:3
    lowName = [stems{k}, '_wilson95_low']; highName = [stems{k}, '_wilson95_high'];
    require_columns(T, {rates{k},lowName,highName}, context);
    expected = double(T.([stems{k}, '_trials']))./N;
    if any(abs(double(T.(rates{k}))-expected) > 1e-8)
        error('KW51:Rates', '%s: rates do not match counts for %s.', context, stems{k});
    end
    lo = double(T.(lowName)); hi = double(T.(highName));
    if any(~isfinite(lo)) || any(~isfinite(hi)) || any(lo < -1e-10) || any(hi > 1+1e-10) || any(lo > hi)
        error('KW51:Intervals', '%s: invalid interval bounds.', context);
    end
end
if ismember('loss', T.Properties.VariableNames) && any(abs(100*double(T.loss)-double(T.severity_pct)) > 1e-8)
    error('KW51:Loss', '%s: loss and severity_pct disagree.', context);
end
end


function row = match_row(T, feature, zone, n, loss)
row = T(string(T.feature) == string(feature) & string(T.zone_name) == string(zone) & ...
    T.n_passages == n & abs(double(T.severity_pct)-loss) < 1e-8, :);
if height(row) ~= 1
    error('KW51:UniqueRow', 'Expected one row for %s / %s / %d passages / %g%% loss; found %d.', ...
        char(feature), char(zone), n, loss, height(row));
end
end


function feature = feature_for(T, zone, passages)
features = unique(string(T.feature(string(T.zone_name) == string(zone) & T.n_passages == passages)));
if numel(features) ~= 1
    error('KW51:SelectedFeature', 'Expected one saved selected family for %s, %d passages.', char(zone), passages);
end
feature = features(1);
end


function label = feature_label(D, feature)
q = D.catalogue(D.catalogue.feature == string(feature), :);
if height(q) ~= 1, error('KW51:FeatureCatalogue', 'Expected one catalogue row for %s.', char(feature)); end
label = string(q.feature_label);
end


function row = get_zone(D, name)
row = D.zones(string(D.zones.zone_name) == string(name), :);
if height(row) ~= 1, error('KW51:ZoneGeometry', 'Expected one definition for %s.', char(name)); end
end


function label = zone_label(row)
name = char(string(row.zone_name));
tok = regexp(name, '^diagonal_([AC])_(\d+)$', 'tokens', 'once');
if ~isempty(tok)
    label = string(sprintf('%s%02d diagonal', tok{1}, str2double(tok{2}))); return;
end
tok = regexp(name, '^main_girder_([AC])_X\d+$', 'tokens', 'once');
if ~isempty(tok)
    label = string(sprintf('Main girder %s, x=%.1f m', tok{1}, row.x_m)); return;
end
label = strrep(string(name), '_', ' ');
end


function q = zone_geometry(D, name)
row = get_zone(D, name);
parts = split(string(row.element_ids), ';'); parts = parts(strlength(strtrim(parts)) > 0);
ids = str2double(parts);
if any(~isfinite(ids)) || numel(unique(ids)) ~= numel(ids)
    error('KW51:ElementIDs', 'Invalid/repeated element IDs in %s.', char(name));
end
q = D.regions(ismember(double(D.regions.elementID), ids), :);
if height(q) ~= numel(ids), error('KW51:ElementGeometry', 'Missing element geometry for %s.', char(name)); end
q = sortrows(q, 'elementID');
end


function fig = new_figure(w, h, opt)
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 w h], ...
    'PaperUnits','inches','PaperPosition',[0 0 w h],'PaperSize',[w h], ...
    'Visible',char(opt.Visible));
end


function set_axes_style(ax)
set(ax, 'FontName','Times New Roman','FontSize',8.5,'LineWidth',0.8, ...
    'Color','w','GridColor',repmat(0.82,1,3),'GridAlpha',0.72,'MinorGridAlpha',0.45);
end


function files = safe_export(fig, directory, stem, dpi)
ensure_dir(directory);
pdf = fullfile(directory, [stem,'.pdf']); png = fullfile(directory, [stem,'.png']);
tmpPdf = fullfile(directory, [stem,'.writing.pdf']); tmpPng = fullfile(directory, [stem,'.writing.png']);
if isfile(tmpPdf), delete(tmpPdf); end
if isfile(tmpPng), delete(tmpPng); end
drawnow;
if exist('exportgraphics','file') == 2
    exportgraphics(fig, tmpPdf, 'ContentType','vector','BackgroundColor','white');
    exportgraphics(fig, tmpPng, 'Resolution',dpi,'BackgroundColor','white');
else
    print(fig, tmpPdf, '-dpdf', '-painters');
    print(fig, tmpPng, '-dpng', sprintf('-r%d',dpi));
end
movefile(tmpPdf, pdf, 'f'); movefile(tmpPng, png, 'f');
files = {pdf;png};
end


function h = draw_segments(ax, frame, gray, width)
if isempty(frame), h = plot3(ax,nan,nan,nan,'Color',repmat(gray,1,3)); return; end
n = height(frame);
x = [frame.x1_m,frame.x2_m,nan(n,1)].';
y = [frame.y1_m,frame.y2_m,nan(n,1)].';
z = [frame.z1_m,frame.z2_m,nan(n,1)].';
h = plot3(ax,x(:),y(:),z(:),'Color',repmat(gray,1,3),'LineWidth',width);
end


function draw_bridge(ax, D, light)
names = ["deck_longitudinal","deck_transverse","main_girder","diagonal"];
grays = [0.88 0.84 0.58 0.67]; widths = [0.34 0.40 0.72 0.78];
for i = 1:numel(names)
    q = D.regions(D.regions.family == names(i), :);
    if light, gray = 0.85; width = 0.30; else, gray = grays(i); width = widths(i); end
    draw_segments(ax,q,gray,width);
end
end


function draw_sensors(ax, D, small)
xyz = double(D.sensors{:,{'x_m','y_m','z_m'}});
[positions,~,groups] = unique(xyz,'rows'); counts = accumarray(groups,1);
if small, base = 12; else, base = 34; end
scatter3(ax,positions(:,1),positions(:,2),positions(:,3),base+10*(counts-1), ...
    [0.75 0 0],'filled','MarkerEdgeColor','w','LineWidth',0.65);
end


function draw_arch_guides(ax, D)
% Visual guide only: the candidate CSV omits full arch beam connectivity.
d = D.regions(D.regions.family == "diagonal", :);
names = unique(d.memberName,'stable'); points = zeros(numel(names),3);
for i = 1:numel(names)
    m = d(d.memberName == names(i),:);
    ends = unique(round([m{:,{'x1_m','y1_m','z1_m'}};m{:,{'x2_m','y2_m','z2_m'}}],8),'rows');
    [~,k] = max(ends(:,3)); points(i,:) = ends(k,:);
end
sides = unique(round(points(:,2),3));
for i = 1:numel(sides)
    P = points(abs(points(:,2)-sides(i)) < 1e-3,:); P = sortrows(P,1);
    [~,ix] = unique(P(:,1),'stable'); P = P(ix,:);
    if size(P,1) < 2, continue; end
    xf = linspace(min(P(:,1)),max(P(:,1)),250);
    if size(P,1) >= 3, method = 'pchip'; else, method = 'linear'; end
    yf = interp1(P(:,1),P(:,2),xf,'linear'); zf = interp1(P(:,1),P(:,3),xf,method);
    plot3(ax,xf,yf,zf,'Color',repmat(0.38,1,3),'LineWidth',1.2);
end
end


function set_bridge_axes(ax, small, az, el)
xlim(ax,[-2 118.5]); ylim(ax,[-1 12.8]); zlim(ax,[0 20]);
xticks(ax,0:20:120); yticks(ax,[0 6 12]); zticks(ax,[0 5 10 15 20]);
if small
    xlabel(ax,'x (m)'); ylabel(ax,'y (m)'); zlabel(ax,'z (m)');
else
    xlabel(ax,'Longitudinal coordinate, x (m)'); ylabel(ax,'y (m)'); zlabel(ax,'Elevation, z (m)');
end
view(ax,az,el); pbaspect(ax,[5.4 1.25 1.55]); grid(ax,'on'); box(ax,'on'); set_axes_style(ax);
if small, ax.FontSize = 7; end
end


function files = plot_overview(D, figureDir, dpi, opt, arches)
fig = new_figure(12.0,5.2,opt); cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
ax = axes(fig); hold(ax,'on'); draw_bridge(ax,D,false);
if arches, draw_arch_guides(ax,D); end
offsets = [0 0 1.3;1 0 2.1;0 0 2.2;-4.5 0 3;0 0 2.2];
handles = gobjects(numel(D.order)+1,1); labels = strings(numel(handles),1);
for i = 1:numel(D.order)
    zone = zone_geometry(D,D.order(i)); draw_segments(ax,zone,0.02,4.4);
    ends = [zone{:,{'x1_m','y1_m','z1_m'}};zone{:,{'x2_m','y2_m','z2_m'}}];
    plot3(ax,ends(:,1),ends(:,2),ends(:,3),'LineStyle','none','Marker','o', ...
        'MarkerSize',4,'MarkerFaceColor','w','MarkerEdgeColor',repmat(0.02,1,3),'LineWidth',0.65);
    centre = mean(ends,1);
    if i <= size(offsets,1), off = offsets(i,:); else, off = [0 0 2.2]; end
    tag = centre+off; tag(3) = min(tag(3),19.3);
    plot3(ax,[centre(1),tag(1)],[centre(2),tag(2)],[centre(3),tag(3)-0.25], ...
        'Color',repmat(0.15,1,3),'LineWidth',0.7);
    scatter3(ax,tag(1),tag(2),tag(3),125,'w','filled','MarkerEdgeColor',repmat(0.05,1,3),'LineWidth',0.9);
    text(ax,tag(1),tag(2),tag(3),num2str(i),'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontWeight','bold','FontSize',8.5,'Color','k','Clipping','off');
    handles(i) = plot3(ax,nan,nan,nan,'-','Color',repmat(0.02,1,3),'LineWidth',3.2);
    labels(i) = string(sprintf('%d  %s',i,char(D.labels(i))));
end
draw_sensors(ax,D,false);
handles(end) = plot3(ax,nan,nan,nan,'o','LineStyle','none','MarkerSize',6.5, ...
    'MarkerFaceColor',[0.75 0 0],'MarkerEdgeColor','w');
labels(end) = "Measured strain-gauge locations";
if arches, set_bridge_axes(ax,false,-45,18); else, set_bridge_axes(ax,false,-62,20); end
title(ax,'Retained stiffness-loss zones and strain-gauge locations','FontSize',11);
lgd = legend(ax,handles,cellstr(labels),'Location','southoutside','Orientation','horizontal','Interpreter','none','Box','off');
if isprop(lgd,'NumColumns'), lgd.NumColumns = 3; end
if arches
    ax.Position = [0.055 0.28 0.89 0.61];
    lgd.Position = [0.08 0.075 0.84 0.13];
    annotation(fig,'textbox',[0.04 0.008 0.92 0.05], ...
        'String','Arch curves are visual guides through diagonal-to-arch connection coordinates, not additional FE elements.', ...
        'EdgeColor','none','HorizontalAlignment','center','FontName','Times New Roman','FontSize',7);
    stem = '01b_retained_locations_3d_arch_guides';
else
    stem = '01_retained_locations_3d';
end
files = safe_export(fig,figureDir,stem,dpi);
end


function files = plot_joint_curves(D, figureDir, dpi, opt, mode)
fig = new_figure(11.8,4.55,opt); cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
axesHandles = gobjects(2,1); legendHandles = gobjects(numel(D.order),1);
grays = [0.05 0.30 0.48 0.65 0.80]; markers = {'o','s','^','d','v'}; styles = {'-','--','-.',':','-'};
rows = cell(numel(D.order)*2,1); cursor = 0;
for p = 1:2
    ax = nexttile(layout,p); axesHandles(p) = ax; hold(ax,'on'); n = D.windows(p);
    for i = 1:numel(D.order)
        if strcmp(mode,'common'), feature = D.commonFeature; else, feature = feature_for(D.selected,D.order(i),n); end
        q = D.all(D.all.feature == feature & D.all.zone_name == D.order(i) & D.all.n_passages == n,:);
        q = sortrows(q,'severity_pct');
        y = 100*double(q.detected_and_correctly_localised_rate);
        lo = 100*double(q.detected_and_correctly_localised_wilson95_low);
        hi = 100*double(q.detected_and_correctly_localised_wilson95_high);
        k = mod(i-1,5)+1;
        h = errorbar(ax,double(q.severity_pct),y,max(y-lo,0),max(hi-y,0),'Color',repmat(grays(k),1,3), ...
            'Marker',markers{k},'LineStyle',styles{k},'MarkerSize',4.7,'MarkerFaceColor','w','LineWidth',1.45,'CapSize',4);
        if p == 1, legendHandles(i) = h; end
        q.display_location_number = repmat(i,height(q),1); cursor = cursor+1; rows{cursor} = q;
    end
    threshold = 100*D.retentionThreshold;
    yline(ax,threshold,':','Color',repmat(0.62,1,3),'LineWidth',0.8);
    if strcmp(mode,'selected'), thresholdText = 'retention threshold'; else, thresholdText = sprintf('%g%% reference',threshold); end
    text(ax,max(D.severities)-0.5,threshold+1.6,thresholdText,'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom','Color',repmat(0.42,1,3),'FontSize',8);
    title(ax,sprintf('%d accumulated passages',n)); xlabel(ax,'Local stiffness loss (%)');
    xticks(ax,D.severities); ylim(ax,[0 103]); yticks(ax,0:20:100); grid(ax,'on'); box(ax,'on'); set_axes_style(ax);
end
ylabel(axesHandles(1),'Detected and correctly localised trials (%)'); linkaxes(axesHandles,'xy');
if strcmp(mode,'selected')
    title(layout,'Selected feature family by location and passage window','FontSize',11);
    stem = '02_accepted_localisation_probability';
else
    title(layout,['Common ',char(feature_label(D,D.commonFeature)),' feature family'],'FontSize',11,'Interpreter','none');
    stem = ['05_common_',char(D.commonFeature),'_localisation_probability'];
end
lgd = legend(axesHandles(1),legendHandles,cellstr(D.labels),'Orientation','horizontal','Interpreter','none','Box','off');
if isprop(lgd,'NumColumns'), lgd.NumColumns = 3; end
lgd.Layout.Tile = 'south';
write_plot_data(vertcat(rows{:}),fullfile(figureDir,'plot_data',[stem,'.csv']));
files = safe_export(fig,figureDir,stem,dpi);
end


function [files, plotData] = plot_low_severity_trials(D, figureDir, dpi, opt, mode)
% Every bar uses FEATURE + ZONE + LOSS + WINDOW, including in common mode.
allTrials = D.all; displayOrder = fliplr(D.order); nZones = numel(displayOrder);
valuesByWindow = zeros(nZones,4,2); rows = cell(nZones*2,1); cursor = 0;
for p = 1:2
    passages = D.windows(p);
    for i = 1:nZones
        if strcmp(mode,'common')
            feature = D.commonFeature; % default: exactly "area_gram_dct8"
        else
            feature = feature_for(D.selected,displayOrder(i),passages);
        end
        % IMPORTANT: explicit family selection is inside the plotting routine.
        row = allTrials( ...
            string(allTrials.feature) == feature & ...
            string(allTrials.zone_name) == displayOrder(i) & ...
            allTrials.severity_pct == D.lossPct & ...
            allTrials.n_passages == passages, :);
        if height(row) ~= 1
            error('KW51:BarSource', 'Expected one bar source for %s / %s / %d passages / %g%% loss; found %d.', ...
                char(feature),char(displayOrder(i)),passages,D.lossPct,height(row));
        end
        counts = double(row{1,count_columns()});
        valuesByWindow(i,:,p) = [counts,double(row.total_trials)-sum(counts)];
        row.neither_trials = valuesByWindow(i,4,p);
        row.display_location_number = find(D.order == displayOrder(i),1);
        row.plot_mode = string(mode); row.source_table = "tables/all_location_feature_trials.csv";
        cursor = cursor+1; rows{cursor} = row;
    end
end
plotData = vertcat(rows{:}); totals = unique(double(plotData.total_trials));
if numel(totals) ~= 1, error('KW51:BarTotals','All stacked bars must share one trial denominator.'); end
N = totals(1);
if strcmp(mode,'common') && (numel(unique(plotData.feature)) ~= 1 || any(plotData.feature ~= D.commonFeature))
    error('KW51:CommonFeature','The common-family plot contains a different feature.');
end
fig = new_figure(11.8,5.1,opt); cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
axesHandles = gobjects(2,1); legendHandles = gobjects(4,1);
categoryLabels = {'Detected + correct','Detected, wrong location','Correct location, not detected','Neither'};
categoryGray = [0.12 0.46 0.76 1.00];
for p = 1:2
    ax = nexttile(layout,p); axesHandles(p) = ax; hold(ax,'on'); passages = D.windows(p);
    values = valuesByWindow(:,:,p);
    bars = barh(ax,1:nZones,values,'stacked','BarWidth',0.64);
    for k = 1:4
        bars(k).FaceColor = repmat(categoryGray(k),1,3);
        bars(k).EdgeColor = repmat(0.08,1,3); bars(k).LineWidth = 0.55;
        if p == 1, legendHandles(k) = bars(k); end
    end
    cumulative = zeros(nZones,1);
    for k = 1:4
        for i = 1:nZones
            count = values(i,k);
            if count >= 18
                textColour = repmat(0.05,1,3); if k == 1, textColour = [1 1 1]; end
                text(ax,cumulative(i)+count/2,i,num2str(count),'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle','FontSize',12,'Color',textColour);
            end
        end
        cumulative = cumulative+values(:,k);
    end
    title(ax,sprintf('%d accumulated passages',passages)); xlabel(ax,sprintf('Trials (out of %d)',N));
    xlim(ax,[0 N]); if N == 300, step = 50; else, step = max(1,ceil(N/6)); end
    xticks(ax,0:step:N); yticks(ax,1:nZones); ylim(ax,[0.5,nZones+0.5]);
    if p == 1
        labels = strings(nZones,1);
        for i = 1:nZones
            originalIndex = find(D.order == displayOrder(i),1);
            labels(i) = string(sprintf('%d  %s',originalIndex,char(D.labels(originalIndex))));
        end
        yticklabels(ax,cellstr(labels));
    else
        yticklabels(ax,repmat({''},nZones,1));
    end
    grid(ax,'on'); ax.YGrid = 'off'; box(ax,'on'); set_axes_style(ax);
    set(ax, 'FontSize', 14);
end
linkaxes(axesHandles,'xy');
if strcmp(mode,'common')
    title(layout,sprintf('Prediction using common %s feature family (%g%% loss)', ...
        char(feature_label(D,D.commonFeature)),D.lossPct),'Interpreter','none');
    stem = sprintf('04_common_%s_%gpct_trial_outcomes',char(D.commonFeature),D.lossPct);
else
    title(layout,sprintf('Trial outcomes at %g%% local stiffness loss',D.lossPct));
    stem = '03_low_severity_trial_outcomes';
end
lgd = legend(axesHandles(1),legendHandles,categoryLabels,'Orientation','horizontal','Interpreter','none','Box','off');
if isprop(lgd,'NumColumns'), lgd.NumColumns = 4; end
lgd.Layout.Tile = 'south';
lgd.FontSize = 14;
write_plot_data(plotData,fullfile(figureDir,'plot_data',[stem,'.csv']));
files = safe_export(fig,figureDir,stem,dpi);
end


function outputs = plot_per_location(D, figureDir, dpi, opt)
outputs = cell(0,1); zones = sortrows(D.zones,'zone_index');
for k = 1:height(zones)
    row = zones(k,:);
    files = plot_one_location(D,row,figureDir,dpi,opt);
    outputs = [outputs;files]; %#ok<AGROW>
    if mod(k,10) == 0 || k == height(zones), fprintf('  per-location plots: %d/%d\n',k,height(zones)); end
end
end


function files = plot_one_location(D, zone, figureDir, dpi, opt)
fig = new_figure(11.8,4.55,opt); cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
axesHandles = gobjects(2,1); legendHandles = gobjects(3,1); rows = cell(2,1);
stems = {'detected','correctly_localised','detected_and_correctly_localised'};
rates = {'detection_rate','localisation_rate','detected_and_correctly_localised_rate'};
labels = {'Detection','Correct zone (ungated)','Detected + correct zone'};
styles = {'--',':','-'}; markers = {'s','^','o'}; grays = [0.10 0.55 0.30];
for p = 1:2
    n = D.windows(p); ax = nexttile(layout,p); axesHandles(p) = ax; hold(ax,'on');
    feature = feature_for(D.best,string(zone.zone_name),n);
    q = D.all(D.all.feature == feature & D.all.zone_name == string(zone.zone_name) & D.all.n_passages == n,:);
    q = sortrows(q,'severity_pct'); rows{p} = q;
    for k = 1:3
        y = 100*double(q.(rates{k})); lo = 100*double(q.([stems{k},'_wilson95_low']));
        hi = 100*double(q.([stems{k},'_wilson95_high']));
        h = errorbar(ax,double(q.severity_pct),y,max(y-lo,0),max(hi-y,0), ...
            'LineStyle',styles{k},'Marker',markers{k},'Color',repmat(grays(k),1,3), ...
            'MarkerFaceColor','w','MarkerSize',4,'LineWidth',1.25,'CapSize',4);
        if p == 1, legendHandles(k) = h; end
    end
    title(ax,sprintf('%d passages | %s',n,char(feature_label(D,feature))),'Interpreter','none');
    xlabel(ax,'Local stiffness loss (%)'); ylim(ax,[0 103]); xticks(ax,D.severities); yticks(ax,0:20:100);
    grid(ax,'on'); box(ax,'on'); set_axes_style(ax);
end
ylabel(axesHandles(1),'Trial fraction (%)'); linkaxes(axesHandles,'xy');
title(layout,[char(zone_label(zone)),': saved selected families and 95% Wilson intervals'],'FontSize',11,'Interpreter','none');
lgd = legend(axesHandles(1),legendHandles,labels,'Orientation','horizontal','Box','off','Interpreter','none');
if isprop(lgd,'NumColumns'), lgd.NumColumns = 3; end
lgd.Layout.Tile = 'south';
stem = sprintf('zone_%02d_probabilities',zone.zone_index);
write_plot_data(vertcat(rows{:}),fullfile(figureDir,'plot_data','per_location',[stem,'.csv']));
files = safe_export(fig,fullfile(figureDir,'per_location'),stem,dpi);
end


function d = load_trial_arrays(path)
if ~isfile(path)
    error('KW51:MissingTrials','Spatial maps need %s. Set ''Spatial'',false to omit spatial maps.',path);
end
names = {'predictedAtom','detected','truthZone','atomZone','atomLoss','n_passages'};
d = load(path,names{:});
for i = 1:numel(names)
    if ~isfield(d,names{i}), error('KW51:TrialVariable','Missing %s in %s.',names{i},path); end
end
d.atomZone = double(d.atomZone(:)); d.truthZone = double(d.truthZone(:)); d.atomLoss = double(d.atomLoss(:));
d.predictedAtom = double(d.predictedAtom);
if any(~ismember(d.detected(:),[0 1])), error('KW51:DetectionMask','detected must be binary in %s.',path); end
d.detected = logical(d.detected); nAtoms = numel(d.atomZone);
if ~ismatrix(d.predictedAtom) || size(d.predictedAtom,2) ~= nAtoms || ~isequal(size(d.predictedAtom),size(d.detected))
    error('KW51:TrialShape','Expected trials-by-atoms matrices with %d columns in %s.',nAtoms,path);
end
if numel(d.atomLoss) ~= nAtoms || ~isequal(d.truthZone,d.atomZone)
    error('KW51:AtomMap','Atom/truth maps disagree in %s.',path);
end
if any(~isfinite(d.predictedAtom(:))) || any(d.predictedAtom(:) < 0) || ...
        any(d.predictedAtom(:) >= nAtoms) || any(d.predictedAtom(:) ~= round(d.predictedAtom(:)))
    error('KW51:AtomIDs','predictedAtom must contain zero-based integer IDs in %s.',path);
end
end


function [histogram, audit] = spatial_counts(D, d, feature, zoneName, passages, lossPct)
zone = get_zone(D,zoneName); z = double(zone.zone_index);
column = find(d.truthZone == z & abs(100*d.atomLoss-lossPct) < 1e-8);
if numel(column) ~= 1, error('KW51:TruthAtom','Expected one truth atom for %s at %g%% loss.',char(zoneName),lossPct); end
if double(d.n_passages) ~= passages, error('KW51:TrialWindow','Saved n_passages differs from the requested window.'); end
% Saved predictedAtom is ZERO-BASED; add one only for MATLAB indexing.
prediction = d.atomZone(d.predictedAtom(:,column)+1);
detected = d.detected(:,column); correct = prediction == z; N = numel(detected);
counts = [sum(detected & correct),sum(detected & ~correct),sum(~detected & correct),sum(~detected & ~correct)];
row = match_row(D.all,feature,zoneName,passages,lossPct);
expected = double(row{1,count_columns()}); expected = [expected,double(row.total_trials)-sum(expected)];
if N ~= double(row.total_trials) || ~isequal(counts,expected)
    error('KW51:MatCsvMismatch','MAT/CSV mismatch for %s / %s / %d passages / %g%% loss.',char(feature),char(zoneName),passages,lossPct);
end
if ~all(ismember(prediction,double(D.zones.zone_index)))
    error('KW51:PredictedZone','A predicted zone has no saved geometry.');
end
histogram = D.zones(:,{'zone_index','zone_name','x_m','y_m','z_m'});
histogram.Properties.VariableNames(1:2) = {'predicted_zone_index','predicted_zone_name'};
histogram.detected_prediction_trials = zeros(height(histogram),1);
for i = 1:height(histogram)
    histogram.detected_prediction_trials(i) = sum(detected & prediction == double(histogram.predicted_zone_index(i)));
end
histogram.fraction_of_all_trials = histogram.detected_prediction_trials/N;
histogram.truth_zone_index = repmat(z,height(histogram),1);
histogram.truth_zone_name = repmat(string(zoneName),height(histogram),1);
histogram.feature = repmat(string(feature),height(histogram),1);
histogram.n_passages = repmat(passages,height(histogram),1);
histogram.severity_pct = repmat(lossPct,height(histogram),1);
histogram.total_trials = repmat(N,height(histogram),1);
histogram.is_correct_zone = histogram.predicted_zone_index == z;
audit = table(string(feature),string(zoneName),z,passages,lossPct,N,sum(detected),sum(correct), ...
    counts(1),counts(2),counts(3),counts(4),true,'VariableNames', ...
    {'feature','zone_name','zone_index','n_passages','severity_pct','total_trials','detected_trials', ...
    'correctly_localised_trials','detected_and_correctly_localised_trials','detected_wrong_location_trials', ...
    'correctly_localised_without_detection_trials','neither_trials','mat_csv_counts_match'});
end


function [items, audit, trialFiles] = prepare_spatial(D)
cache = containers.Map('KeyType','char','ValueType','any');
items = struct([]); auditRows = cell(0,1); cursor = 0;
for mode = ["selected","common"]
    for i = 1:numel(D.order)
        zoneName = D.order(i);
        for n = D.windows
            if mode == "common", feature = D.commonFeature; else, feature = feature_for(D.selected,zoneName,n); end
            path = fullfile(D.root,'trials',sprintf('%s_%dpass.mat',char(feature),n));
            if ~isKey(cache,path), cache(path) = load_trial_arrays(path); end
            d = cache(path); panels = struct([]);
            for j = 1:numel(D.severities)
                [h,a] = spatial_counts(D,d,feature,zoneName,n,D.severities(j));
                a.plot_mode = mode; a.source_mat = string(path);
                panels(j).histogram = h; panels(j).audit = a; %#ok<AGROW>
                auditRows{end+1,1} = a; %#ok<AGROW>
            end
            cursor = cursor+1;
            items(cursor).mode = mode; items(cursor).zone = zoneName; items(cursor).passages = n;
            items(cursor).feature = feature; items(cursor).panels = panels;
        end
    end
end
audit = vertcat(auditRows{:}); trialFiles = keys(cache).';
end


function outputs = plot_spatial(D, items, figureDir, dpi, opt)
outputs = cell(0,1);
for i = 1:numel(items)
    outputs = [outputs; plot_one_spatial(D,items(i),figureDir,dpi,opt)]; %#ok<AGROW>
    fprintf('  spatial plots: %d/%d\n',i,numel(items));
end
end


function files = plot_one_spatial(D, item, figureDir, dpi, opt)
nPanels = numel(item.panels); nRows = ceil((nPanels+1)/3);
fig = new_figure(15.6,4.2*nRows,opt); cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
layout = tiledlayout(fig,nRows,3,'TileSpacing','compact','Padding','compact');
zone = get_zone(D,item.zone); histograms = cell(nPanels,1);
for j = 1:nPanels
    ax = nexttile(layout,j); hold(ax,'on'); draw_bridge(ax,D,true);
    draw_segments(ax,zone_geometry(D,item.zone),0.02,2.8); draw_sensors(ax,D,true);
    h = item.panels(j).histogram; a = item.panels(j).audit;
    q = h(h.detected_prediction_trials > 0,:);
    if ~isempty(q)
        scatter3(ax,q.x_m,q.y_m,q.z_m,850*q.fraction_of_all_trials,'o', ...
            'MarkerFaceColor','none','MarkerEdgeColor',repmat(0.03,1,3),'LineWidth',1.0);
    end
    scatter3(ax,zone.x_m,zone.y_m,zone.z_m,45,'k','x','LineWidth',1.3);
    set_bridge_axes(ax,true,-55,18);
    title(ax,sprintf('%g%% loss | %d/%d detected',a.severity_pct,a.detected_trials,a.total_trials),'FontSize',9);
    histograms{j} = h;
end
notes = nexttile(layout,nPanels+1); axis(notes,'off');
noteText = sprintf(['Bold member / x: prescribed stiffness-loss zone\n' ...
    'Red dots: measured strain-gauge positions\nOpen circles: detected predicted zones\n' ...
    'Circle area: fraction of ALL trials\nCoordinates: saved FE zone centroids\n\n' ...
    'Each MAT histogram is checked against\nthe same feature/window/loss row in the CSV.\n' ...
    'Same-model, semi-synthetic assessment.']);
text(notes,0.03,0.90,noteText,'Units','normalized','VerticalAlignment','top','FontName','Times New Roman','FontSize',9,'Interpreter','none');
if item.mode == "common", subtitle = 'Common feature'; else, subtitle = 'Saved selected feature'; end
label = zone_label(zone);
title(layout,{sprintf('%s | %d accumulated passages',char(label),item.passages), ...
    sprintf('%s: %s',subtitle,char(feature_label(D,item.feature)))},'Interpreter','none','FontSize',12);
stem = sprintf('zone_%02d_%dpass_predictions',zone.zone_index,item.passages);
write_plot_data(vertcat(histograms{:}),fullfile(figureDir,'plot_data','spatial_predictions',char(item.mode),[stem,'.csv']));
files = safe_export(fig,fullfile(figureDir,'spatial_predictions',char(item.mode)),stem,dpi);
end


function write_tables(D, figureDir, selectedLow, commonLow)
tableDir = fullfile(figureDir,'tables'); ensure_dir(tableDir);
featureRows = cell(0,7); locationRows = cell(numel(D.order),5);
for i = 1:numel(D.order)
    zone = get_zone(D,D.order(i));
    locationRows(i,:) = {i,char(D.labels(i)),sprintf('(%.1f, %.2f, %.2f)',zone.x_m,zone.y_m,zone.z_m), ...
        sprintf('%.1f',zone.physical_length_m),double(zone.n_elements)};
    for n = D.windows
        feature = feature_for(D.selected,D.order(i),n);
        featureRows(end+1,:) = {i,D.order(i),D.labels(i),n,feature,feature_label(D,feature), ...
            "tables/retained_location_trials.csv"}; %#ok<AGROW>
    end
end
features = cell2table(featureRows,'VariableNames',{'display_location_number','zone_name','location_label', ...
    'n_passages','feature','feature_label','source_table'});
write_plot_data(features,fullfile(tableDir,'selected_feature_by_location_and_window.csv'));
write_tex_table(fullfile(tableDir,'retained_locations.tex'), ...
    'Retained physical stiffness-loss zones. Coordinates are saved zone centroids.', ...
    'tab:kw51-retained-zones',{'ID','Zone','Centroid $(x,y,z)$ (m)','Length (m)','Elements'},locationRows);
featureTex = [featureRows(:,1),featureRows(:,3),featureRows(:,4),featureRows(:,6)];
write_tex_table(fullfile(tableDir,'selected_features.tex'), ...
    'Saved selected feature families, listed separately for each passage window.', ...
    'tab:kw51-selected-features',{'ID','Zone','Passages','Selected feature'},featureTex);
for n = D.windows
    rows = cell(numel(D.order),2+numel(D.severities)); headers = {'ID','Zone'};
    for s = 1:numel(D.severities), headers{end+1} = sprintf('%g\\%%',D.severities(s)); end %#ok<AGROW>
    for i = 1:numel(D.order)
        feature = feature_for(D.selected,D.order(i),n); rows(i,1:2) = {i,char(D.labels(i))};
        for s = 1:numel(D.severities)
            r = match_row(D.all,feature,D.order(i),n,D.severities(s));
            rows{i,s+2} = sprintf('%d/%d',r.detected_and_correctly_localised_trials,r.total_trials);
        end
    end
    write_tex_table(fullfile(tableDir,sprintf('accepted_trials_%d.tex',n)), ...
        sprintf('Detected and correctly localised trials for %d accumulated passages, using saved selected families.',n), ...
        sprintf('tab:kw51-accepted-%d',n),headers,rows);
end
for mode = ["selected","common"]
    if mode == "selected"
        low = selectedLow; filename = 'low_severity_metrics'; detail = 'saved selected families';
    else
        low = commonLow; filename = 'common_feature_low_severity_metrics';
        detail = ['the common ',char(feature_label(D,D.commonFeature)),' family'];
    end
    rows = cell(numel(D.order)*2,6); records = cell(numel(D.order)*2,1); cursor = 0;
    for i = 1:numel(D.order)
        for n = D.windows
            r = low(low.zone_name == D.order(i) & low.n_passages == n,:);
            cursor = cursor+1;
            rows(cursor,:) = {i,char(D.labels(i)),n,sprintf('%d/%d',r.detected_trials,r.total_trials), ...
                sprintf('%d/%d',r.correctly_localised_trials,r.total_trials), ...
                sprintf('%d/%d',r.detected_and_correctly_localised_trials,r.total_trials)};
            records{cursor} = r;
        end
    end
    caption = sprintf(['Detection (D), ungated correct-zone selection (L), and their joint outcome at ' ...
        '%g\\%% stiffness loss, using %s.'],D.lossPct,tex_escape(detail));
    write_tex_table(fullfile(tableDir,[filename,'.tex']),caption,['tab:kw51-',strrep(filename,'_','-')], ...
        {'ID','Zone','Passages','D','L','D+L'},rows);
    write_plot_data(vertcat(records{:}),fullfile(tableDir,[filename,'.csv']));
end
q = D.all(D.all.feature == D.commonFeature & D.all.severity_pct == D.lossPct,:);
q.shown_in_retained_zone_figures = ismember(q.zone_name,D.order);
write_plot_data(q,fullfile(tableDir,'common_feature_all_zones_at_low_loss.csv'));
end


function write_tex_table(path, caption, label, headers, rows)
columns = ['@{}',repmat('l',1,numel(headers)),'@{}'];
lines = {'\begin{table}[htbp]';'\centering';'\small';'\setlength{\tabcolsep}{4pt}'; ...
    ['\caption{',caption,'}'];['\label{',label,'}'];['\begin{tabular}{',columns,'}']; ...
    '\toprule';[strjoin(headers,' & '),' \\'];'\midrule'};
for i = 1:size(rows,1)
    values = cell(1,size(rows,2));
    for j = 1:size(rows,2), values{j} = tex_escape(rows{i,j}); end
    lines{end+1,1} = [strjoin(values,' & '),' \\']; %#ok<AGROW>
end
lines = [lines;{'\bottomrule';'\end{tabular}';'\end{table}';''}];
write_lines(path,lines);
end


function out = tex_escape(value)
out = char(string(value));
out = strrep(out,'&','\&'); out = strrep(out,'%','\%');
out = strrep(out,'_','\_'); out = strrep(out,'#','\#');
end


function write_plot_data(T, path)
ensure_dir(fileparts(path)); writetable(T,path);
end


function write_lines(path, lines)
ensure_dir(fileparts(path)); fid = fopen(path,'w','n','UTF-8');
if fid < 0, error('KW51:WriteFile','Cannot open %s for writing.',path); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
for i = 1:numel(lines), fprintf(fid,'%s\n',lines{i}); end
end


function ensure_dir(path)
if ~isempty(path) && ~isfolder(path), mkdir(path); end
end
