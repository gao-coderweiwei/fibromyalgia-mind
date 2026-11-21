function run_thresholds()
% Purpose: Run thresholds.
% Inputs: Baseline NBS MAT files and design matrices.
% Outputs: Threshold and corrected component table.
% Manuscript: S3.

stage = getenv('THRESHOLD_STAGE');
addpath(getenv('NBS_HOME'));
original = fullfile(getenv('NBS_STAGE'), 'outputs');
base_threshold = str2double(getenv('NBS_BASE_THRESHOLD'));
pairs = {'SalVentAttn-Default', 'SalVentAttn-DorsAttn'};
directions = {'FM_gt_HC', 'HC_gt_FM'};
contrasts = {'[0 1 0]', '[0 -1 0]'};
rows = cell(0, 10);

global nbs
for threshold = sscanf(getenv('NBS_THRESHOLDS'), '%f')'
    for pair_index = 1:numel(pairs)
        pair = pairs{pair_index};
        input_dir = fullfile(original, 'inputs', pair);
        for direction_index = 1:numel(directions)
            direction = directions{direction_index};
            baseline = load(fullfile(original, 'raw', [pair '_' direction '.mat']));
            if threshold == base_threshold
                components = baseline.component_matrices;
                p_values = baseline.component_pvalues;
                stats = baseline.test_stat;
            else
                UI = baseline.UI;
                UI.thresh.ui = num2str(threshold);
                UI.contrast.ui = contrasts{direction_index};
                UI.design.ui = fullfile(input_dir, 'design.txt');
                UI.matrices.ui = fullfile(input_dir, 'matrices.mat');
                nbs = [];
                rng(str2double(getenv('PERM_SEED')), 'twister');
                evalc('NBSrun(UI, struct());');
                components = nbs.NBS.con_mat;
                p_values = nbs.NBS.pval;
                stats = nbs.NBS.test_stat;
            end

            if isempty(components)
                rows(end + 1, :) = {pair, direction, threshold, NaN, 0, 0, NaN, NaN, NaN, NaN};
            end
            for component_index = 1:numel(components)
                adjacency = components{component_index};
                [edge_from, edge_to] = find(adjacency);
                edges = find(adjacency);
                nodes = unique([edge_from; edge_to]);
                mass = sum(stats(edges) - threshold);
                p_fwe = p_values(component_index);
                overlap = NaN;
                if direction_index == 1
                    reference_edges = find(baseline.component_matrices{1});
                    overlap = numel(intersect(edges, reference_edges));
                end
                % Correct across two masks and two directions at this threshold.
                rows(end + 1, :) = {pair, direction, threshold, component_index, ...
                    numel(edges), numel(nodes), mass, p_fwe, min(1, str2double(getenv('NBS_CORRECTION_TESTS')) * p_fwe), overlap};
            end
        end
    end
end

results = cell2table(rows, 'VariableNames', {'network_pair', 'direction', 'threshold', ...
    'component_id', 'n_edges', 'n_nodes', 'mass', 'p_fwe', 'p_bonferroni_within_threshold', ...
    'n_edges_in_original_positive_component'});
writetable(results, fullfile(stage, 'nbs_threshold_results.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
end
