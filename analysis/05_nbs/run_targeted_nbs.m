function run_targeted_nbs(stage_dir, network_pair, direction)
% Purpose: Run targeted nbs.
% Inputs: Network-pair design and MAT inputs.
% Outputs: Baseline NBS component MAT files.
% Manuscript: S3.

toolbox_dir = getenv('NBS_HOME');
addpath(toolbox_dir);
input_dir = fullfile(stage_dir, 'outputs', 'inputs', network_pair);
output_dir = fullfile(stage_dir, 'outputs', 'raw');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
if strcmp(direction, 'FM_gt_HC')
    contrast = '[0 1 0]';
elseif strcmp(direction, 'HC_gt_FM')
    contrast = '[0 -1 0]';
else
    error('Unknown direction');
end
rng(str2double(getenv('PERM_SEED')), 'twister');
UI.method.ui = 'Run NBS';
UI.test.ui = 't-test';
UI.size.ui = 'Intensity';
UI.thresh.ui = getenv('NBS_BASE_THRESHOLD');
UI.perms.ui = getenv('NBS_PERMUTATIONS');
UI.alpha.ui = '1';
UI.contrast.ui = contrast;
UI.design.ui = fullfile(input_dir, 'design.txt');
UI.matrices.ui = fullfile(input_dir, 'matrices.mat');
UI.exchange.ui = '';
UI.node_coor.ui = '';
UI.node_label.ui = '';
global nbs
clear global nbs
global nbs
S = struct();
NBSrun(UI, S);
test_stat = nbs.NBS.test_stat;
component_matrices = nbs.NBS.con_mat;
component_pvalues = nbs.NBS.pval;
save(fullfile(output_dir, [network_pair '_' direction '.mat']), 'test_stat', 'component_matrices', 'component_pvalues', 'UI', '-v7');
end
