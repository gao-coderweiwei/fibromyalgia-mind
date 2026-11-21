function segment(manifest, spm_path)
% Purpose: Segment.
% Inputs: T1w manifest and SPM TPM.
% Outputs: Smoothed 2 mm grey matter maps.
% Manuscript: S6.

addpath(spm_path);
spm('defaults', 'FMRI');
spm_jobman('initcfg');
inputs = readlines(manifest);
inputs = inputs(strlength(inputs) > 0);
for input = inputs'
    [folder, filename] = fileparts(input);
    output = fullfile(folder, 'gm_s8.nii');
    if isfile(output)
        continue;
    end
    modulated = fullfile(folder, 'mwc1' + filename + ".nii");
    if ~isfile(modulated)
        batch = struct();
        batch.channel.vols = {char(input + ",1")};
        batch.channel.biasreg = 0.001;
        batch.channel.biasfwhm = 60;
        batch.channel.write = [0 0];
        gaussians = [1 1 2 3 4 2];
        for tissue = 1:6
            batch.tissue(tissue).tpm = {sprintf('%s,%d', fullfile(spm_path, 'tpm', 'TPM.nii'), tissue)};
            batch.tissue(tissue).ngaus = gaussians(tissue);
            batch.tissue(tissue).native = [0 0];
            batch.tissue(tissue).warped = [0 tissue == 1];
        end
        batch.warp = struct('mrf', 1, 'cleanup', 1, 'reg', [0 .001 .5 .05 .2], ...
            'affreg', 'mni', 'fwhm', 0, 'samp', 3, 'write', [0 0], ...
            'vox', str2double(getenv('VBM_VOXEL_MM')), 'bb', [-90 -126 -72; 90 90 108]);
        matlabbatch = {struct('spm', struct('spatial', struct('preproc', batch)))};
        spm_jobman('run', matlabbatch);
    end
    temporary = fullfile(folder, 'gm_s8_partial.nii');
    spm_smooth(char(modulated), char(temporary), repmat(str2double(getenv('VBM_SMOOTH_MM')), 1, 3), 16);
    movefile(temporary, output);
end
end
