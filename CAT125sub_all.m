function CAT125sub_all(sub_raw, xsub, TMPdir, targetdir, force)
%%%
%       CAT12 segmentation using SHOOTING & DARTEL
%
%       script based on:
%       $Id: cat_defaults_expert.m 895 2016-03-11 16:07:09Z gaser $
%       $Id: cat_defaults.m 1183 2017-09-08 16:45:08Z dahnke $
%
%       Felix Hoffstaedter
%       INM-7 - Brain and Behaviour
%       Data and Platforms
%       March 2018 - Research Centre Juelich
%
%       CATsub_all( <subjectfolder>, <subjectname>, <TEMPdirectory>, <targetdirectory>, <force processing> )
%%%

warning off all

addpath /data/BnB2/TOOLS/spm12_cat12.5

% --------------- start script only if T1w exists ---------------
if exist(sub_raw,'file')

% --- how many subjects are processed in parallel [in window = 0 ] --
jobs  = 0;    %  a job will use min 4 cores for i.e. denoising
% --- how many cores are use for fractal dimension estimation      --
cores  = 0;   %  1 if multiple subjects are processed
% ------- smoothing of modulated volume data im mm [5 & 8 default] --
vbmFWHM = [5 8];
% ---------- estimate fractal dimensions - SLOW !  ------------------
FD      = 1;
% ---------- smoothing of cortical thickness in mm [>=15] -----------
thkFWHM = [0 15];
% ---------- smoothing of folding parameters in mm [>=25] -----------
srfFWHM = [0 25];
%--------------- load cat12 expert defaults -------------------------
def = fullfile(spm_file(which('CAT125sub_all'), 'path'), 'cat125_defaults_bnb');
try global cat; if cat.extopts.expertgui == 0; spm fmri; cat12(def); end
catch;  spm fmri; cat12(def); end

sub_nii  = fullfile(TMPdir, [xsub '.nii']);          %   Subjects T1 file for VBM

% --- delete old data if force = true ---
if  force == 1
    try rmdir(TMPdir,'s'); end
end

if ~exist(sub_nii,'file')
    % ----------  copy & extract nifti file to CAT dir -----------------
    mkdir(TMPdir);

    if strfind(sub_raw,'run-02')
    % --- realignment and mean image creation if multiple T1w acquired ---
    % ----- write loop for more than 2 acquisitions per scan session -----
        copyfile(sub_raw, fullfile(TMPdir, [xsub '-1.nii.gz']));
        gunzip(fullfile(TMPdir, [xsub '-1.nii.gz']));
        delete(fullfile(TMPdir, [xsub '-1.nii.gz']));
        sub_raw2 = strrep(sub_raw,'run-01','run-02');
        copyfile(sub_raw2, fullfile(TMPdir, [xsub '-2.nii.gz']));
        gunzip(fullfile(TMPdir, [xsub '-2.nii.gz']));
        delete(fullfile(TMPdir, [xsub '-2.nii.gz']));
        matlabbatch = {};
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fullfile(TMPdir, [xsub '-1.nii'])};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {fullfile(TMPdir, [xsub '-2.nii'])};
        spm_jobman('run', matlabbatch)
        matlabbatch = {};
        matlabbatch{1}.spm.util.imcalc.input = {fullfile(TMPdir, [xsub '-1.nii'])
                                                fullfile(TMPdir, [xsub '-2.nii'])};
        matlabbatch{1}.spm.util.imcalc.output = [xsub '.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {TMPdir};
        matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
        spm_jobman('run', matlabbatch)
    else
        if strcmp(spm_file(sub_raw,'ext'),'gz')
            copyfile(sub_raw, fullfile(TMPdir, [xsub '.nii.gz']));
            gunzip(fullfile(TMPdir, [xsub '.nii.gz']));
            delete(fullfile(TMPdir, [xsub '.nii.gz']));
        elseif strcmp(spm_file(sub_raw,'ext'),'nii')
            copyfile(sub_raw, sub_nii);
        else
            return
        end
    end

    % for robustness set (0 0 0) coordinate to middle of image
    V = spm_vol(sub_nii);
    % pre-estimated COM of MNI template
    com_reference = [0 -20 -15];
    Affine = eye(4);
    vol = spm_read_vols(V);
    avg = mean(vol(:));
    avg = mean(vol(vol>avg));
    % don't use background values
    [x,y,z] = ind2sub(size(vol),find(vol>avg));
    com = V.mat(1:3,:)*[mean(x) mean(y) mean(z) 1]';
    com = com';
    M = spm_get_space(V.fname);
    Affine(1:3,4) = (com - com_reference)';
    spm_get_space(V.fname,Affine\M);

end

% ---------- high-dimensional SHOOTING & DARTEL normalization -----------------
if ~exist(fullfile(TMPdir, 'mri', ['m0wp1' xsub '.nii']), 'file') || force == 1
    matlabbatch = {};
    matlabbatch{1}.spm.tools.cat.estwrite.nproc = jobs;
    % data to process
    matlabbatch{1}.spm.tools.cat.estwrite.data{1} = [sub_nii ',1'];
    % run things
    spm_jobman('run', matlabbatch)
end

% ---------- Smooth modulated gray & white matter segments --------------------
if ~exist(fullfile(TMPdir, 'mri', ['s' num2str(vbmFWHM(1)) 'm0wp1' xsub '.nii']),'file') ...
        || force == 1
    matlabbatch = {};
    matlabbatch{1}.spm.spatial.smooth.data = {
                               fullfile(TMPdir, 'mri', ['mwp1' xsub '.nii'])
                               fullfile(TMPdir, 'mri', ['m0wp1' xsub '.nii'])};
    for i = 1:numel(vbmFWHM)
        matlabbatch{1}.spm.spatial.smooth.fwhm = [vbmFWHM(i) vbmFWHM(i) vbmFWHM(i)];
        matlabbatch{1}.spm.spatial.smooth.prefix = ['s' num2str(vbmFWHM(i))];
        spm_jobman('run', matlabbatch)
    end
end

% ---------- Estimate total intracranial volume (TIV) & -------------------
% ---------- global GM, WM, CSF, WM-hyperintensity vol --------------------
if (~exist(fullfile(TMPdir, ['TIV_' xsub '.txt']),'file') || force == 1) && ...
  exist(fullfile(TMPdir, 'report', ['cat_' xsub '.xml']), 'file');   % only with existing 'report/cat_subj.mat'
    matlabbatch = {};
    matlabbatch{1}.spm.tools.cat.tools.calcvol.data_xml = ...
    {fullfile(TMPdir, 'report', ['cat_' xsub '.xml'])};
    matlabbatch{1}.spm.tools.cat.tools.calcvol.calcvol_TIV = 0;
    matlabbatch{1}.spm.tools.cat.tools.calcvol.calcvol_name = [ 'TIV_' xsub '.txt'];
    spm_jobman('run', matlabbatch);
    movefile(fullfile(pwd, ['TIV_' xsub '.txt']), fullfile(TMPdir, ['TIV_' xsub '.txt']))
end

if exist(fullfile(TMPdir, 'surf', ['lh.central.' xsub '.gii']),'file')
    % ---------- Extract gyrification, cortical complexity and sulcus depth  --------------------
    if ~exist(fullfile(TMPdir, 'surf', ['rh.sqrtsulc.' xsub]),'file') || force == 1
        matlabbatch = {};
        matlabbatch{1}.spm.tools.cat.stools.surfextract.data_surf = {
                                fullfile(TMPdir, 'surf', ['lh.central.' xsub '.gii'])};
        matlabbatch{1}.spm.tools.cat.stools.surfextract.GI = 1;
        matlabbatch{1}.spm.tools.cat.stools.surfextract.FD = FD;
        matlabbatch{1}.spm.tools.cat.stools.surfextract.SD = 1;
        matlabbatch{1}.spm.tools.cat.stools.surfextract.nproc = cores;
        spm_jobman('run', matlabbatch)
    end

    % ---------- Resample and smooth cortical thickness  --------------------
    if ~exist(fullfile(TMPdir, 'surf', ['s' num2str(thkFWHM(1)) '.rh.thickness.resampled.' xsub '.gii']),'file') ...
            || force == 1
        matlabbatch = {};
        matlabbatch{1}.spm.tools.cat.stools.surfresamp.data_surf = {
                                fullfile(TMPdir, 'surf', ['lh.thickness.' xsub])};
        matlabbatch{1}.spm.tools.cat.stools.surfresamp.nproc = cores;
        for i = 1:numel(thkFWHM)
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.fwhm_surf = thkFWHM(i);
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.merge_hemi = 0;
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.mesh32k = 0;
            spm_jobman('run', matlabbatch)
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.mesh32k = 1;
            spm_jobman('run', matlabbatch)
        end
    end

    % ---------- Resample and smooth gyrification and sulcus depth  --------------------
    if ~exist(fullfile(TMPdir, 'surf', ['s' num2str(srfFWHM(1)) '.rh.sqrtsulc.resampled.' xsub '.gii']),'file') ...
            || force == 1
        matlabbatch = {};
        if FD == 0
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.data_surf = {
                                    fullfile(TMPdir, 'surf', ['lh.gyrification.' xsub])
                                    fullfile(TMPdir, 'surf', ['lh.sqrtsulc.' xsub])};
        else
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.data_surf = {
                                    fullfile(TMPdir, 'surf', ['lh.gyrification.' xsub])
                                    fullfile(TMPdir, 'surf', ['lh.sqrtsulc.' xsub])
                                    fullfile(TMPdir, 'surf', ['lh.fractaldimension.' xsub])};
        end
        matlabbatch{1}.spm.tools.cat.stools.surfresamp.nproc = cores;
        for i = 1:numel(srfFWHM)
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.fwhm_surf = srfFWHM(i);
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.merge_hemi = 0;
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.mesh32k = 0;
            spm_jobman('run', matlabbatch)
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
            matlabbatch{1}.spm.tools.cat.stools.surfresamp.mesh32k = 1;
            spm_jobman('run', matlabbatch)
        end
    end
end

if exist(fullfile(TMPdir, 'surf', ['s' num2str(srfFWHM(numel(srfFWHM))) ...
        '.mesh.sqrtsulc.resampled_32k.' xsub '.gii']),'file')
    mkdir(targetdir);
    movefile(fullfile(TMPdir,'*'),targetdir);
    rmdir(TMPdir);
end

end
exit
end
