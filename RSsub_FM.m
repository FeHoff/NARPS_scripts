function RSsub_FM(sub_raw, FMdir_raw, xsub, cat_dir, RSdir, FMdir, targetdir, TR, FWHM, force)
cwd = pwd;

addpath /data/BnB2/TOOLS/spm12
addpath /data/BnB2/TOOLS/DVARS-master/


origdir   = fullfile(RSdir,'Orig');              %   Unprocessed resting-state EPI folder
stdir     = fullfile(RSdir,'ST');                %   Slicetimed data after realigment
normdir   = fullfile(RSdir,'Norm');              %   Normalized resting-state EPI folder, 3D
smoothdir = fullfile(RSdir,['Smooth_' int2str(FWHM) 'mm']);    %   Normalized resting-state EPI folder
T1dir     = cat_dir;   % 3D CAT folder

% --- delete old data if force = true ---
if  force == 1
    try rmdir(RSdir,'s'); end
end

%--------------- get RAW ready for SPM preprocessing---------------------
if (numel(dir(origdir)) <3 && ~exist(fullfile(RSdir,['Counfounds_' xsub '.mat']),'file')) ...
    || force == 1
    % ---------- split 4D NIFTI ---------------
    mkdir(origdir);
    if strcmp(spm_file(sub_raw, 'ext'), 'gz')
        if ~exist(fullfile(RSdir, [xsub '.nii']), 'file')
            copyfile(sub_raw, fullfile(RSdir, [xsub '.nii.gz']))
            gunzip(fullfile(RSdir, [xsub '.nii.gz']))
        end
        spm_file_split(fullfile(RSdir, [xsub '.nii']), origdir)
        delete(fullfile(RSdir, [xsub '.nii.gz']))
        delete(fullfile(RSdir, [xsub '.nii']))
        copyfile(strrep(sub_raw, 'bold', 'sbref'), fullfile(RSdir, [xsub '_sbref.nii.gz']))
        gunzip(fullfile(RSdir, [xsub '_sbref.nii.gz']))
        delete(fullfile(RSdir, [xsub '_sbref.nii.gz']))
    elseif strcmp(spm_file(sub_raw, 'ext'), 'nii')
        copyfile(sub_raw, fullfile(RSdir, [xsub '.nii']))
        spm_file_split(fullfile(RSdir, [xsub '.nii']), origdir)
        delete(fullfile(RSdir, [xsub '.nii']))
    elseif exist(sub_raw,'dir')
        copyfile(fullfile(sub_raw, '*.nii'), origdir)
    else
        return
    end
    % ---------- delete dummys --------------------
    fils = dir(fullfile(origdir,'*.nii')); fils = fils(~[fils.isdir]);
    for ex = 1:4
        delete(fullfile(origdir, fils(ex).name))
    end
end

if numel(dir(FMdir_raw)) > 3
    if numel(dir(FMdir)) <3 || force == 1
        % ---------- split 4D NIFTI ---------------
        mkdir(FMdir);
        fils = dir(FMdir_raw); fils = fils(~[fils.isdir]);
        for i = 1:numel(fils)
            if strcmp(spm_file(fullfile(FMdir_raw, fils(i).name),'ext'),'gz')
                copyfile(fullfile(FMdir_raw, fils(i).name), FMdir);
                gunzip(fullfile(FMdir, fils(i).name));
                delete(fullfile(FMdir, fils(i).name));
            end
        end
    end
end
%%% -------- Proceed with pre-processing -------------------
if ~exist(fullfile(RSdir, [xsub '_meanEPI.nii']), 'file') || force == 1

    if exist(FMdir)
        if isempty(dir(fullfile(FMdir, 'vdm*'))) || force == 1
    % ----- create field map (vdm file) -----
            if force == 1; try
                delete(fullfile(FMdir, 'fpm*'));
                delete(fullfile(FMdir, 'vdm*'));
                delete(fullfile(FMdir, 'sc*'));
            end; end
            % Get Magnitude Difference
            matlabbatch = {};
            matlabbatch{1}.spm.util.imcalc.expression = 'i1-i2';
            fils = dir(fullfile(FMdir, '*magnitude*')); fils = fils(~[fils.isdir]);
            for i = 1:numel(fils)
              matlabbatch{1}.spm.util.imcalc.input{i,1} = fullfile(FMdir, fils(i).name);
            end
            matlabbatch{1}.spm.util.imcalc.output = fullfile(FMdir, strrep(fils(1).name, 'tude1', 'tudediff'));
            spm_jobman('run', matlabbatch)

            matlabbatch = {};
            % Echo times [short TE   long TE]
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [4.92 7.38];
            % Blip direction
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1;
            % Total EPI readout time
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = 29.15;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'vdm5';
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 1;
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 1;
            % CHECK THE FIELDMAP INPUT & OUTPUT
            fils = dir(fullfile(FMdir, '*.*')); fils = fils(~[fils.isdir]);
            if numel(fils)==3
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase{1} = ...
                    fullfile(FMdir, fils(3).name);
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude{1} = ...
                    fullfile(FMdir, fils(1).name);
            elseif numel(fils)==4
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortphase{1} = ...
                    fullfile(FMdir, fils(3).name);
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortmag{1} = ...
                    fullfile(FMdir, fils(1).name);
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longphase{1} = ...
                    fullfile(FMdir, fils(4).name);
                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longmag{1} = ...
                    fullfile(FMdir, fils(2).name);
            end
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template{1} = ...
                fullfile(spm('dir'), 'toolbox', 'Fieldmap', 'T1.nii');
            % fils = dir(fullfile(origdir,'*.nii')); fils = fils(~[fils.isdir]);
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi{1} = ...
                fullfile(RSdir, [xsub '_sbref.nii']);
            matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat{1} = fullfile(T1dir, [xsub '.nii']);
            spm_jobman('run', matlabbatch)
        else
            fprintf('\n --- Fieldmap already calculated --- \n')
            % apply fieldmap to single band reference EPI
            matlabbatch = {};
            matlabbatch{1}.spm.tools.fieldmap.applyvdm.data.scans{1} = fullfile(RSdir, [xsub '_sbref.nii']);
            Fieldmap = dir(fullfile(FMdir,'vdm5*.nii'));
            matlabbatch{1}.spm.tools.fieldmap.applyvdm.data.vdmfile{1} = fullfile(FMdir, Fieldmap.name);
            matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.pedir = 2;
            matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.which = [2 0];
            matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.rinterp = 7;
            spm_jobman('run', matlabbatch)
        end

    % ----- realign & unwarp (using the vdm file) -----
        try
            matlabbatch = {};
            matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.95;
            matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 3;
            matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
            matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 1;
            matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 7;
            matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
            matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
            matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
            matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
            matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
            matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
            matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
            matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
            matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
            matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
            matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
            matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [1 1];
            matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 7;
            matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
            matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
            fils = dir(fullfile(origdir, '*.nii')); fils = fils(~[fils.isdir]);
            for file = 1:size(fils,1)
                matlabbatch{1}.spm.spatial.realignunwarp.data.scans{file,1} = ...
                    fullfile(origdir, fils(file).name);
            end
            Fieldmap = dir(fullfile(FMdir,'vdm5*.nii'));
            matlabbatch{1}.spm.spatial.realignunwarp.data.pmscan{1} = fullfile(FMdir, Fieldmap.name);
            spm_jobman('run', matlabbatch)

            mEPI = dir(fullfile(origdir, 'mean*.nii'));
            movefile(fullfile(origdir, mEPI.name), fullfile(RSdir, [xsub '_meanEPI.nii']))
            rptxt = dir(fullfile(origdir, 'rp*.txt'));
            copyfile(fullfile(origdir, rptxt.name), fullfile(RSdir, ['rp_' xsub '.txt']))
            delete(fullfile(origdir, [xsub '_*.nii']));  % delete unwarped Niftis
            % movefile(fullfile(origdir, 'w*.nii'),RSdir);  % move warped magnitude image


        catch
            fprintf('Problem with Fieldmap of %s: %s \n', xsub, err.message);
        end
    else
        % ---------- realignment --------------------
        matlabbatch = {};
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.95;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 3;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 7;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 7;
        fils = dir(fullfile(origdir, '*.nii')); fils = fils(~[fils.isdir]);
        for file = 1:size(fils,1)
            matlabbatch{1}.spm.spatial.realign.estwrite.data{1}{file,1} = fullfile(origdir, fils(file).name);
        end
        spm_jobman('run', matlabbatch)
        mEPI = dir(fullfile(origdir, 'mean*.nii'));
        movefile(fullfile(origdir, mEPI.name), fullfile(RSdir, [xsub '_meanEPI.nii']))
        rptxt = dir(fullfile(origdir, 'rp*.txt'));
        copyfile(fullfile(origdir, rptxt.name), fullfile(RSdir, ['rp_' xsub '.txt']))
    end
end

%%% !!!! CAREFUL: incorporate slice order manually for now !!!!
if TR > 1.5 && ~exist(fullfile(RSdir, [xsub '_st.nii.gz']), 'file') || TR > 1.5 && force == 1
  %--------------- slicetiming ----------------
    % merge 3D to 4D Nifti
    matlabbatch = {};
    fils = dir(fullfile(origdir, '*.nii')); fils = fils(~[fils.isdir]);
    for file = 1:size(fils,1)
        matlabbatch{1}.spm.util.cat.vols{file,1} = fullfile(origdir, fils(file).name);
    end
    matlabbatch{1}.spm.util.cat.name = fullfile(RSdir, [xsub '.nii']);
    matlabbatch{1}.spm.util.cat.dtype = 4;
    spm_jobman('run', matlabbatch);
    system(['slicetimer -i ' fullfile(RSdir, [xsub '.nii']) ' -o ' fullfile(RSdir, [xsub '_st.nii']) ' -r ' num2str(TR)],'-echo')
    gunzip(fullfile(RSdir, [xsub '_st.nii.gz']));
    mkdir(stdir);
    % delete(fullfile(RSdir, [xsub '_st.nii.gz']));
    spm_file_split(fullfile(RSdir, [xsub '_st.nii']), stdir);
    delete(fullfile(RSdir, [xsub '_st.nii']));
end

if ~exist(fullfile(RSdir, ['w' xsub '_meanEPI.nii']), 'file') || force == 1
    % ---------- Normalise EPIs using Unified Segment -----
    if ~exist(fullfile(RSdir, ['u' xsub '_sbref_seg_sn.mat']), 'file') || force == 1
        % ---------- Affine registration into MNI space --------------------
        % coreg single band EPI to mean EPI
        matlabbatch = {};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = ...
            [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate.source{1} = fullfile(RSdir, ['u' xsub '_sbref.nii']);
        matlabbatch{1}.spm.spatial.coreg.estimate.ref{1} = fullfile(RSdir, [xsub '_meanEPI.nii']);
        spm_jobman('run', matlabbatch)
        % coreg single band EPI to gray matter PM
        matlabbatch = {};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = ...
            [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate.source{1} = fullfile(RSdir, ['u' xsub '_sbref.nii']);
        matlabbatch{1}.spm.spatial.coreg.estimate.ref{1} = fullfile(spm('dir'), 'toolbox', 'OldSeg', 'grey.nii');
        if exist(stdir,'dir')
            fils = dir(fullfile(stdir, '*.nii')); fils = fils(~[fils.isdir]);
            for file = 1:size(fils,1)
                matlabbatch{1}.spm.spatial.coreg.estimate.other{file,1} = fullfile(stdir, fils(file).name);
            end
        else
            fils = dir(fullfile(origdir, '*.nii')); fils = fils(~[fils.isdir]);
            for file = 1:size(fils,1)
                matlabbatch{1}.spm.spatial.coreg.estimate.other{file,1} = fullfile(origdir, fils(file).name);
            end
        end
        matlabbatch{1}.spm.spatial.coreg.estimate.other{file+1,1} = fullfile(RSdir, [xsub '_meanEPI.nii']);
        spm_jobman('run', matlabbatch)

        % --- Normalise using unified segmentation with cutoff of DCT bases (4x5x4) ---
        matlabbatch = {};
        matlabbatch{1}.spm.tools.oldseg.data{1} = fullfile(RSdir, ['u' xsub '_sbref.nii']);
        matlabbatch{1}.spm.tools.oldseg.output.GM = [0 0 0];
        matlabbatch{1}.spm.tools.oldseg.output.WM = [0 0 0];
        matlabbatch{1}.spm.tools.oldseg.output.biascor = 0;
        matlabbatch{1}.spm.tools.oldseg.opts.warpco = 45;   % cutoff
        matlabbatch{1}.spm.tools.oldseg.opts.samp = 2;
        spm_jobman('run', matlabbatch)
    end

    if ~exist(fullfile(RSdir, ['w' xsub '_meanEPI.nii']),'file') || force == 1
     % ---------- Apply deformations --------------------
        mkdir(normdir);
        matlabbatch = {};
        matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname{1} = fullfile(RSdir, ['u' xsub '_sbref_seg_sn.mat']);
        % matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname{1} = fullfile(RSdir,[xsub '_meanEPI_seg_sn.mat']);
        matlabbatch{1}.spm.util.defs.comp{1}.sn2def.vox = [2 2 2];
        if exist(stdir,'dir')
            fils = dir(fullfile(stdir,'*.nii')); fils = fils(~[fils.isdir]);
            for file = 1:size(fils,1)
                matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{file,1} = fullfile(stdir,fils(file).name);
            end
        else
            fils = dir(fullfile(origdir,'*.nii')); fils = fils(~[fils.isdir]);
            for file = 1:size(fils,1)
                matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{file,1} = fullfile(origdir,fils(file).name);
            end
         end
        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr{1} = normdir;
        matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 3;
        spm_jobman('run',matlabbatch)
        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {};
        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{1,1} = fullfile(RSdir,['u' xsub '_sbref.nii']);
        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{2,1} = fullfile(RSdir,[xsub '_meanEPI.nii']);
        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr{1} = RSdir;
        spm_jobman('run',matlabbatch)
    end
end

if exist(normdir,'dir') && ~exist(smoothdir,'dir') || force == 1
% ---------- Smooth --------------------
    matlabbatch = {};
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    matlabbatch{1}.spm.spatial.smooth.fwhm = [FWHM FWHM FWHM];
    fils = dir(fullfile(normdir,'w*.nii')); fils = fils(~[fils.isdir]);
    for file = 1:size(fils,1)
        matlabbatch{1}.spm.spatial.smooth.data{file,1} = fullfile(normdir,fils(file).name);
    end
    spm_jobman('run',matlabbatch)
    mkdir(smoothdir);
    movefile(fullfile(normdir,'sw*.nii'),smoothdir)
end

if ~exist(fullfile(RSdir,['w' xsub '.nii']),'file')
    % merge 3D to 4D Nifti
    matlabbatch = {};
    fils = dir(fullfile(normdir,'w*.nii'));
    for file = 1:size(fils,1)
        matlabbatch{1}.spm.util.cat.vols{file,1} = fullfile(normdir,fils(file).name);
    end
    matlabbatch{1}.spm.util.cat.name = fullfile(RSdir,['w' xsub '.nii']);
    matlabbatch{1}.spm.util.cat.dtype = 4;
    spm_jobman('run',matlabbatch);
    matlabbatch = {};
    fils = dir(fullfile(smoothdir,'sw*.nii'));
    for file = 1:size(fils,1)
        matlabbatch{1}.spm.util.cat.vols{file,1} = fullfile(smoothdir,fils(file).name);
    end
    matlabbatch{1}.spm.util.cat.name = fullfile(RSdir,['s' num2str(FWHM) 'w' xsub '.nii']);
    matlabbatch{1}.spm.util.cat.dtype = 4;
    spm_jobman('run',matlabbatch);
end

%%% redo smoothing with 8mm -- MAYBE IMPLEMENT A LOOP with previous
try rmdir(smoothdir,'s'); end
FWHM = FWHM + 3;
smoothdir = fullfile(RSdir,['Smooth_' int2str(FWHM) 'mm']);    %   smooothed resting-state EPI folder
if exist(normdir,'dir') && ~exist(smoothdir,'dir') || force == 1
% ---------- Smooth --------------------
    matlabbatch = {};
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    matlabbatch{1}.spm.spatial.smooth.fwhm = [FWHM FWHM FWHM];
    fils = dir(fullfile(normdir,'w*.nii')); fils = fils(~[fils.isdir]);
    for file = 1:size(fils,1)
        matlabbatch{1}.spm.spatial.smooth.data{file,1} = fullfile(normdir,fils(file).name);
    end
    spm_jobman('run',matlabbatch)
    mkdir(smoothdir);
    movefile(fullfile(normdir,'sw*.nii'),smoothdir)
end

if ~exist(fullfile(RSdir,['s' num2str(FWHM) 'w' xsub '.nii']),'file')
    % merge 3D to 4D Nifti
    matlabbatch = {};
    fils = dir(fullfile(smoothdir,'sw*.nii'));
    for file = 1:size(fils,1)
        matlabbatch{1}.spm.util.cat.vols{file,1} = fullfile(smoothdir,fils(file).name);
    end
    matlabbatch{1}.spm.util.cat.name = fullfile(RSdir,['s' num2str(FWHM) 'w' xsub '.nii']);
    matlabbatch{1}.spm.util.cat.dtype = 4;
    spm_jobman('run',matlabbatch);
end

if exist(fullfile(T1dir,'mri',['wp0' xsub '.nii']),'file') && ...
        ~exist(fullfile(RSdir,['wp0' xsub '.nii']),'file')
    % resample partial volume label image to 2x2x2mm
    matlabbatch = {};
    matlabbatch{1}.spm.spatial.coreg.write.ref = ...
        {fullfile(RSdir,['w' xsub '_meanEPI.nii'])};
    matlabbatch{1}.spm.spatial.coreg.write.source = ...
        {fullfile(T1dir,'mri',['wp0' xsub '.nii'])};
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch)
    movefile(fullfile(T1dir,'mri',['rwp0' xsub '.nii']),fullfile(RSdir,['wp0' xsub '.nii']))
end

if exist(fullfile(T1dir,'mri',['mwp1' xsub '.nii']),'file') && ...
        ~exist(fullfile(RSdir,['mwp1' xsub '.nii']),'file')
    % resample modulated GM segment image to 2x2x2mm
    matlabbatch = {};
    matlabbatch{1}.spm.spatial.coreg.write.ref = ...
        {fullfile(RSdir,['w' xsub '_meanEPI.nii'])};
    matlabbatch{1}.spm.spatial.coreg.write.source = ...
        {fullfile(T1dir,'mri',['mwp1' xsub '.nii'])};
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 7; % 7th degree B-Spline
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch)
    movefile(fullfile(T1dir,'mri',['rmwp1' xsub '.nii']),fullfile(RSdir,['mwp1' xsub '.nii']))
end

if ~exist(fullfile(RSdir,['Counfounds_' xsub '.mat']),'file') || force == 1 || force == 2
    clear pXYZ tX tY tZ gx gx2 rp rpp A dat msk gx2 reg y
    % getting & computing movement parameter (FD & sqrt) from realignment parameter
    txt  = dir(fullfile(RSdir,'rp_*.txt'));
    rp   = load(fullfile(RSdir,txt.name));
    [FDts,FD_Stat]=FDCalc(rp);
    rp   = rp-repmat(mean(rp),size(rp,1),1);
    drp = [zeros(1,6); diff(rp)];
    move = sqrt(sum(drp(:,1:3).^2,2));
    % compute root mean squared movement & frame wise displacement
    euler = zeros(size(rp,1),1);
    for ii=1:size(rp,1)
        euler(ii) = 50*acos((cos(drp(ii,4))*cos(drp(ii,5)) + cos(drp(ii,4))*cos(drp(ii,6)) + ...
            cos(drp(ii,5))*cos(drp(ii,6)) + sin(drp(ii,4))*sin(drp(ii,6))*sin(drp(ii,5)) - 1)/2);
    end
    RMS  = sqrt(nanmean([move+euler].^2));
    tmp  = sum(abs([drp(:,1:3) drp(:,4:6)*50]),2);
    FD   = sqrt(nanmean(tmp.^2));
    % data file
    % Vi  = spm_vol(fullfile(fsmoothdir,['sw' xsub '.nii']));
    Vi  = spm_vol(fullfile(RSdir,['w' xsub '.nii']));
    % computing global, GM, WM, CSF, mean signals
    % getting partial volume mask
    msk = spm_read_vols(spm_vol(fullfile(RSdir,['wp0' xsub '.nii'])));
    % msk: min 0.0157; max 2.9961; 1=CSF, 2=GM, 3=WM
    idxGM  = find(msk>1.99 & msk<2.01);
    % GMmsk = msk;
    % GMmsk(msk<1.99 | msk>2.01) = 0;
    % GMmsk = spm_erode(GMmsk);            % erode once
    % % GMmsk = spm_erode(spm_erode(GMmsk)); % erode twice
    % idxGM  = find(GMmsk>0);
    % idxWM  = find(msk>2.99);
    idxCSF = find(msk>0.99 & msk<1.01);
    gx     = nan(4,numel(Vi));
    %   i.e. aCompCor does >.99 WM 2x erode & >.99 CSF NN clustering criteria
%     WMmsk = msk;
    WMmsk(msk<2.99) = 0;
    WMmsk = spm_erode(WMmsk);            % erode once
    % WMmsk = spm_erode(spm_erode(WMmsk)); % erode twice
    idxWM  = find(WMmsk>0);
    datWM  = zeros(length(idxWM),numel(Vi));
    datCSF = zeros(length(idxCSF),numel(Vi));
    datALL = zeros(length(find(msk>0)),numel(Vi));
    for i=1:numel(Vi)
        dat     = spm_read_vols(Vi(i));
        gx(1,i) = mean(dat(idxGM));
        datWM(:,i)  = dat(idxWM);    % WM voxels only
        gx(2,i) = mean(datWM(:,i));
        datCSF(:,i) = dat(idxCSF);   % CSF voxels only
        gx(3,i) = mean(datCSF(:,i));
        datALL(:,i) = dat(msk>0);    % All brain voxels
        gx(4,i) = mean(dat(msk>0));
    end

    PCA_WM  = pca(datWM - repmat(mean(datWM),size(datWM,1),1));
    PCA_CSF = pca(datCSF - repmat(mean(datCSF),size(datCSF,1),1));

    gx2(1,:) = gx(1,:)-mean(gx(1,:));
    gx2(2,:) = gx(2,:)-mean(gx(2,:));
    gx2(3,:) = gx(3,:)-mean(gx(3,:));
    gx2(4,:) = gx(4,:)-mean(gx(4,:));

    reg = [gx2' gx2'.^2 rp rp.^2 drp drp.^2 PCA_WM(:,1:7) PCA_CSF(:,1:7)];
    reg = reg-repmat(mean(reg),size(reg,1),1);
    reg = reg./repmat(std(reg),size(reg,1),1);

    [DVARS,DVARS_Stat]=DVARSCalc(datALL,'scale',1/100,'RDVARS');
    badTP=find(DVARS_Stat.pvals<0.05./(numel(DVARS)));
    [V,DSE_Stat]=DSEvars(datALL,'scale',1/100);

    %   FIGURE with BOLD intensity carpet-plot
    %     fMRIDiag_plot(V,DVARS_Stat,'BOLD',datALL,'ColRng',[-3 3],'FD',FDts,'AbsMov',[FD_Stat.AbsRot FD_Stat.AbsTrans]);
    %   Diagnosis plot based on stored variables
    %     f_hdl=figure('position',[50,500,600,600]);
    %     fMRIDiag_plot(V,DVARS_Stat,'FD',FDts,'AbsMov',[FD_Stat.AbsRot FD_Stat.AbsTrans],'TickScaler',0.5,'figure',f_hdl);

    save(fullfile(RSdir,['Counfounds_' xsub '.mat']),'reg','rp','drp','gx',...
        'gx2','TR','RMS','FD','DVARS','DVARS_Stat','V','DSE_Stat','FDts','FD_Stat','badTP')

end

if exist(fullfile(RSdir,['w' xsub '.nii']),'file')
    try rmdir(origdir,'s'); end
    try rmdir(stdir,'s'); end
    try rmdir(normdir,'s'); end
    try rmdir(smoothdir,'s'); end
end

if exist(fullfile(RSdir,['Counfounds_' xsub '.mat']),'file')
    mkdir(targetdir);
    movefile(fullfile(RSdir,'*'),targetdir);
    rmdir(RSdir);
end

% exit
end
