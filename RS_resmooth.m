function RS_resmooth(sub_raw, FMdir_raw, xsub, cat_dir, targetdir, FMdir, RSdir, TR, FWHM, force)
cwd = pwd;

addpath /data/BnB2/TOOLS/spm12
addpath /data/BnB2/TOOLS/DVARS-master/

FWHM = 5;

origdir   = fullfile(RSdir,'Orig');              %   Unprocessed resting-state EPI folder
stdir     = fullfile(RSdir,'ST');                %   Slicetimed data after realigment
normdir   = fullfile(RSdir,'Norm');              %   Normalized resting-state EPI folder, 3D
smoothdir = fullfile(RSdir,['Smooth_' int2str(FWHM) 'mm']);    %   Normalized resting-state EPI folder
T1dir     = cat_dir;   % 3D CAT folder


if exist(fullfile(RSdir,['sw' xsub '.nii']),'file')
    movefile(fullfile(RSdir,['sw' xsub '.nii']),fullfile(RSdir,['s5w' xsub '.nii']))
end

if ~exist(normdir,'dir') || force == 1
     % ---------- split 4D NIFTI ---------------
    mkdir(normdir);
    matlabbatch = {};
    matlabbatch{1}.spm.util.split.vol = {fullfile(RSdir,['w' xsub '.nii'])};
    matlabbatch{1}.spm.util.split.outdir = {normdir};
    spm_jobman('run',matlabbatch)
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
    movefile(fullfile(normdir,'sw*.nii'), smoothdir)
end

if ~exist(fullfile(RSdir,['s' num2str(FWHM) 'w' xsub '.nii']),'file') && ...
        exist(smoothdir,'dir')
    % merge 3D to 4D Nifti
    matlabbatch = {};
    fils = dir(fullfile(smoothdir,['sw*.nii']));
    for file = 1:size(fils,1)
        matlabbatch{1}.spm.util.cat.vols{file,1} = fullfile(smoothdir,fils(file).name);
    end
    matlabbatch{1}.spm.util.cat.name = fullfile(RSdir,['s' num2str(FWHM) 'w' xsub '.nii']);
    matlabbatch{1}.spm.util.cat.dtype = 4;
    spm_jobman('run',matlabbatch);
end

% if exist(fullfile(RSdir,['s' num2str(FWHM) 'w' xsub '.nii']),'file')
%     try rmdir(smoothdir,'s'); end
% end

% if exist(fullfile(T1dir,'mri',['wp0' xsub '.nii']),'file') && ...
%         ~exist(fullfile(RSdir,['wp0' xsub '.nii']),'file')
%     % resample partial volume label image to 2x2x2mm
%     matlabbatch = {};
%     matlabbatch{1}.spm.spatial.coreg.write.ref = ...
%         {fullfile(RSdir,['w' xsub '_meanEPI.nii'])};
%     matlabbatch{1}.spm.spatial.coreg.write.source = ...
%         {fullfile(T1dir,'mri',['wp0' xsub '.nii'])};
%     matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
%     matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
%     matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
%     matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
%     spm_jobman('run',matlabbatch)
%     movefile(fullfile(T1dir,'mri',['rwp0' xsub '.nii']),fullfile(RSdir,['wp0' xsub '.nii']))
% end
% % 
% if exist(fullfile(T1dir,'mri',['mwp1' xsub '.nii']),'file') && ...
%         ~exist(fullfile(RSdir,['mwp1' xsub '.nii']),'file')
%     % resample modulated GM segment image to 2x2x2mm
%     matlabbatch = {};
%     matlabbatch{1}.spm.spatial.coreg.write.ref = ...
%         {fullfile(RSdir,['w' xsub '_meanEPI.nii'])};
%     matlabbatch{1}.spm.spatial.coreg.write.source = ...
%         {fullfile(T1dir,'mri',['mwp1' xsub '.nii'])};
%     matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 7; % 7th degree B-Spline
%     matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
%     matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
%     matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
%     spm_jobman('run',matlabbatch)
%     movefile(fullfile(T1dir,'mri',['rmwp1' xsub '.nii']),fullfile(RSdir,['mwp1' xsub '.nii']))
% end
% 
% if ~exist(fullfile(RSdir,['Counfounds_' xsub '.mat']),'file') || force == 1 || force == 2
%     clear pXYZ tX tY tZ gx gx2 rp rpp A dat msk gx2 reg y
%     % getting & computing movement parameter (FD & sqrt) from realignment parameter
%     txt  = dir(fullfile(RSdir,'rp_*.txt'));
%     rp   = load(fullfile(RSdir,txt.name));
%     [FDts,FD_Stat]=FDCalc(rp);
%     rp   = rp-repmat(mean(rp),size(rp,1),1);
%     drp = [zeros(1,6); diff(rp)];
%     move = sqrt(sum(drp(:,1:3).^2,2));
%     % compute root mean squared movement & frame wise displacement
%     euler = zeros(size(rp,1),1);
%     for ii=1:size(rp,1)
%         euler(ii) = 50*acos((cos(drp(ii,4))*cos(drp(ii,5)) + cos(drp(ii,4))*cos(drp(ii,6)) + ...
%             cos(drp(ii,5))*cos(drp(ii,6)) + sin(drp(ii,4))*sin(drp(ii,6))*sin(drp(ii,5)) - 1)/2);
%     end
%     RMS  = sqrt(nanmean([move+euler].^2));
%     tmp  = sum(abs([drp(:,1:3) drp(:,4:6)*50]),2);
%     FD   = sqrt(nanmean(tmp.^2));
%     % data file
%     % Vi  = spm_vol(fullfile(fsmoothdir,['sw' xsub '.nii']));
%     Vi  = spm_vol(fullfile(RSdir,['w' xsub '.nii']));
%     % computing global, GM, WM, CSF, mean signals
%     % getting partial volume mask
%     msk = spm_read_vols(spm_vol(fullfile(RSdir,['wp0' xsub '.nii'])));
%     % msk: min 0.0157; max 2.9961; 1=CSF, 2=GM, 3=WM
%     idxGM  = find(msk>1.99 & msk<2.01);
%     % GMmsk = msk;
%     % GMmsk(msk<1.99 | msk>2.01) = 0;
%     % GMmsk = spm_erode(GMmsk);            % erode once
%     % % GMmsk = spm_erode(spm_erode(GMmsk)); % erode twice
%     % idxGM  = find(GMmsk>0);
%     % idxWM  = find(msk>2.99);
%     idxCSF = find(msk>0.99 & msk<1.01);
%     gx     = nan(4,numel(Vi));
%     %   i.e. aCompCor does >.99 WM 2x erode & >.99 CSF NN clustering criteria
%     WMmsk = msk;
%     WMmsk(msk<2.99) = 0;
%     WMmsk = spm_erode(WMmsk);            % erode once
%     % WMmsk = spm_erode(spm_erode(WMmsk)); % erode twice
%     idxWM  = find(WMmsk>0);
%     datWM  = zeros(length(idxWM),numel(Vi));
%     datCSF = zeros(length(idxCSF),numel(Vi));
%     datALL = zeros(length(find(msk>0)),numel(Vi));
%     for i=1:numel(Vi)
%         dat     = spm_read_vols(Vi(i));
%         gx(1,i) = mean(dat(idxGM));
%         datWM(:,i)  = dat(idxWM);    % WM voxels only
%         gx(2,i) = mean(datWM(:,i));
%         datCSF(:,i) = dat(idxCSF);   % CSF voxels only
%         gx(3,i) = mean(datCSF(:,i));
%         datALL(:,i) = dat(msk>0);    % All brain voxels
%         gx(4,i) = mean(dat(msk>0));
%     end
% 
%     PCA_WM  = pca(datWM - repmat(mean(datWM),size(datWM,1),1));
%     PCA_CSF = pca(datCSF - repmat(mean(datCSF),size(datCSF,1),1));
% 
%     gx2(1,:) = gx(1,:)-mean(gx(1,:));
%     gx2(2,:) = gx(2,:)-mean(gx(2,:));
%     gx2(3,:) = gx(3,:)-mean(gx(3,:));
%     gx2(4,:) = gx(4,:)-mean(gx(4,:));
% 
%     reg = [gx2' gx2'.^2 rp rp.^2 drp drp.^2 PCA_WM(:,1:7) PCA_CSF(:,1:7)];
%     reg = reg-repmat(mean(reg),size(reg,1),1);
%     reg = reg./repmat(std(reg),size(reg,1),1);
% 
%     [DVARS,DVARS_Stat]=DVARSCalc(datALL,'scale',1/100,'RDVARS');
%     badTP=find(DVARS_Stat.pvals<0.05./(numel(DVARS)));
%     [V,DSE_Stat]=DSEvars(datALL,'scale',1/100);
% 
%     %   FIGURE with BOLD intensity carpet-plot
%     %     fMRIDiag_plot(V,DVARS_Stat,'BOLD',datALL,'ColRng',[-3 3],'FD',FDts,'AbsMov',[FD_Stat.AbsRot FD_Stat.AbsTrans]);
%     %   Diagnosis plot based on stored variables
%     %     f_hdl=figure('position',[50,500,600,600]);
%     %     fMRIDiag_plot(V,DVARS_Stat,'FD',FDts,'AbsMov',[FD_Stat.AbsRot FD_Stat.AbsTrans],'TickScaler',0.5,'figure',f_hdl);
% 
%     save(fullfile(RSdir,['Counfounds_' xsub '.mat']),'reg','rp','drp','gx',...
%         'gx2','TR','RMS','FD','DVARS','DVARS_Stat','V','DSE_Stat','FDts','FD_Stat','badTP')
% 
% end

% if exist(fullfile(RSdir,['w' xsub '.nii']),'file')
%     try rmdir(origdir,'s'); end
%     try rmdir(stdir,'s'); end
%     try rmdir(normdir,'s'); end
%     try rmdir(smoothdir,'s'); end
% end

% if exist(fullfile(RSdir,['Counfounds_' xsub '.mat']),'file')
%     mkdir(targetdir);
%     movefile(fullfile(RSdir,'s*'),targetdir);
%     rmdir(RSdir);
% end

% exit
end
