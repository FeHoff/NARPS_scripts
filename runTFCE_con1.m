%------run TFCE for NARPS -----------------------------------------------------------------
matlabbatch = {};
matlabbatch{1}.spm.tools.tfce_estimate.spmmat = {'/data/BnB_TEMP/Data_NARPS/ANOVA/NARPS_Gain1st_censBadTP/SPM.mat'};
matlabbatch{1}.spm.tools.tfce_estimate.mask = '';
matlabbatch{1}.spm.tools.tfce_estimate.conspec.titlestr = '';
matlabbatch{1}.spm.tools.tfce_estimate.conspec.contrasts = 1;
matlabbatch{1}.spm.tools.tfce_estimate.conspec.n_perm = 5000;
matlabbatch{1}.spm.tools.tfce_estimate.nuisance_method = 2;
matlabbatch{1}.spm.tools.tfce_estimate.tbss = 0;
matlabbatch{1}.spm.tools.tfce_estimate.E_weight = 0.5;
matlabbatch{1}.spm.tools.tfce_estimate.singlethreaded = 0;
spm_jobman('run', matlabbatch)
