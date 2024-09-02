warning off

addpath /data/BnB2/TOOLS/spm12/
clear all

global defaults,
spm('Defaults','FMRI')

cwd = pwd;

workdir        = '/data/BnB_TEMP/Data_NARPS/';        

SAMPLE_SELECTION_TYPE = 3; 

%1st level
MODELL_1st = 'NARPS_loss1st';
%MODELL_1st = 'NARPS_loss1st';

%2nd level
MODELL = 'NARPS_Loss1st_varequal';
%MODELL = 'NARPS_Loss1st'; 


%load setup.mat
load setup_subvarequal.mat  

% specify number of groups and conditions

groups = 2;        

betas = [5:8]; 

datapath = fullfile(workdir, 'SingleSubjectAnalysis', MODELL_1st);
statpath = fullfile(workdir, 'ANOVA', MODELL);

mkdir(statpath);

if SAMPLE_SELECTION_TYPE == 3
    
    Si = dir(datapath); Si = Si([Si.isdir]); Si = Si(3:end);
    for i=1:size(Si,1)
        subs{i} = Si(i).name;
    end
    for i=1:size(subs,2)
        xsub{i,1} = subs{i};
        xsub{i,2} = mod(str2num(xsub{i,1}(1,5:end)),2);
    end            
    
    D1 = xsub(find(cell2mat(xsub(:,2))==1),1);  % 1 = odd no. = equalIndifference 
    D2 = xsub(find(cell2mat(xsub(:,2))==0),1);  % 0 = even no. = equalRange 
    
    for i=1:numel(D1)
        subjects{i} = D1{i};
    end
    
    for i=1:numel(D2)
        subjects2{i} = D2{i};
    end
    
end % SAMPLE_SELECTION_TYPE


for sub= 1:numel(subjects);
    
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sub).scans = {};
    
    for scan = 1:numel(betas)
        
        if betas(scan)<10
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sub).scans{scan,1} = ...
                fullfile(datapath, subjects{sub}, ['con_000' int2str(betas(scan)) '.nii']);
        else
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sub).scans{scan,1} = ...
                fullfile(datapath, subjects{sub}, [con_00' int2str(betas(scan)) '.nii']);
        end
    end
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sub).conds = [1:numel(betas)];
end


if groups >=2 
    
    for sub= 1:numel(subjects2);
        
    sub2= numel(subjects)+sub;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sub2).scans = {};
    for scan = 1:numel(betas)
        
        if betas(scan)<10
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sub2).scans{scan,1} = ...
                fullfile(datapath, subjects2{sub}, ['con_000' int2str(betas(scan)) '.nii']);
        else
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sub2).scans{scan,1} = ...
                fullfile(datapath, subjects2{sub}, ['con_00' int2str(betas(scan)) '.nii']);
        end
    end
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sub2).conds = ...
        [(numel(betas))+1:2*(numel(betas))];
end
end


matlabbatch{1}.spm.stats.factorial_design.dir{1} = [statpath];


save( [statpath filesep 'setupANOVA.mat'],'matlabbatch');
spm_jobman('run',matlabbatch)
    
load('estimate.mat');
jobs{1}.stats{1}.fmri_est.spmmat{1} = [statpath filesep 'SPM.mat'];
spm_jobman('run',jobs)

