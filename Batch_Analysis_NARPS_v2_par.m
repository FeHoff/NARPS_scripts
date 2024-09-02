clear all
addpath /data/BnB2/TOOLS/spm12

global defaults,
spm('defaults','FMRI')

defaults.stats.fmri.t  = 16;
defaults.stats.fmri.t0 = 1;


cwd = pwd;


workdir         = '/data/BnB_TEMP/Data_NARPS/';
datapath        = '/data/BnB_TEMP/Data_NARPS/Derivatives/';


Si = dir(fullfile(datapath,'MGT1')); Si = Si([Si.isdir]); Si = Si(3:end);
for i=1:size(Si,1)
    subjects{i} = Si(i).name;
end

MODELL = 'NARPS_gain1st_censBadTP'; 
%MODELL = 'NARPS_loss1st';

labels = {'MGT1', 'MGT2', 'MGT3', 'MGT4'};  
results_folder = fullfile(workdir, '/SingleSubjectAnalysis', MODELL);
mkdir(results_folder);

TR          =  1.0;

for sub= 1:numel(subjects)
	  if ~exist(fullfile(results_folder,subjects{sub},'SPM.mat'),'file')
            setupSPM(subjects{sub}, workdir, results_folder, labels, cwd, TR);
    end
end


% parallelize model estimation
for sub= 1:numel(subjects)
    if ~exist(fullfile(results_folder,subjects{sub},'RPV.nii'),'file')
        % try
          estimateSPM(fullfile(results_folder, subjects{sub}));
        % end
    end
end

cd(cwd)

%quit;
