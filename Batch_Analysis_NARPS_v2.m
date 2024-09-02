clear all
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

supra_results_folder = [pwd filesep 'SingleSubjectAnalysis'];
mkdir (supra_results_folder, MODELL)
results_folder = ['/data/BnB_TEMP/Data_NARPS/SingleSubjectAnalysis' filesep MODELL];

TR          =  1.0;
define      = 'secs'; 

for sub= 1:numel(subjects)
  work_dir{sub}=fullfile(results_folder,subjects{sub});
  try
  if ~exist(fullfile(work_dir{sub},'SPM.mat'),'file')
    err = 0;
    mkdir(work_dir{sub}); delete([work_dir{sub} '/*.*']);
    datarun=[];

    runs = labels;
    runs = {runs{[1 2 3 4]}};
    
    nruns           = numel(runs);
    for j = 1:nruns; datarun{j}=fullfile(datapath,runs{j},subjects{sub}); end % bei 4-D Nifti


    
   clear soas
    load(fullfile(workdir, '/Behavior/Logfiles/SOAfiles/', [subjects{sub} '.mat']));
  
    for i=1:nruns        
        indx = nan(1,numel(labels));
        Labs = {};
        for ii=1:numel(labels)

            for iii=1:numel(soas)
            
                try
                if strcmp(soas(iii).name{1},labels{ii})
                    indx(ii) = iii;
                    Labs{ii} = [labels{ii}];
                end
                end
                
            end
        end
        
        SESS{i,1} = Labs{i};
        SESS{i,2} = indx(i);
        
    end

    if any(isnan([SESS{:,2}]))
        huhu Error!!
    end
    
    
    clear nscans
     
    P = '';
    for run = 1:nruns
        files = spm_select('expand', fullfile(datarun{run},['s5w' subjects{sub} '.nii'])); % 4-D Nifti 
        P  = [P; files];
        nscans(run)  = size(files,1);
    end

    clear SPM; SPM       = struct('nscan',nscans);
    SPM.xY.P  = P;
    SPM.xY.RT = TR;

    SPM.xBF.UNITS  = define;  SPM.xBF.name  = 'hrf (with time derivative)';  % Option: 'hrf'
    SPM.xBF.T      = 16;       SPM.xBF.T0    = 8;       SPM.xBF.dt = SPM.xY.RT/SPM.xBF.T;
    SPM.xBF.length = 32.2;     SPM.xBF.order = 2;       SPM.xBF = spm_get_bf(SPM.xBF);      SPM.xBF.Volterra   = 1;

    for ses = 1:numel(nscans)
        
            SPM.Sess(ses).U.ons  = soas(SESS{ses,2}).ons;
            SPM.Sess(ses).U.dur  = soas(SESS{ses,2}).dur;
            SPM.Sess(ses).U.name = {SESS{ses,1}};
            if numel(soas(SESS{ses,2}).par1')==0
                SPM.Sess(ses).U.P(1).name = 'none';
            else
                SPM.Sess(ses).U.P(1).name = 'Gain';
                SPM.Sess(ses).U.P(1).P = soas(SESS{ses,2}).par1';
                SPM.Sess(ses).U.P(1).h = 1;  % 1 = lin. Modulation
            end
                
            if numel(soas(SESS{ses,2}).par2')==0
                SPM.Sess(ses).U.P(2).name = 'none';                
            else
                SPM.Sess(ses).U.P(2).name = 'Loss';
                SPM.Sess(ses).U.P(2).P = soas(SESS{ses,2}).par2';
                SPM.Sess(ses).U.P(2).h = 1;  % 1 = lin. Modulation
            end
                
            if numel(soas(SESS{ses,2}).par3')==0
                SPM.Sess(ses).U.P(3).name = 'none';                
            else
                SPM.Sess(ses).U.P(3).name = 'RT';
                SPM.Sess(ses).U.P(3).P = soas(SESS{ses,2}).par3';
                SPM.Sess(ses).U.P(3).h = 1;  % 1 = lin. Modulation
            end                  
             
             
            if any(SPM.Sess(ses).U.dur < 0) 
                err = 1;
            elseif any(diff(round(SPM.Sess(ses).U.ons)) < 1)  % war mal 5
                err = 1;
            end
                
            rp = dir(fullfile(datapath,runs{ses},subjects{sub},'rp_*.txt'));
            [r1,r2,r3,r4,r5,r6] = textread(fullfile(datapath,runs{ses},subjects{sub},rp.name),'%f%f%f%f%f%f');

            
            Conf   = load(fullfile(datapath,runs{ses},subjects{sub},['Counfounds_' subjects{sub} '.mat']));
            censTP = zeros(numel(r1),1);
            censTP(Conf.badTP,1) = 1;
            [r1,r2,r3,r4,r5,r6] = textread(fullfile(datapath,runs{ses},subjects{sub},smoothedfolder,rp.name),'%f%f%f%f%f%f');

            if sum(censTP) > 0
                SPM.Sess(ses).C.C    = [r1,r2,r3,r4,r5,r6,censTP]; 
                SPM.Sess(ses).C.name = {'x', 'y', 'z', 'yaw', 'pitch', 'roll', 'badTP'};
            else
                SPM.Sess(ses).C.C    = [r1,r2,r3,r4,r5,r6]; 
                SPM.Sess(ses).C.name = {'x', 'y', 'z', 'yaw', 'pitch', 'roll'};
            end
    end
        
    SPM.xGX.iGXcalc = 'none';
    SPM.xVi.form = 'AR(1) + w'; 
    SPM.xX.K(1).HParam  = 128;

    cd(work_dir{sub});
    SPM = spm_fmri_spm_ui(SPM);

    save SPM.mat SPM;

    cd(cwd)
  end
  end
end


% parallelize model estimation

parfor sub= 1:numel(subjects)
    try
        estimateSPM(work_dir{sub});
    end
end

cd(cwd)

%quit;