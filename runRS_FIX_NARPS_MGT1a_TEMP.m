clear, clc; warning off all
%%%
%       EPI preprocessing
%
%       Felix Hoffstaedter
%       Data and Platforms
%       INM-7 - Brain and Behaviour
%       March 2018 - Research Centre Juelich
%
%       calling: RSsub_FIXed(sub_raw{sub},xsubs{sub},targetdir{sub},cat_dir{sub},TR,FWHM,force)
%%%

% ----------  define sample to process   ----------------------------
sample = 'NARPS';
% ---------- Repition time for band-pass filter TR  -----------------------
TR = 2;   % get from json !!!
% ---------- smoothing kernel in mm  --------------------------------------
FWHM = 5;
% ---------- Define RS scan  -----------------------
scan = 'rest_1';   % get from json !!!
% ---------- Define [out '_fsl']put folder  -----------------------
out   = 'RS'; clean = 'FIX';
% force  --- reprocessing of ready subjects [0 = NOT reprocess] -----
force = 0;

% ------  define slurm call with max RAM & wall clock time   -------
call_slurm = ['srun -p large -N 1 -n 1 -c 2 --mem-per-cpu=7500 --job-name=' sample];
% call_slurm = ['srun -p large -N 1 -n 1 --time=6:00:00 --exclude=brainb32a --job-name=' sample];
call_matlab = ' /nfsusr/local/MATLAB/R2016a/bin/matlab -nodisplay -nodesktop -nosplash -r';

% ---------- get subject and sessions  ------------------------------
Dirs = fullfile('/data/BnB1/Raw_Data', sample);
% scaninfo = jsondecode(fileread(fullfile(Dirs,[scan '.json'])));
% ---------- CAT [out '_fsl']put folder  -------------------------------------
target = fullfile('/data/BnB2/Derivatives/func', sample);
target2 = fullfile('/data/BnB_TEMP/Derivatives/func', sample);
cat   = fullfile('/data/BnB2/Derivatives/CAT/12.5', sample);
Si = dir(fullfile(Dirs,'*')); Si = Si([Si.isdir]); Si = Si(3:end);

if strfind(Si(3).name,'sub')
    Di = Si;
    fprintf('Sample "%s" contains ~%d subjects \n', sample, size(Di,1));
    cnt=1;
    for i=1:size(Di,1)
        if strfind(Di(i).name,'sub')
            subdir = fullfile(Dirs, Di(i).name);
            Dii = dir(fullfile(subdir,'*'));
            try Dii = Dii([Dii.isdir]); Dii = Dii(3:end); end
            if strfind(Dii(1).name,'ses')
                for ii=1:size(Dii,1)
                    if exist(fullfile(target, [out '_fsl'], Di(i).name, Dii(ii).name, [clean '_func_data.nii.gz']),'file')
                        sub_raw{cnt} = fullfile(target, [out '_fsl'], Di(i).name, Dii(ii).name, [clean '_func_data.nii.gz']);
                        xsubs{cnt} = [Di(i).name '_' Dii(ii).name];
                        targetdir{cnt} = fullfile(target2, [out '_' clean], Di(i).name, Dii(ii).name);
                        cat_dir{cnt} = fullfile(cat, Di(i).name, Dii(ii).name);
                        cnt=cnt+1;
                    end
                end
            else
                if exist(fullfile(target, [out '_fsl'], Di(i).name, [clean '_func_data.nii.gz']),'file')
                    xsubs{cnt} = Di(i).name;
                    sub_raw{cnt} = fullfile(target, [out '_fsl'], Di(i).name, [clean '_func_data.nii.gz']);
                    targetdir{cnt} = fullfile(target2, [out '_' clean], Di(i).name);
                    cat_dir{cnt} = fullfile(cat, Di(i).name);
                    cnt=cnt+1;
                end
            end
        end
    end
else
    fprintf('%s Sites found \n', num2str(size(Si,1)));
    cnt=1;
    for s = 3 %1:size(Si,1)
        site = Si(s).name;
        Di = dir(fullfile(Dirs, site, '*')); Di = Di([Di.isdir]); Di = Di(3:end);
        fprintf('%s: %d subjects \n', site, size(Di,1));
        for i = 1:size(Di,1)
            if strfind(Di(i).name,'sub')
                subdir = fullfile(Dirs, site , Di(i).name);
                Dii = dir(fullfile(subdir,'*'));
                try Dii = Dii([Dii.isdir]); Dii = Dii(3:end);
                    if strfind(Dii(1).name,'ses')
                        for ii = 1:size(Dii,1)
                            if exist(fullfile(target, site, [out '_fsl'], Di(i).name, Dii(ii).name, [clean '_func_data.nii.gz']),'file')
                              sub_raw{cnt} = fullfile(target, site, [out '_fsl'], Di(i).name, Dii(ii).name, [clean '_func_data.nii.gz']);
                              xsubs{cnt} = [Di(i).name '_' Dii(ii).name];
                              targetdir{cnt} = fullfile(target2, site, [out '_' clean], Di(i).name, Dii(ii).name);
                              cat_dir{cnt} = fullfile(cat, site, Di(i).name, Dii(ii).name);
                              cnt=cnt+1;
                            end
                        end
                    else
                        if exist(fullfile(target, site, [out '_fsl'], Di(i).name, [clean '_func_data.nii.gz']),'file')
                            xsubs{cnt} = Di(i).name;
                            sub_raw{cnt} = fullfile(target, site,  [out '_fsl'], Di(i).name, [clean '_func_data.nii.gz']);
                            targetdir{cnt} = fullfile(target2, site,  [out '_' clean], Di(i).name);
                            cat_dir{cnt} = fullfile(cat, site, Di(i).name);
                            cnt=cnt+1;
                        end
                    end
                end
            end
        end
    end
end

% ----------  check for already processed subjects  -----------------------
idx= []; i=0;
if force
    idx = 1 : numel(xsubs);
else
  for sub = 1:numel(xsubs)
      if exist(sub_raw{sub},'file') && ...
        ~exist(fullfile(targetdir{sub}, ['Counfounds_' xsubs{sub} '.mat']),'file')
              i = i + 1;
              idx(i) = sub;
      end
  end
end
if numel(idx) == 0
    fprintf('No subjects of %i to process! Aborting.', numel(xsubs)); return; 
else
    fprintf('Found %i unprocessed of %i subjects..\n', numel(idx), numel(xsubs));
end
sub_raw = sub_raw(idx);
targetdir = targetdir(idx);
xsubs = xsubs(idx);
cat_dir = cat_dir(idx);
%  
pause

% ----------  submit slurm jobs  -----------------------
for sub = 1:numel(xsubs)
  fprintf('Submit EPI job for subject %s \n', xsubs{sub});
    system([call_slurm ' --error=slurm_log2/RS_FIX_SSD_' xsubs{sub} ...
        '_err.%j --output=slurm_log2/RS_FIX_SSD_' xsubs{sub} '_out.%j' call_matlab ...
        ' "RSsub_FIXed_SSD(''' sub_raw{sub} ''' ,' '''' ...
        xsubs{sub} ''' ,' '''' targetdir{sub} ''',' '''' cat_dir{sub} ''' ,' '' num2str(TR) ...
        '  ,' '' num2str(FWHM) '  ,' '' num2str(force) ' )" &' ],'-echo');
end

for sub = 1:numel(xsubs)
  fprintf('Submit EPI job for subject %s \n', xsubs{sub});
    system([call_matlab ' "RSsub_FIXed_SSD(''' sub_raw{sub} ''' ,' '''' ...
        xsubs{sub} ''' ,' '''' targetdir{sub} ''',' '''' cat_dir{sub} ''' ,' '' num2str(TR) ...
        '  ,' '' num2str(FWHM) '  ,' '' num2str(force) ' )" &' ],'-echo');
end


% parfor sub = 1:numel(xsubs)
%     try RSsub_FIXed_SSD(sub_raw{sub},xsubs{sub},targetdir{sub},cat_dir{sub},TR,FWHM,force); end
% end