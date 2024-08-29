clear, clc; warning off all
%%%
%       EPI preprocessing
%
%       Felix Hoffstaedter
%       Data and Platforms
%       INM-7 - Brain and Behaviour
%       March 2018 - Research Centre Juelich
%
%       calling: CATsub_all( <path2raw.nii>, <subjectname(_session)>, <targetdir>, force )
%%%

% ----------  define sample to process   ----------------------------
sample = 'NARPS';
% ---------- Define T1w scan  -----------------------
scan = 'T1w';
% ---------- Define output folder  -----------------------
out  = 'MGT1';
% force  --- reprocessing of ready subjects [0 = NOT reprocess] -----
force = 0;

% ------  define slurm call with max RAM & wall clock time   -------
call_slurm = ['srun -p large -N 1 -n 1 --mem-per-cpu=5500 --time=9:00:00 --job-name=' sample];
% call_slurm = ['srun -p large -N 1 -n 1 --time=6:00:00 --exclude=brainb32a --job-name=' sample];
call_matlab = ' /nfsusr/local/MATLAB/R2016a/bin/matlab -nodisplay -nodesktop -nosplash -r';

% ---------- get subject and sessions  ------------------------------
RAWdirs = fullfile('/data/BnB_TEMP/Data_NARPS/', sample);
% scaninfo = jsondecode(fileread(fullfile(Dirs,[scan '.json'])));
% ---------- CAT output folder  -------------------------------------
target = fullfile('/data/BnB_TEMP/Data_NARPS/Derivatives/CAT/12.5');
TMP    = fullfile('/data/BnB_TEMP/Data_NARPS/Derivatives/CATtmp/12.5');
Si = dir(fullfile(RAWdirs,'*')); Si = Si([Si.isdir]); Si = Si(3:end);

if strfind(Si(3).name,'sub')
    Di = Si;
    fprintf('Sample "%s" contains ~%d subjects \n', sample, size(Di,1));
    cnt=1;
    for i=1:size(Di,1)
        if strfind(Di(i).name,'sub')
            subdir = fullfile(RAWdirs, Di(i).name);
            Dii = dir(fullfile(subdir,'*'));
            try Dii = Dii([Dii.isdir]); Dii = Dii(3:end); end
            if strfind(Dii(1).name,'ses')
                for ii=1:size(Dii,1)
                    if exist(fullfile(subdir, Dii(ii).name, 'anat'),'dir')
                        fil = dir(fullfile(subdir, Dii(ii).name, 'anat', ...
                                [Di(i).name '_' Dii(ii).name '_' scan '.nii.gz']));
                        if numel(fil) > 0
                            sub_raw{cnt} = fullfile(subdir, Dii(ii).name, 'anat', fil(1).name);
                            xsubs{cnt} = [Di(i).name '_' Dii(ii).name];
                            TMPdir{cnt} = fullfile(TMP, Di(i).name, Dii(ii).name);
                            targetdir{cnt} = fullfile(target, Di(i).name, Dii(ii).name);
                            cnt=cnt+1;
                        end

                    end
                end
            else
                if exist(fullfile(subdir, 'anat'),'dir')
                    fil = dir(fullfile(subdir, 'anat', [Di(i).name '_' scan '.nii.gz']));
                    if numel(fil) > 0
                       sub_raw{cnt} = fullfile(subdir, 'anat', fil(1).name);
                       xsubs{cnt} = Di(i).name;
                       TMPdir{cnt} = fullfile(TMP, Di(i).name);
                       targetdir{cnt} = fullfile(target, Di(i).name);
                       cnt=cnt+1;
                    end
                end
            end
        end
    end
else
    fprintf('%s Sites found \n', num2str(size(Si,1)));
    cnt=1;
    for s = 1:size(Si,1)
        site = Si(s).name;
        Di = dir(fullfile(RAWdirs, site, '*')); Di = Di([Di.isdir]); Di = Di(3:end);
        fprintf('%s: %d subjects \n', site, size(Di,1));
        for i = 1:size(Di,1)
            if strfind(Di(i).name,'sub')
                subdir = fullfile(RAWdirs, site , Di(i).name);
                Dii = dir(fullfile(subdir,'*'));
                try Dii = Dii([Dii.isdir]); Dii = Dii(3:end); end
                if strfind(Dii(1).name,'ses')
                    for ii=1:size(Dii,1)
                        if exist(fullfile(subdir, Dii(ii).name, 'anat'),'dir')
                            fil = dir(fullfile(subdir, Dii(ii).name, 'anat', ...
                                    [Di(i).name '_' Dii(ii).name '_' scan '.nii.gz']));
                            if numel(fil) > 0
                                sub_raw{cnt} = fullfile(subdir, Dii(ii).name, 'anat', fil(1).name);
                                xsubs{cnt} = [Di(i).name '_' Dii(ii).name];
                                TMPdir{cnt} = fullfile(TMP, site, Di(i).name, Dii(ii).name);
                                targetdir{cnt} = fullfile(target, site, Di(i).name, Dii(ii).name);
                                cnt=cnt+1;
                            end
                         end
                    end
                else
                    if exist(fullfile(subdir, 'anat'),'dir')
                        fil = dir(fullfile(subdir, 'anat', [Di(i).name '_' scan '.nii.gz']));
                        if numel(fil) > 0
                            sub_raw{cnt} = fullfile(subdir, 'anat', fil(1).name);
                            xsubs{cnt} = Di(i).name;
                            TMPdir{cnt} = fullfile(TMP, site, Di(i).name);
                            targetdir{cnt} = fullfile(target, site, Di(i).name);
                            cat_dir{cnt} = fullfile(CAT, site, Di(i).name);
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
%         if exist(sub_raw{sub},'file') && exist(cat_dir{sub},'dir') && ...
        if exist(sub_raw{sub},'file') && ...
                (~exist(fullfile(targetdir{sub}, 'surf', ['s25.rh.sqrtsulc.resampled.' xsubs{sub} '.gii']),'file') || ...
                            ~exist(fullfile(targetdir{sub}, ['TIV_' xsubs{sub} '.txt']),'file'))
                i = i + 1;
                idx(i) = sub;
        end
    end
end

if numel(idx) == 0
    fprintf('No subjects to process! Aborting.'); return;
else
    fprintf('Found %i unprocessed of %i subjects..\n', numel(idx), numel(xsubs));
end

sub_raw   = sub_raw(idx);
xsubs     = xsubs(idx);
TMPdir    = TMPdir(idx);
targetdir = targetdir(idx);

pause
% for sub = 1:numel(xsubs)
%   fprintf('Submit CAT12 job for subject %s \n', xsubs{sub});
%     system([call_slurm ' --error=slurm_log/CAT_' xsubs{sub} call_matlab ...
%         ' "CAT125sub_all(''' sub_raw{sub} ''' ,' '''' xsubs{sub} ''' ,' '''' ...
%         TMPdir{sub} ''',' '''' targetdir{sub} ''', ' '' num2str(force) ' )" &' ],'-echo');
% end

for sub = 1:numel(xsubs)
  fprintf('Submit EPI job for subject %s \n', xsubs{sub});
    system([call_matlab ' "CAT125sub_all(''' sub_raw{sub} ''' ,' '''' xsubs{sub} ''' ,' '''' ...
        TMPdir{sub} ''',' '''' targetdir{sub} ''', ' '' num2str(force) ' )" &' ],'-echo');
end

% for sub = 1:numel(xsubs)
%     try
%         CAT125sub_all(sub_raw{sub}, xsubs{sub}, TMPdir{sub}, targetdir{sub}, force);
%     end
% end
