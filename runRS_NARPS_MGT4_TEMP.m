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
% ---------- Repition time for band-pass filter TR  -----------------------
TR = 1;   % get from json !!!
% ---------- smoothing kernel in mm  --------------------------------------
FWHM = 5;
% ---------- Define RS scan  -----------------------
scan = 'MGT_run-04';   % get from json !!!
% ---------- Define output folder  -----------------------
out  = 'MGT4';
% force  --- reprocessing of ready subjects [0 = NOT reprocess] -----
force = 0;

% ------  define slurm call with max RAM & wall clock time   -------
call_slurm = ['srun -p large -N 1 -n 1 -c 2 --mem-per-cpu=5500 --time=13:00:00 --job-name=' sample];
% call_slurm = ['srun -p large -N 1 -n 1 --time=6:00:00 --exclude=brainb32a --job-name=' sample];
call_matlab = ' /nfsusr/local/MATLAB/R2016a/bin/matlab -nodisplay -nodesktop -nosplash -r';

% ---------- get subject and sessions  ------------------------------
RAWdirs = fullfile('/data/BnB_TEMP/Data_NARPS/', sample);
% scaninfo = jsondecode(fileread(fullfile(Dirs,[scan '.json'])));
% ---------- CAT output folder  -------------------------------------
target = fullfile('/data/BnB_TEMP/Data_NARPS/Derivatives');
TMP    = fullfile('/data/BnB_TEMP/Data_NARPS/Derivatives/tmp');
CAT    = fullfile('/data/BnB_TEMP/Data_NARPS/Derivatives/CAT/12.5');
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
                    if exist(fullfile(subdir, Dii(ii).name, 'func'),'dir')
                        fil = dir(fullfile(subdir, Dii(ii).name, 'func', ...
                                [Di(i).name '_' Dii(ii).name '_task-' scan '_bold.nii.gz']));
                        if numel(fil) > 0
                            sub_raw{cnt} = fullfile(subdir, Dii(ii).name, 'func', fil(1).name);
                            xsubs{cnt} = [Di(i).name '_' Dii(ii).name];
                            TMPdir{cnt} = fullfile(TMP, out, Di(i).name, Dii(ii).name);
                            cat_dir{cnt} = fullfile(CAT, Di(i).name, Dii(ii).name);
                            targetdir{cnt} = fullfile(target, out, Di(i).name, Dii(ii).name);
                            if exist(fullfile(subdir, Dii(ii).name, 'fmap'),'dir')
                                FMdir_raw{cnt} = fullfile(subdir, Dii(ii).name, 'fmap');
                                FMdir{cnt} = fullfile(target, 'FM', Di(i).name, Dii(ii).name);
                            else
                                FMdir_raw{cnt} = '';
                                FMdir{cnt} = '';                   
                            end
                            cnt=cnt+1;
                        end

                    end
                end
            else
                if exist(fullfile(subdir, 'func'),'dir')
                    fil = dir(fullfile(subdir, 'func', [Di(i).name '_task-' scan '_bold.nii.gz']));
                    if numel(fil) > 0
                        sub_raw{cnt} = fullfile(subdir, 'func', fil(1).name);
                        xsubs{cnt} = Di(i).name;
                        TMPdir{cnt} = fullfile(TMP, out, Di(i).name);
                        cat_dir{cnt} = fullfile(CAT, Di(i).name);
                        targetdir{cnt} = fullfile(target, out, Di(i).name);
                        if exist(fullfile(subdir, 'fmap'),'dir')
                            FMdir_raw{cnt} = fullfile(subdir, 'fmap');
                            FMdir{cnt} = fullfile(target, 'FM', Di(i).name);
                        else
                            FMdir_raw{cnt} = '';
                            FMdir{cnt} = '';                   
                        end
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
                        if exist(fullfile(subdir, Dii(ii).name, 'func'),'dir')
                            fil = dir(fullfile(subdir, Dii(ii).name, 'func', ...
                                    [Di(i).name '_' Dii(ii).name '_task-' scan '_bold.nii.gz']));
                            if numel(fil) > 0
                                sub_raw{cnt} = fullfile(subdir, Dii(ii).name, 'func', fil(1).name);
                                xsubs{cnt} = [Di(i).name '_' Dii(ii).name];
                                TMPdir{cnt} = fullfile(TMP, site, out, Di(i).name, Dii(ii).name);
                                cat_dir{cnt} = fullfile(CAT, site, Di(i).name, Dii(ii).name);
                                targetdir{cnt} = fullfile(target, site, out, Di(i).name, Dii(ii).name);
                                if exist(fullfile(subdir, Dii(ii).name, 'fmap'),'dir')
                                    FMdir_raw{cnt} = fullfile(subdir, Dii(ii).name, 'fmap');
                                    FMdir{cnt} = fullfile(target, site, 'FM', Di(i).name, Dii(ii).name);
                                else
                                    FMdir_raw{cnt} = '';
                                    FMdir{cnt} = '';                   
                                end
                                cnt=cnt+1;
                            end
                         end
                    end
                else
                    if exist(fullfile(subdir, 'func'),'dir')
                        fil = dir(fullfile(subdir, 'func', [Di(i).name '_task-' scan '_bold.nii.gz']));
                        if numel(fil) > 0
                            sub_raw{cnt} = fullfile(subdir, 'func', fil(1).name);
                            xsubs{cnt} = Di(i).name;
                            TMPdir{cnt} = fullfile(TMP, site, out, Di(i).name);
                            cat_dir{cnt} = fullfile(CAT, site, Di(i).name);
                            targetdir{cnt} = fullfile(target, site, out, Di(i).name);
                            if exist(fullfile(subdir, 'fmap'),'dir')
                                FMdir_raw{cnt} = fullfile(subdir, 'fmap');
                                FMdir{cnt} = fullfile(target, site, 'FM', Di(i).name);
                            else
                                FMdir_raw{cnt} = '';
                                FMdir{cnt} = '';                   
                            end
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
        if exist(sub_raw{sub},'file') && exist(cat_dir{sub},'dir') && ...
          ~exist(fullfile(targetdir{sub}, ['Counfounds_' xsubs{sub} '.mat']),'file')
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
FMdir_raw = FMdir_raw(idx);
xsubs     = xsubs(idx);
cat_dir   = cat_dir(idx);
TMPdir    = TMPdir(idx);
FMdir     = FMdir(idx);
targetdir = targetdir(idx);

pause
for sub = 1:numel(xsubs)
  fprintf('Submit EPIprep job for subject %s \n', xsubs{sub});
    system([call_slurm ' --error=slurm_log/RS_' xsubs{sub} call_matlab ...
        ' "RSsub_FM(''' sub_raw{sub} ''' ,' '''' FMdir_raw{sub} ''' ,' '''' ...
        xsubs{sub} ''' ,' '''' cat_dir{sub} ''' ,' '''' TMPdir{sub} ''',' '''' ...
        FMdir{sub} ''',' '''' targetdir{sub} ''', ' '' num2str(TR) ...
        '  ,' '' num2str(FWHM) '  ,' '' num2str(force) ' )" &' ],'-echo');
end

% for sub = 1:numel(xsubs)
%   fprintf('Submit EPI job for subject %s \n', xsubs{sub});
%     system([call_matlab ' "RSsub_FM(''' sub_raw{sub} ''' ,' '''' FMdir_raw{sub} ''' ,' '''' ...
%         xsubs{sub} ''' ,' '''' cat_dir{sub} ''' ,' '''' TMPdir{sub} ''',' '''' ...
%         FMdir{sub} ''',' '''' targetdir{sub} ''', ' '' num2str(TR) ...
%         '  ,' '' num2str(FWHM) '  ,' '' num2str(force) ' )" &' ],'-echo');
% end

% for sub = 1:numel(xsubs)
%     try
%         RSsub_FM(sub_raw{sub}, FMdir_raw{sub}, xsubs{sub}, cat_dir{sub}, ...
%             TMPdir{sub}, FMdir{sub}, targetdir{sub}, TR, FWHM, force);
%     end
% end
