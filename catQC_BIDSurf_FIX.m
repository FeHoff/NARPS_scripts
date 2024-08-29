clear, clc; warning off all

%--------------- readout of path-to-data ----------------------------
datadir = fullfile('/data/BnB_TEMP/Data_NARPS/Derivatives/');

%--------------- select sample directory containing preprocessed data ----------------------------
inpt = spm_input('Input directories or Lookupfile?','!+1','b','Directory|Lookup.mat',[1 2],1,'Select mode');

% --------------- select type of data to do covariance analysis ----------------------------
data = spm_input('Choose data to QC ->','!+1','b','mri|surf|MGT1|MGT2|MGT3|MGT4',[1 2 3 4 5 6],1,'Select mode');
datafolder = {'' '' '' '' '' ''};
datatype = {'CAT/12.5' 'CAT/12.5' 'MGT1_FIX' 'MGT2_FIX' 'MGT3_FIX' 'MGT4_FIX'};
cnt = 0; Pcnt = 0;

%--------------- Gap between voxels to be considered in covariance: default = 3 -------------------
%--------------- lower values take very long -------------------
if data ~= 2
    matlabbatch{1}.spm.tools.cat.tools.check_cov.gap = 2;
end

if inpt == 1
%     Dirs = spm_select(1,'dir','Select sample directory to check homogeniety','',fullfile(datadir,datatype{data}));
    Dirs = fullfile(datadir,datatype{data});
    Di = dir(fullfile(Dirs,datafolder{data},'*')); Di = Di([Di.isdir]); Di = Di(3:end);
    for i=1:size(Di,1)
        subjects{i} = fullfile(Dirs,datafolder{data},Di(i).name);
        xsub{i} = Di(i).name;
    end    
elseif inpt == 2
    Dirs = spm_select(1,'.*\.mat$','Select Lookup file to check homogeniety','');
    load(Dirs);
    
    cova = spm_input('Include nuisance variable?','!+1','b','--|Age|ICV|GMV|FD',[0 1 2 3 4],1,'Select mode');
    
    for i=1:numel(sub)
        subjects{i} = fullfile(datadir,datatype{data},SubDir{i},datafolder{data},sub{i});
    end
    xsub = sub;
end

for sub = 1:numel(xsub)
    % ----- define files checke covariance: 3D = modulated / RS = normalized -----
    datafile = {['mri/m0wp1' xsub{sub} '.nii'] ['surf/mesh.thickness.resampled_32k.' xsub{sub} '.gii'] ...
        ['w' xsub{sub} '_meanEPI.nii'] ['w' xsub{sub} '_meanEPI.nii'] ['w' xsub{sub} '_meanEPI.nii'] ['w' xsub{sub} '_meanEPI.nii']};
%         ['wu' xsub{sub} '_sbref.nii'] ['wu' xsub{sub} '_sbref.nii'] ['wu' xsub{sub} '_sbref.nii'] ['wu' xsub{sub} '_sbref.nii']};

    if inpt == 1
        if exist(fullfile(subjects{sub},datafile{data}),'file')
            cnt = cnt+1;
            if data == 2
                matlabbatch{1}.spm.tools.cat.stools.check_mesh_cov.data_surf{1}{cnt,1} = ...
                    fullfile(subjects{sub},datafile{data});
            else
                matlabbatch{1}.spm.tools.cat.tools.check_cov.data_vol{1}{cnt,1} = ...
                    fullfile(subjects{sub},datafile{data});
            end
    % ----- for cat preprocessing include quality parameters in QC -----
            if data == 1
                matlabbatch{1}.spm.tools.cat.tools.check_cov.data_xml{cnt,1} = ...
                    fullfile(subjects{sub},'report',['cat_' xsub{sub} '.xml']);
            elseif data == 2
                matlabbatch{1}.spm.tools.cat.stools.check_mesh_cov.data_xml{cnt,1} = ...
                    fullfile(subjects{sub},'report',['cat_' xsub{sub} '.xml']);
            end
        else
            fprintf(1,'%s\n',[datafolder{data} ' missing for ' xsub{sub}]);
        end
    elseif inpt == 2
        if exist(fullfile(subjects{sub},datafile{data}),'file')
            if isPat(sub) == 0
                cnt = cnt+1;
                if data == 2
                    matlabbatch{1}.spm.tools.cat.stools.check_mesh_cov.data_surf{1}{cnt,1} = ...
                        fullfile(subjects{sub},datafile{data});
                else
                    matlabbatch{1}.spm.tools.cat.tools.check_cov.data_vol{1}{cnt,1} = ...
                        fullfile(subjects{sub},datafile{data});
                end
    % ----- for cat data include quality parameters in QC -----
                if data == 1
                    matlabbatch{1}.spm.tools.cat.tools.check_cov.data_xml{cnt,1} = ...
                        fullfile(subjects{sub},'report',['cat_' xsub{sub} '.xml']);
                elseif data == 2
                    matlabbatch{1}.spm.tools.cat.stools.check_mesh_cov.data_xml{cnt,1} = ...
                        fullfile(subjects{sub},'report',['cat_' xsub{sub} '.xml']);
                end
                            
    % ----- gather nuisance variables for controls -----
                if cova == 1
                    c(Pcnt+cnt,1) = Cov(sub,find(strcmpi(Covariates,'age')==1));
                elseif cova == 2
                    c(Pcnt+cnt,1) = Cov(sub,find(strcmpi(Covariates,'icv')==1));
                elseif cova == 3
                    c(Pcnt+cnt,1) = Cov(sub,find(strcmpi(Covariates,'gmv')==1));
                elseif cova == 4
                    c(Pcnt+cnt,1) = Cov(sub,find(strcmpi(Covariates,'fd')==1));
                end
            
            elseif isPat(sub) == 1
                Pcnt = Pcnt+1;
                if data == 2
                    matlabbatch{1}.spm.tools.cat.stools.check_mesh_cov.data_surf{2}{Pcnt,1} = ...
                        fullfile(subjects{sub},datafile{data});
                else
                    matlabbatch{1}.spm.tools.cat.tools.check_cov.data_vol{2}{Pcnt,1} = ...
                        fullfile(subjects{sub},datafile{data});
                end
    % ----- for cat data include quality parameters in QC -----
                if data == 1
                    matlabbatch{1}.spm.tools.cat.tools.check_cov.data_xml{Pcnt+cnt,1} = ...
                        fullfile(subjects{sub},'report',['cat_' xsub{sub} '.xml']);
                elseif data == 2
                    matlabbatch{1}.spm.tools.cat.stools.check_mesh_cov.data_xml{Pcnt+cnt,1} = ...
                        fullfile(subjects{sub},'report',['cat_' xsub{sub} '.xml']);
                end               
    % ----- gather nuisance variables for patients -----
                if cova == 1
                    c(Pcnt+cnt,1) = Cov(sub,find(strcmpi(Covariates,'age')==1));
                elseif cova == 2
                    c(Pcnt+cnt,1) = Cov(sub,find(strcmpi(Covariates,'icv')==1));
                elseif cova == 3
                    c(Pcnt+cnt,1) = Cov(sub,find(strcmpi(Covariates,'gmv')==1));
                elseif cova == 4
                    c(Pcnt+cnt,1) = Cov(sub,find(strcmpi(Covariates,'dvars')==1));
                end
            end
        else
            fprintf(1,'%s\n',[datafolder{data} ' missing for ' xsub{sub}]);
        end    
    end
end    

%--------------- with dircectory as input no nuisance variables are considered ----------------------------
if inpt == 1 || cova == 0
    if data == 2
        matlabbatch{1}.spm.tools.cat.stools.check_mesh_cov.c = {};        
    else
        matlabbatch{1}.spm.tools.cat.tools.check_cov.c = {};
    end
    fprintf(1,'\n%s\n',[num2str(cnt) ' of ' num2str(sub) ' subjects found with ' datafolder{data}]);

else
    if data == 2
        matlabbatch{1}.spm.tools.cat.stools.check_mesh_cov.c{1} = c(:);        
    else
        matlabbatch{1}.spm.tools.cat.tools.check_cov.c{1} = c(:);
    end
    fprintf(1,'\n%s\n',[num2str(cnt) ' controls and ' num2str(Pcnt) ' patients found with ' datafolder{data}]);
end


try
    spm_jobman('run',matlabbatch)
catch
    if data == 1
        matlabbatch{1}.spm.tools.cat.tools.check_cov.data_xml = {};
        spm_jobman('run',matlabbatch)
    elseif data == 2
        matlabbatch{1}.spm.tools.cat.stools.check_mesh_cov.data_xml = {};
    end
end
