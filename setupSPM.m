function setupSPM(subject, workdir, results_folder, labels, cwd, TR)



  addpath /data/BnB2/TOOLS/spm12

  work_dir=fullfile(results_folder,subject);


    mkdir(work_dir); delete([work_dir '/*.*']);
    datarun=[];

%     xrun = dir(['Data/' subject '/Orig/*']);
%
%     cnt = 1;
%     runs = {};
%     for i=3:numel(xrun)
%         if xrun(i).isdir
%             runs{cnt} = xrun(i).name;
%             cnt = cnt+1;
%         end
%     end


   % runs = {runs{[1 2 3 4 5 6 7]}};  %  alle 7 Sessions (Reihenfolge: alphabetisch nach Beding.namen)
   % runs = {runs{[2 3 4 5 6 7]}};  %  6 Sessions (ohne erste) (Reihenfolge: alphabetisch nach Beding.namen)
   % runs = {runs{[7]}};  %  alle Sessions einzeln
    runs = labels;
    runs = {runs{[1 2 3 4]}};

    nruns           = numel(runs);
    %for j = 1:nruns; datarun{j}=fullfile(workdir, 'Derivatives',subject,smoothedfolder,runs{j}); end
    %for j = 1:nruns; datarun{j}=fullfile(workdir, 'Derivatives',runs{j},subject,smoothedfolder); end % bei 3-D Niftis
    for j = 1:nruns;
        datarun{j}=fullfile(workdir, 'Derivatives', runs{j},subject);
    end % bei 4-D Nifti



   clear soas
    %load(fullfile('S:\D\Arbeit\ExtKoop\NARPS\data\Behavior\Logfiles\SOAfiles\',[subject '.mat']));
    load(fullfile(workdir, '/Behavior/Logfiles/SOAfiles/', [subject '.mat']));

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

%         SESS{i,1} = Labs;
%         SESS{i,2} = indx;
        SESS{i,1} = Labs{i};
        SESS{i,2} = indx(i);

    end

    if any(isnan([SESS{:,2}]))
        huhu Error!!
    end


    clear nscans

%       P = [];
%     for run = 1:nruns
%         [Files,Dirs] = spm_list_files(datarun{run},'s*.img');
%         P            = [P; [repmat([datarun{run} filesep],size(Files,1),1) Files]];
%         nscans(run)  = size(Files,1);
%     end;

    %size(spm_select('expand', fullfile('S:\D\Arbeit\ExtKoop\NARPS\data\fMRI\MGT1\sub-001\','s8wsub-001.nii')),1)

    P = '';
    for run = 1:nruns
        %files = dir(fullfile(datarun{run},'s*.nii')); % 3-D Niftis
        files = spm_select('expand', fullfile(datarun{run},['s5w' subject '.nii'])); % 4-D Nifti
        %P  = [P; [repmat([datarun{run} filesep], size(files,1), 1) files(i).name]]; % 3-D Niftis
        P  = [P; files];
        nscans(run)  = size(files,1);
    end

    clear SPM; SPM       = struct('nscan',nscans);
    SPM.xY.P  = P;
    SPM.xY.RT = TR;

    SPM.xBF.UNITS  = 'secs';  SPM.xBF.name  = 'hrf (with time derivative)';  % Option: 'hrf'
    SPM.xBF.T      = 16;       SPM.xBF.T0    = 8;       SPM.xBF.dt = SPM.xY.RT/SPM.xBF.T;
    SPM.xBF.length = 32.2;     SPM.xBF.order = 2;       SPM.xBF = spm_get_bf(SPM.xBF);      SPM.xBF.Volterra   = 1;

%     QAfile = dir(fullfile(workdir, 'Derivatives',subject,'Orig','QA.mat'));
%     load([workdir, 'Derivatives', subject, filesep, 'Orig', filesep, QAfile(1).name]);
%     % qaruns = [6 8 4 9 5 10];  % Re-Gruppierung (zu alphabetisch) der Run-Reihenfolge aus der Qualitaetsanylase AQUA
%       (Original-Reihenfolge laut Script aqua_ale:
        % ap_simult00, avp_mix0000, av_simult00, kontr_piezo, piezo000000, audi0000000, avp_simult0, kontr_audi0,...
%         kontr_visu0, visu0000000)
%      qaruns = [1 2 3 4 5 6 7];

    for ses = 1:numel(nscans)

            SPM.Sess(ses).U.ons  = soas(SESS{ses,2}).ons;
            SPM.Sess(ses).U.dur  = soas(SESS{ses,2}).dur;
            SPM.Sess(ses).U.name = {SESS{ses,1}};
            %if numel(soas(SESS{ses,2}(xi)).par2')==0
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

            rp = dir(fullfile(workdir, 'Derivatives',runs{ses},subject,'rp_*.txt'));
            [r1,r2,r3,r4,r5,r6] = textread(fullfile(workdir, 'Derivatives',runs{ses},subject,rp.name),'%f%f%f%f%f%f');


            Conf   = load(fullfile(workdir, 'Derivatives',runs{ses},subject,['Counfounds_' subject '.mat']));
            censTP = zeros(numel(r1),1);
            censTP(Conf.badTP,1) = 1;

            % if sum(censTP) > 0
                %    SPM.Sess(ses).C.C    = [];  % Modell ohne  zusaetzliche Regressoren
                SPM.Sess(ses).C.C    = [r1,r2,r3,r4,r5,r6,censTP];
                %    SPM.Sess(ses).C.name = {};  % Modell ohne zusaetzliche Regressoren
                SPM.Sess(ses).C.name = {'x', 'y', 'z', 'yaw', 'pitch', 'roll', 'badTP'};
            % else
            %     %    SPM.Sess(ses).C.C    = [];  % Modell ohne  zusaetzliche Regressoren
            %     SPM.Sess(ses).C.C    = [r1,r2,r3,r4,r5,r6];
            %     %    SPM.Sess(ses).C.name = {};  % Modell ohne zusaetzliche Regressoren
            %     SPM.Sess(ses).C.name = {'x', 'y', 'z', 'yaw', 'pitch', 'roll'};
            % end
    end

    SPM.xGX.iGXcalc = 'none';
  %  SPM.xVi.form = 'none';
    SPM.xVi.form = 'AR(1) + w';
    SPM.xX.K(1).HParam  = 128;

    cd(work_dir);
    SPM = spm_fmri_spm_ui(SPM);

    save SPM.mat SPM;

    cd(cwd)

end
