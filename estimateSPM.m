function estimateSPM(work_dir)

    addpath /data/BnB2/TOOLS/spm12

	load(fullfile(work_dir, 'SPM.mat'))
    SPM.swd = work_dir;

    SPM = spm_spm(SPM);

    SPM = rmfield(SPM,'xCon');
    cnt = 1;
    for ses = 1:numel(SPM.Sess)
        for con = 1:numel(SPM.Sess(ses).Fc)

            c = zeros(1,size(SPM.xX.X,2));
            c(SPM.Sess(ses).col(SPM.Sess(ses).Fc(con).i(1))) = 1;

            if cnt == 1
                SPM.xCon           = spm_FcUtil('Set',[SPM.Sess(ses).Fc(con).name],'T','c',c',SPM.xX.xKXs);
            else
                SPM.xCon(end+1)    = spm_FcUtil('Set',[SPM.Sess(ses).Fc(con).name],'T','c',c',SPM.xX.xKXs);
            end

            cnt = cnt+1;
        end
    end

    spm_contrasts(SPM);

    % Sum contrasts for regressors of interest
    c = zeros(4,64);
    c(1,1) = 1; c(1,16) = 1; c(1,31) = 1; c(1,46) = 1; % Task (main)
    c(2,3) = 1; c(2,18) = 1; c(2,33) = 1; c(2,48) = 1; % Gain (PM 1)
    c(3,5) = 1; c(3,20) = 1; c(3,35) = 1; c(3,50) = 1; % Loss (PM 2)
    c(4,7) = 1; c(4,22) = 1; c(4,37) = 1; c(4,52) = 1; % RT (PM 3)

    % c(1,1) = 1; c(1,15) = 1; c(1,29) = 1; c(1,43) = 1; % Task (main)
    % c(2,3) = 1; c(2,17) = 1; c(2,31) = 1; c(2,45) = 1; % Gain (PM 1)
    % c(3,5) = 1; c(3,19) = 1; c(3,33) = 1; c(3,47) = 1; % Loss (PM 2)
    % c(4,7) = 1; c(4,21) = 1; c(4,35) = 1; c(4,49) = 1; % RT (PM 3)
    SPM.xCon(5)    = spm_FcUtil('Set','Task','T','c',c(1,:)',SPM.xX.xKXs);
    SPM.xCon(6)    = spm_FcUtil('Set','Gain','T','c',c(2,:)',SPM.xX.xKXs);
    SPM.xCon(7)    = spm_FcUtil('Set','Loss','T','c',c(3,:)',SPM.xX.xKXs);
    SPM.xCon(8)    = spm_FcUtil('Set','RT',  'T','c',c(4,:)',SPM.xX.xKXs);

    spm_contrasts(SPM,[5:8]);

end
