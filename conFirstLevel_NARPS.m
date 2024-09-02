%spm
clear all

global defaults,
spm('defaults','FMRI');
%spm_defaults


cwd = pwd;


MODELL = 'NARPS_gain1st';  


statpath = (fullfile(cwd,'SingleSubjectAnalysis', MODELL))

subs = spm_select(Inf,'dir','Select subjects','',[statpath filesep]);    
for i=1:size(subs,1)
    subjects{i} = subs(i,end-6:end);
end



for sub= 1:numel(subjects);
    
    cd(fullfile(statpath,subjects{sub}))
    try
    load('SPM.mat')
    
    c = zeros(4,60);
    c(1,1) = 1; c(1,15) = 1; c(1,29) = 1; c(1,43) = 1; % Task (main)
    c(2,3) = 1; c(2,17) = 1; c(2,31) = 1; c(2,45) = 1; % Gain (PM 1)
    c(3,5) = 1; c(3,19) = 1; c(3,33) = 1; c(3,47) = 1; % Loss (PM 2)
    c(4,7) = 1; c(4,21) = 1; c(4,35) = 1; c(4,49) = 1; % RT (PM 3)
    SPM.xCon(5)    = spm_FcUtil('Set','Task','T','c',c(1,:)',SPM.xX.xKXs);
    SPM.xCon(6)    = spm_FcUtil('Set','Gain','T','c',c(2,:)',SPM.xX.xKXs);
    SPM.xCon(7)    = spm_FcUtil('Set','Loss','T','c',c(3,:)',SPM.xX.xKXs);
    SPM.xCon(8)    = spm_FcUtil('Set','RT','T','c',c(4,:)',SPM.xX.xKXs);
    
    spm_contrasts(SPM,[5:8]);
    
    end
    continue

end

cd(cwd)