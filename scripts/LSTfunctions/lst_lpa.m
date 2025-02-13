function lst_lpa(inputFLAIR, report)

tStartBeginScript = tic;

% Initialisation of SPM
spm('defaults','PET');
spm_jobman('initcfg')

% specify LST LPA job
matlabbatch{1}.spm.tools.LST.lpa.data_F2 = {[inputFLAIR, ',1']};
matlabbatch{1}.spm.tools.LST.lpa.data_coreg = {''};
matlabbatch{1}.spm.tools.LST.lpa.html_report = report;

% Specify inputs (already in batch)
% inputs = cell(0, 1);

% Run batch
% spm_jobman('run', matlabbatch, inputs{:});
spm_jobman('run', matlabbatch);

disp(['TotalScriptTime = ' datestr(datenum(0,0,0,0,0,toc(tStartBeginScript)),'HH:MM:SS')])

end

