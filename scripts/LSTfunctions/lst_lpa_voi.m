function lst_lpa_voi(input_ples, th)

tStartBeginScript = tic;

th = str2double(th);

% Initialisation of SPM
spm('defaults','PET');
spm_jobman('initcfg');

% Extract values of interest
matlabbatch{1}.spm.tools.LST.tlv.data_lm = {[input_ples, ',1']};
matlabbatch{1}.spm.tools.LST.tlv.bin_thresh = th;

% Specify inputs (already in batch)
inputs = cell(0, 1);

% Run batch
spm_jobman('run', matlabbatch, inputs{:});

disp(['TotalScriptTime = ' datestr(datenum(0,0,0,0,0,toc(tStartBeginScript)),'HH:MM:SS')])

end
