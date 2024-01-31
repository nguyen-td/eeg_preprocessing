load metadata

data_dir = '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Tutorials/BuenosAires/data/EEG/';
results_dir = '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Tutorials/BuenosAires_project/results/';

% collect_sensorspace_spectral(1:51, data_dir, results_dir);

% univ_sensorspace_spectral(1:51, data_dir, results_dir);


% for isbj = 1:nfiles
%    mgsub({}, @save_sensorspace_spectral5, {isbj, data_dir, results_dir}, 'qsub_opts', '-l h_vmem=25G');%, 'toolboxes', {'statistics_toolbox'}
% end

% for isbj = 1:nfiles
%    mgsub({}, @save_eloreta_spectral6, {isbj, data_dir, results_dir}, 'qsub_opts', '-l h_vmem=25G');%, 'toolboxes', {'statistics_toolbox'}
% end
% 
% for isbj = 1:nfiles
%    mgsub({}, @save_lcmv_spectral6, {isbj, data_dir, results_dir}, 'qsub_opts', '-l h_vmem=25G');%, 'toolboxes', {'statistics_toolbox'}
% end

for isbj = 1:nfiles
   mgsub({}, @save_eloreta_spectral6_indiv, {isbj, data_dir, results_dir}, 'qsub_opts', '-l h_vmem=25G');%, 'toolboxes', {'statistics_toolbox'}
end

for isbj = 1:nfiles
   mgsub({}, @save_lcmv_spectral6_indiv, {isbj, data_dir, results_dir}, 'qsub_opts', '-l h_vmem=25G');%, 'toolboxes', {'statistics_toolbox'}
end
