function save_lcmv_spectral6_indiv(sbjs, data_dir, results_dir)

load metadata

fs = 128;
fres = 2*fs;

frqs = sfreqs(fres, fs);
maxfreq_ind = max(find(frqs <= 45));
frqs = frqs(1:maxfreq_ind);


[b_notch, a_notch]=butter(3,[49 51]/fs*2,'stop'); % 60Hz line noise
[hpnum,hpdenom]=butter(3, 1/fs*2,'high'); % drift removal
a_notch = poly([roots(hpdenom);roots(a_notch)]);
b_notch = conv(hpnum,b_notch);

eloc128 = readlocs('biosemi-layout-128.loc');
clab = {eloc128.labels};

M = length(clab);
H = eye(M)-ones(M) ./ M;

sa = prepare_sourceanalysis(clab, 'icbm152b_sym_biosemi128');
% load in_cort
L = tprod(H, [1 -1], sa.cortex75K.V_bem(:, sa.cortex1K.in_from_cortex75K(sa.cortex1K.in_cort), :), [-1 2 3]);
nvox = size(L, 2);
clear sa

for isbj = intersect(1:length(fnames), sbjs)
  isbj
  load([data_dir fnames{isbj}]);

  data = resample(double(EEG.data)', fs, EEG.srate);
  data = filter(b_notch, a_notch, data);

%   data = data(10*fs + (1:60*fs*2), :);
%   data = data(130*fs + (1:60*fs*2), :);
  data = data(10*fs + (1:60*fs*4), :);
  
  data = tprod(data, [1 -1 3], H, [-1 2]);
  
  [power_lcmv_1D, P_lcmv_1D, data_lcmv_1D, info_lcmv] = lcmv(data, L);


  [TRGC, GC, PS, ~, nlags] = data2strgc(data_lcmv_1D, fres*2, [], maxfreq_ind, 0, []);

  
%   semilogy(frqs, mean(PS, 3)')

  mkdir([results_dir 'lcmv_spectral6_indiv'])
  save([results_dir 'lcmv_spectral6_indiv/sbj' num2str(isbj)], ...
    'fs', 'fres', 'frqs', 'TRGC', 'GC', 'PS', 'nlags', 'P_lcmv_1D', 'maxfreq_ind', '-v7.3');
end
