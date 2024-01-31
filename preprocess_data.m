function preprocess_data(sbjs, data_dir, results_dir)

addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Tutorials/BuenosAires/'));
addpath('/Users/nguyentiendung/Desktop/Studium/Charite/matlab/eeglab')

data_dir = '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Tutorials/BuenosAires/data/EEG/';
results_dir = '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Tutorials/BuenosAires_project/results/';

%%% Creating result folder (Plots will be saved here)
resultDir = ['../analysis_output/preprocessing/'];
if ~isdir(resultDir)
    mkdir(resultDir)
end


%%% Creating output data (matfile) folder
outputDataDir = [resultDir, 'data/'];
if ~isdir(outputDataDir)
    mkdir(outputDataDir)
end
    
load metadata

fs = 128;
fres = fs;

frqs = sfreqs(fres, fs);
% maxfreq_ind = max(find(frqs <= 45));
% frqs = frqs(1:maxfreq_ind);

%% Preprocessing Settings


saveOutputData = 1;   % 1: To save preprocessd data and preprocessing details  0: Not save

%% %%%%%%%%%%%%%%%%%
lpFilter =   63;       % low-pass filter cut-off
hpFilter =   1;      % high-pass filter cut-off
bsFilter =   [49 51];       % band-stop filter range
filterOrder = 2;    % Butterworth filter order
dsRatio =  4;       % downsampling rate
epochLeng = 2;   % epoch length in seconds
mincomp = 30;

% Updating 'info' variable
c = fix(clock);
info = [];
info.prepStartTime = [' ',num2str(c(3)),'.',num2str(c(2)),'.',num2str(c(1)),'     ',num2str(c(4)),':',num2str(c(5))];
info.prep.lpFilter = lpFilter;
info.prep.hpFilter  = hpFilter;
info.prep.filterOrder = filterOrder;


[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

for isbj = intersect(1:length(fnames), sbjs)
  isbj
  % set isbj = 14;
  load([data_dir fnames{isbj}]);
  
  EEG.epoch = 1;


  % Updating "info" variable
  info.prep.Nchan0 = EEG.nbchan; info.data.fs0 = EEG.srate;
  info.data.lengthPoints = EEG.pnts;
  info.prep.ref0 = 'common'; info.prep.chanLabels = {EEG.chanlocs.labels};



  %% Plotting Channel locations
  EEG = eeg_checkset( EEG );
  figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
  saveas(gcf,[resultDir,'topo_labels'],'jpg');
  EEG = eeg_checkset( EEG );
  figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'numpoint', 'chaninfo', EEG.chaninfo);
  saveas(gcf,[resultDir,'topo_indices'],'jpg');


  %% Plotting Spectra of raw data
  EEG = eeg_checkset( EEG );
  figure; pop_spectopo(EEG, 1, [0 EEG.xmax*1000], 'EEG' , 'percent', 15, 'freq', [10 20], 'freqrange',[0 63],'electrodes','on');
  saveas(gcf,[resultDir,'spectra_raw'],'jpg');


  %% zero-phase filtering
  [b_low, a_low] = butter(filterOrder, lpFilter/(EEG.srate/2), 'low');
  [b_high, a_high] = butter(filterOrder, hpFilter/(EEG.srate/2), 'high');
  [b_notch, a_notch] = butter(filterOrder, bsFilter/(EEG.srate/2),'stop');
  a_all = poly([roots(a_low);roots(a_high); roots(a_notch)]);
  b_all = conv(conv(b_low,b_high), b_notch);
  %     fvtool(b_all, a_all)
  EEG.data = filtfilt(b_all, a_all, double(EEG.data)')';


  %% plain downsampling to 100 Hz. 
  % that's not 100 Hz but 128 Hz
  EEG.data = EEG.data(:, 1:dsRatio:end);
  EEG.srate = EEG.srate/dsRatio;
  EEG.pnts    = size(EEG.data,2);
  EEG.xmax    = EEG.xmin + (EEG.pnts-1)/EEG.srate;
  EEG.times   = linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts);
  EEG = eeg_checkset( EEG );


  %% Epoching: Trasforming data from Nchan*NtpAll  to Nchan*NtpEpoch*Nepoch
  %epochLeng = 2; % epoch length in seconds
  [Nchan, Ntp]=size(EEG.data);
  EEG = pop_select( EEG,'point',[1 (Ntp - mod(Ntp,EEG.srate*epochLeng))] );
  Nepoch = EEG.pnts / (EEG.srate*epochLeng);
  EEG = eeg_checkset(EEG);
  %EEG.data = EEG.data(:, 1:(Ntp - mod(Ntp,fsNew*8)) );
  for ievent = 1: Nepoch
      EEG.event(ievent).type = num2str(epochLeng);
      EEG.event(ievent).latency = (ievent-1)*(EEG.srate*epochLeng)+1;
      EEG.event(ievent).duration = epochLeng;
  end
  EEG = eeg_checkset( EEG );
  EEG = pop_epoch( EEG, {  num2str(epochLeng)  }, [0  epochLeng], 'epochinfo', 'yes');


  %% Plotting Spectra of filtered and downsampled data
  EEG = eeg_checkset( EEG );
  figure; pop_spectopo(EEG, 1, [0      EEG.xmax*1000*Nepoch], 'EEG' , 'percent', 15, 'freq', [10 20 ], 'freqrange',[0 63],'electrodes','on');
  saveas(gcf,[resultDir,'spectra_prep1'],'jpg');


  % Updating 'info' variable
  % info.subject.ID = isbj;
  % info.subject.eyeCondition = EYE;
  % info.subject.medicalCondition = MED;
  % info.subject.fileDirectory = dataDir;
  info.subject.outputDirectory = resultDir;
  info.prep.fs = EEG.srate;  % sampling frequency after down-sampling
  info.prep.Nepochs0 = size(EEG.data,3);
  info.prep.epochPoints = size(EEG.data, 2);
  info.prep.Nchannels0 = size(EEG.data, 1);
  info.prep.dsRation = dsRatio;
  info.prep.ref = '';


  %% Detecting Bad channels/epochs
  par.HFAnalysis.run = 1 ;
  par.HFAnalysis.channelRejection = 1;
  info.outlier = k1_detect_bad_epoch_channel(EEG,par);

  %% Detecting channels/epochs rejected because of Strong alpha activity
  [b, a] = butter(2, [5 12]/(EEG.srate/2),'stop');
  XNoAlpha = reshape(filtfilt(b, a, reshape(double(EEG.data), Nchan, [])')',Nchan,[],Nepoch);
  par = [];
  par.srate = EEG.srate;
  par.deviationAnalysis.run = 1;
  par.HFAnalysis.run = 1;
  par.deviationAnalysis.threshold = 4.5;
  par.HFAnalysis.channelRejection =1;
  infoAlpha = k1_detect_bad_epoch_channel(XNoAlpha, par);

  info.outlier.epochRejectFinal = info.outlier.epochRejectFinal & infoAlpha.epochRejectFinal;  % |



  %% Plotting Rejceted channels/epochs
  epoch = find(info.outlier.epochRejectFinal);
  epochchanind = cell(1, length(epoch)); %{[],[3,9,12],[6:12,16,18]};
  rejepochcol =  [.6 .8 .9];
  rejepoch = zeros(1, EEG.trials);
  rejepoch(epoch) = ones(1, length(epoch));
  rejepochE = zeros(EEG.nbchan, EEG.trials);
  for i = 1:length(find(rejepoch))
      rejepochE(epochchanind{i},epoch(i))=ones(size(epochchanind{i}));
  end
  winrej = trial2eegplot(rejepoch, rejepochE, EEG.pnts, rejepochcol);
  colors = cell(1,Nchan); colors(1,:) = {'k'};
  colors(1,info.outlier.chanRejectFinal) = {'r'};
  eegplot(EEG.data, 'srate', EEG.srate, 'title', 'Rejected Channels Epochs', ...
      'limits', [EEG.xmin EEG.xmax]*1000, 'color', colors, 'winrej',winrej, 'eloc_file', EEG.chanlocs);


  %% Additional plots for testing channel/epoch rejection
  if 1
      figure, imagesc(info.outlier.epoch.epochDeviation);
      title('Epoch-Epoch-Deviation')
      xlabel('Epoch','fontweight','b'); ylabel('Channel','fontweight','b');
  %         xticks(5:5:Nepoch);
  %         yticks(1:2:Nchan);
      saveas(gcf,[resultDir,'epoch_epoch_deviation'],'jpg');

      figure, imagesc(info.outlier.epoch.chanDeviation);
      title('Epoch-Chan-Deviation')
      xlabel('Epoch','fontweight','b'); ylabel('Channel','fontweight','b');
  %         xticks(5:5:Nepoch);
  %         yticks(1:2:Nchan);
      saveas(gcf,[resultDir,'epoch_chan_deviation'],'jpg');

      figure, imagesc(info.outlier.epoch.epochHFNoise);
      title('Epoch-Epoch-HFNOISE')
      xlabel('Epoch','fontweight','b'); ylabel('Channel','fontweight','b');
  %         xticks(5:5:Nepoch);
  %         yticks(1:2:Nchan);
      saveas(gcf,[resultDir,'epoch_epoch_hfnoise'],'jpg');

      figure, imagesc(info.outlier.epoch.chanHFNoise);
      title('Epoch-Chan-HFNOISE')
      xlabel('Epoch','fontweight','b'); ylabel('Channel','fontweight','b');
  %         xticks(5:5:Nepoch);
  %         yticks(1:2:Nchan);
      saveas(gcf,[resultDir,'epoch_chan_hfnoise'],'jpg');

      chanEpochMatrix = zeros(Nchan, Nepoch);
      chanEpochMatrix(info.outlier.chanRejectFinal, :) = 1;
      chanEpochMatrix(:, info.outlier.epochRejectFinal) = 1;
      figure, imagesc(chanEpochMatrix);
      title('Detected outlier channels epochs')
      xlabel('Epoch','fontweight','b'); ylabel('Channel','fontweight','b');
  %         xticks(5:5:Nepoch);
  %         yticks(1:2:Nchan);
      saveas(gcf,[resultDir,'outlier_chans_epochs'],'jpg');
  end



  %% Rejecting bad channels/epochs before ICA
  info.originalChanlocs = EEG.chanlocs;
  EEG = pop_select( EEG,'notrial',find(info.outlier.epochRejectFinal) ,'nochannel',find(info.outlier.chanRejectFinal));


  %% Computing ICA weights accoding to INFOMAX implemented in EEGLAB
  npca = sum(max(svd(EEG.data(:, :)))./sort(svd(EEG.data(:, :)), 'descend') < 1e4);
  
  if npca < mincomp
    warning(['Effective rank lower than ' num2str(mincomp) '. Lowering minum number of components to ' num2str(npca) '.'])
  end
  
  EEG = pop_runica(EEG, 'extended',1,'interupt','off', 'pca', npca);
  %TO DO %[weights,sphere,meanvar,bias,signs,lrates,data,y] = runica(EEG.data);


  %% find subset of electrodes that is close to those defined by 10/20 system on which MARA is trained
  load electrode_montages_MARA
  
  d = eucl([EasyCap128.Channel.Loc]', [BioSemi128.Channel.Loc]');
  [di mi] = min(d, [], 2);
  [so inso] = sort(di);
  nmatch = sum(so*1000 < 8);

  
  clab_BioSemi128 = {BioSemi128.Channel(mi(inso(1:nmatch))).Name};
  clab_EasyCap128 = {EasyCap128.Channel(inso(1:nmatch)).Name};
  
  [clabs_int, ia, ib] = intersect(clabs, clab_BioSemi128);
  
  chanlocs_bak = EEG.chanlocs;
  
  for ii = 1:length(EEG.chanlocs)
    EEG.chanlocs(ii).labels = ['X' num2str(ii)];
  end
  
  for ii = 1:length(clabs_int)
    EEG.chanlocs(ia(ii)).labels = clab_EasyCap128{ib(ii)};
  end
  
  %% Detecting bad ICs based on MARA
  [ALLEEG,EEG,CURRENTSET] = processMARA (ALLEEG,EEG,CURRENTSET );
  
  EEG.chanlocs = chanlocs_bak;

  %% Ensuring that at least mincomp ICs are retained (According to MARA artefact probabiliy)
  % unless npca is even lower
  if sum(~EEG.reject.gcompreject) < min(mincomp, npca)
      [~, indxICascending] = sort(EEG.reject.MARAinfo.posterior_artefactprob);
      EEG.reject.gcompreject(indxICascending(1:min(mincomp, npca))) = 0;
  end
  
  
  
  %% Visualizing ICs and saving plots
  %pop_visualizeMARAfeatures_modified(EEG.reject.gcompreject, EEG.reject.MARAinfo);
  EEG.saveDir = resultDir;
%   pop_selectcomps_MARA_save(EEG);

  % Updating 'info' variable
  info.rejectedIC = EEG.reject.gcompreject;
  c = fix(clock);
  info.prepEndTime = [' ',num2str(c(3)),'.',num2str(c(2)),'.',num2str(c(1)),'   ',num2str(c(4)),':',num2str(c(5))];


  %% Rejecting bad ICs and back-projecting the retained ICs to Sensor space
  [EEG, com] = pop_subcomp(EEG,find(EEG.reject.gcompreject));


  %% Plotting Spectra after ICA analysis
  EEG = eeg_checkset( EEG );
  figure; pop_spectopo(EEG, 1, [0      EEG.xmax*1000*Nepoch], 'EEG' , 'percent', 15, 'freq', [10 20 ], 'freqrange',[0 63],'electrodes','on');
  saveas(gcf,[resultDir,'spectra_prep2'],'jpg');


  %% Detecting Bad epochs, again, after ICA
  par = [];
  par.HFAnalysis.run = 1 ;
  par.HFAnalysis.channelRejection = 0;
  info.outlier2 = k1_detect_bad_epoch_channel(EEG,par);

  %% Plotting Detected Bad epochs after ICA
  chanEpochMatrix = zeros(Nchan, Nepoch);
  chanEpochMatrix(info.outlier2.chanRejectFinal, :) = 1;
  chanEpochMatrix(:, info.outlier2.epochRejectFinal) = 1;
  figure, imagesc(chanEpochMatrix);
  title('Detected outlier channels epochs after ICA')
  xlabel('Epoch','fontweight','b'); ylabel('Channel','fontweight','b');
  %     xticks(5:5:Nepoch);
  %     yticks(1:2:Nchan);
  saveas(gcf,[resultDir,'outlier2_chans_epochs'],'jpg');

  %% Rejecting bad epochs found after ICA
  %     info.originalChanlocs = EEG.chanlocs;
  EEG = pop_select( EEG,'notrial',find(info.outlier2.epochRejectFinal) );


  %% Saving data before channel interpolation
  if saveOutputData
      % in EEGLAB format
      EEG = eeg_checkset( EEG );
      pop_saveset( EEG, 'filename',['prep_sub_',num2str(isbj),'.set'],'filepath',outputDataDir);
      % in .mat format
      %save(['prep_sub_',num2str(isbj)], 'EEG')
  end

  %% Interpolating rejected channels
  EEG = pop_interp(EEG, info.originalChanlocs, 'spherical');

  %% Plotting Spectra after Interpolation
  EEG = eeg_checkset( EEG );
  figure; pop_spectopo(EEG, 1, [0      EEG.xmax*1000*Nepoch], 'EEG' , 'percent', 15, 'freq', [10 20], 'freqrange',[0 63],'electrodes','on');
  saveas(gcf,[resultDir,'spectra_prep3'],'jpg');

  %% Saving Interpolated data
  if saveOutputData
      % in EEGLAB format
      EEG = eeg_checkset( EEG );
      pop_saveset( EEG, 'filename',['interp_prep_sub_',num2str(isbj),'.set'],'filepath',outputDataDir);
      % in .mat format
      %save(['prep_sub_',num2str(isbj)], 'EEG')
  end

  %% Saving preprocessing details
  if saveOutputData
      save([outputDataDir,'info_preprocessing'], 'info')
  end


  %% Making report

  k1_mk_reprot_latex(info,resultDir);


  %% Converting tex report file to PDF
  %     command = ['/usr/local/texlive/2017/bin/x86_64-darwin/pdflatex  ',...
  %         '-output-directory '  resultDir, ' ',resultDir 'report.tex;'];
  command = ['cd ', resultDir,'; ', '/usr/local/texlive/2017/bin/x86_64-darwin/pdflatex  ',...
      'report.tex;'];
  system(command);

  %% Opening Report file on MAC
  system(['open -a Preview ',resultDir,'report.pdf'])
    
    

%   mkdir([results_dir 'sensorspace_spectral6'])
%   save([results_dir 'sensorspace_spectral6/sbj' num2str(isbj)], ...
%     'fs', 'fres', 'frqs', 'TRGC', 'GC', 'PS', 'CS', 'nlags', ...
%     'maxfreq_ind', '-v7.3');
end
