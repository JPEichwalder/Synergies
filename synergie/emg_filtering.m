history -c
clear
clc
close all

pkg load signal %signal package is needed for filter functions-->has to be downloaded (pkg install -forge signal)

% add functions path
functionsPath = addpath(genpath('D:\octave\')); %folder with skript and functions

% directory of the folder with the raw c3d and loasol data for each participant
data_folder = ['D:/master_data/raw/'];
data_ev_folder = ['D:/master_data/events/'];

##save_folder = ['D:/studienassisdent/sumo_dl/results/'];

%captuered muscles via sEMG
muscles = {'tib_ant'; 'per_long'; 'soleus'; 'gast_med'; 'vast_lat'; 'rect_fem'; 'bic_fem'; 'sem_tend'; ...
           'add_long'; 'grac';'glut_med'; 'glut_max'; 'rec_abd'; 'ext_obli'; 'mulitfid'; 'erec_spin'};
nr_mus  = length(muscles);

%gets probandnames
probands_1  = readdir([data_folder]);
prb = 1;
for file = 1:length(probands_1)
  switch probands_1{file}(1)
    case 'P'
      probands{prb,1} = probands_1{file};
      prb = prb +1;
  end%switch
end%for
nr_prb      = length(probands);
clear probands_1;

for prb = 1:nr_prb
  proband = probands{prb}
  cd ([data_folder, proband, '/c3d']);

  c3dnames_1  = readdir([data_folder, proband, '/c3d']);%gets names of c3d c3d files
  c3dnames_1  = sort(c3dnames_1);
  nr_c3d      = length(c3dnames_1);

  if nr_c3d > 2
    a = 1;
    for c3 = 1:nr_c3d
      trial = c3dnames_1{c3};
      if sum(strfind(trial, 'MVC')) > 0 || sum(strfind(trial, 'EMG'))>0% || sum(strfind(trial, 'static'))>0
      c3dnames_2{a,1} = c3dnames_1{c3};
      a = a+1;
      end%if
    end%for
    c3dnames_3  = readdir([data_ev_folder, proband]);
    c3dnames_4  = c3dnames_3((3:length(c3dnames_3)),1);
    c3d_names   = vertcat(c3dnames_2, c3dnames_4);

    nr_c3d      = length(c3d_names);
    c3dnames.(probands{prb}) = c3d_names;
  for c3 = 1:nr_c3d
    trial = c3d_names{c3};
    if sum(strfind(trial, 'MVC')) == 0 && sum(strfind(trial, 'EMG'))==0 && sum(strfind(trial, 'static'))==0
      cd ([data_ev_folder, proband]);
      walk = 1;
    else
      cd ([data_folder, proband, '/c3d']);
      filrms = 1;
    end%if
    c3d   = c3d_names{c3};
    file  = ezc3dRead(c3d);
    trial = c3d(1:end-4);
    trialnames.(probands{prb}){c3,1} = trial;
    if exist('walk') && walk == 1
      walkfiles.(proband).(trial) = file;

      [an_fR, an_fF, an_lF]     =analogsInfo(file);%gets basic information for the analogs
      [pts_fR, pts_fF, pts_lF]  =pointsInfo(file); %gets basic information for the points (f.example the markers)
      anmulti = an_fR/pts_fR;
      line1= collectPointsSignals(file, 'lin1');
      line1y = line1(2,1,:);
      nonan = ~isnan(line1y);
      line1y = line1y(nonan>0);
      linm = mean(line1y);
      if linm > 0
        line1= collectPointsSignals(file, 'lin2');
        line1y = line1(2,1,:);
        nonan = ~isnan(line1y);
        line1y = line1y(nonan>0);
        linm = mean(line1y);
      endif
      rtoe = collectPointsSignals(file,'RTOE');
      rtoey = rtoe(2,1,:);
      ltoe = collectPointsSignals(file,'LTOE');
      ltoey = ltoe(2,1,:);
      [xr,yr] = find(rtoey < linm);
      [xl,yl] = find(ltoey < linm);

      marl(1) = max(yr);
      marl(2) = max(yl);

      all_events = getEventTimes(file, 'Foot Strike', 'Right', 2);
      all_events = round(all_events*pts_fR)+1-pts_fF+1;

      firstcyc = all_events(all_events>max(marl))(1:2);
      if length(firstcyc)<2
        disp(['no cycle for trial:', (proband), (trial)]);
      end%if
      an_firstcyc = firstcyc*anmulti;
      cycle_1_an.(proband).(trial) = an_firstcyc;
      clear an_fR an_fF an_lF pts_fR pts_fF pts_lF anmulti line1 line1y linm rtoe rtoey ltoe ltoey xr yr xl yl marl all_events firstcyc an_firstcyc;
    end%if

    trial = c3d(1:end-4);
    trialnames.(probands{prb}){c3,1} = trial;
    if sum(strfind(trial, 'MVC')) == 0 && sum(strfind(trial, 'EMG'))==0 && sum(strfind(trial, 'static'))==0
      walktrialnames.(probands{prb}){c3,1} = trial;

    end%if
    [fR, fF, lF] = analogsInfo(file);

    for mus = 1:nr_mus
      if prb == 1 && mus == 12
        if (strcmp(trial(1:5), 'ebeam')) == 1 || (strcmp(trial(1:5), 'eline')) == 1 || (strcmp(trial(1:5), 'tr_37')) == 1 || (strcmp(trial(1:5), 'tr_39')) == 1 || (strcmp(trial(1:5), 'tr_40')) == 1 || (strcmp(trial(1:5), 'tr_41')) == 1
          EMG_raw1          = collectANALOGSignals(file, 'tens_fas')';
        else
          EMG_raw1          = collectANALOGSignals(file, muscles{mus})';
        end%if
      else
        EMG_raw1          = collectANALOGSignals(file, muscles{mus})'; %gets EMG signal of muscle
      end%if

      EMG_filtered1     = filterSignal_butter(EMG_raw1, 'high', fR, 'order', 4, 'cutoff', 25);
      EMG_demeaned      = EMG_filtered1 - mean(EMG_filtered1); % demeaned
      EMG_rectified     = sqrt(EMG_demeaned.^2); %rectified
      EMG_filtered_low  = filterSignal_butter(EMG_rectified, 'low', fR, 'order', 4, 'cutoff', 7); % 4th order low-pass Butterworth filter 6 Hz
      EMG_filtered_low(EMG_filtered_low<0) = 0;

      if exist('filrms') && filrms == 1
      EMG_filtered_low2 = filterSig_rms(EMG_filtered_low, fR, 0.05);
      EMG_filtered2.(probands{prb}).(trial)(mus,:) = EMG_filtered_low2;
      end%if

      EMG_raw.(probands{prb}).(trial)(mus,:)      = EMG_raw1;
      EMG_filtered.(probands{prb}).(trial)(mus,:) = EMG_filtered_low;
##      slopedifference.(probands{prb}).(trial){mus,1}=difference;

      clear EMG_raw1 EMG_filtered1 EMG_demeande EMG_rectified EMG_filtered_low;
    end%for

  if sum(strfind(trial, 'EMG')) == 0
    EMG_baseclear.(probands{prb}).(trial) = EMG_filtered.(probands{prb}).(trial)- EMG_baseline.(probands{prb}).rm;
    EMG_baseclear.(probands{prb}).(trial)(EMG_baseclear.(probands{prb}).(trial)<0) = 0;
    if sum(strfind(trial, 'MVC')) > 0
      EMG_baseclear2.(probands{prb}).(trial) = EMG_filtered2.(probands{prb}).(trial)- EMG_baseline.(probands{prb}).rm;
      EMG_baseclear2.(probands{prb}).(trial)(EMG_baseclear2.(probands{prb}).(trial)<0) = 0;
    end%if
  end%if

  if exist('walk') && walk == 1
    for mus=1:nr_mus
      muscle(mus,:) = EMG_baseclear.(probands{prb}).(trial)(mus,:);
      mvc = MVC_normalizer2.(probands{prb}).values(mus,1);
      amplNorm  = 1/mvc;
      muscle_Norm(mus,:)=muscle(mus,:)*amplNorm;
      bigger=find(muscle_Norm>1);
##      if length(bigger)>0
##        disp(['bigger values than normalized in: ', probands{prb}, trial]);
##      end%if
      if length(cycle_1_an.(probands{prb}).(trial))==2
        timenorm(mus,:)  = timeNormalizeSig(muscle_Norm(mus,cycle_1_an.(probands{prb}).(trial)(1):cycle_1_an.(probands{prb}).(trial)(2)), 101);
      end%if
    end%for
    EMG_amplnorm.(probands{prb}).(trial) = muscle_Norm;
    EMG_timenorm.(probands{prb}).(trial) = timenorm;
    clear muscle max_value amplNorm muscle_Norm mvc timenorm;
  end%if
  if c3 == 1
  MVC_normalizer.(probands{prb}).values(16,1) = 0;
  MVC_normalizer2.(probands{prb}).values(16,1) = 0;
  MVC_end.(probands{prb}).values(1:16,1) = 0;
  MVC_end2.(probands{prb}).values(1:16,1) = 0;
  end%if
  switch trial(1:4)
    case 'EMG_'
      if prb == 4%Prb5
        EMG_baseline.(probands{prb}).me = mean(EMG_filtered.(probands{prb}).(trial)(:,3000:end-500),2);
        EMG_baseline.(probands{prb}).rm = rms(EMG_filtered.(probands{prb}).(trial)(:,3000:end-500),2);
      elseif prb == 5%prb6
        EMG_baseline.(probands{prb}).me = mean(EMG_filtered.(probands{prb}).(trial)(:,500:end-11000),2);
        EMG_baseline.(probands{prb}).rm = rms(EMG_filtered.(probands{prb}).(trial)(:,500:end-11000),2);
      elseif prb == 7%prb8
        EMG_baseline.(probands{prb}).me = mean(EMG_filtered.(probands{prb}).(trial)(:,5000:end-500),2);
        EMG_baseline.(probands{prb}).rm = rms(EMG_filtered.(probands{prb}).(trial)(:,5000:end-500),2);
      else
        EMG_baseline.(probands{prb}).me = mean(EMG_filtered.(probands{prb}).(trial)(:,1000:end-1000),2);
        EMG_baseline.(probands{prb}).rm = rms(EMG_filtered.(probands{prb}).(trial)(:,1000:end-1000),2);
      end%if
      if max(EMG_baseline.(probands{prb}).rm) > 0.01
        display(['EMG baseline too big for ', probands{prb}])
      end%if
    case 'MVC_'
      maxval = max(EMG_baseclear.(probands{prb}).(trial)(:,500:end-500),[],2);
      if (strfind(trial, 'abdu'))>0
        mvcmuscles = [11,12];
      elseif (strfind(trial, 'addu'))>0
        mvcmuscles = [7,8,9,10];
      elseif (strfind(trial, 'bauc'))>0
        mvcmuscles = [6,13,14];
      elseif  (strfind(trial, 'dors'))>0
        mvcmuscles = [1,2];
      elseif (strfind(trial, 'glut'))>0
        mvcmuscles = [11,12,15,16];
      elseif (strfind(trial, 'knieex'))>0
        mvcmuscles = [5,6];
      elseif (strfind(trial, 'kniefl'))>0
        mvcmuslces = [4,7,8,9,10];
      elseif (strfind(trial, 'plan'))>0
        mvcmuscles = [2,3,4];
      elseif (strfind(trial, 'ruec'))>0
        mvcmuscles = [7,8,11,12,15,16];
      end%if
      for mus = 1:nr_mus
        if maxval(mus)> MVC_normalizer.(probands{prb}).values(mus)
          MVC_normalizer.(probands{prb}).values(mus)  = maxval(mus);
          MVC_normalizer.(probands{prb}).trial{mus,1}   = trial;
        end%if
      end%for
      maxval2 = max(EMG_baseclear2.(probands{prb}).(trial)(:,500:end-500),[],2);
      for mus = 1:nr_mus
        if maxval2(mus)> MVC_normalizer2.(probands{prb}).values(mus)
          MVC_normalizer2.(probands{prb}).values(mus)  = maxval2(mus);
          MVC_normalizer2.(probands{prb}).trial{mus,1}   = trial;
        end%if
      end%for
    case 'MVCe'
      maxval = max(EMG_baseclear.(probands{prb}).(trial)(:,500:end-500),[],2);
      for mus = 1:nr_mus
        if maxval(mus)> MVC_end.(probands{prb}).values(mus)
          MVC_end.(probands{prb}).values(mus)  = maxval(mus);
          MVC_end.(probands{prb}).trial{mus,1}   = trial;
        end%if
      end%for
      maxval2 = max(EMG_baseclear2.(probands{prb}).(trial)(:,500:end-500),[],2);
      for mus = 1:nr_mus
        if maxval(mus)> MVC_end2.(probands{prb}).values(mus)
          MVC_end2.(probands{prb}).values(mus)  = maxval2(mus);
          MVC_end2.(probands{prb}).trial{mus,1}   = trial;
        end%if
      end%for

  end%switch
  clear c3d file fR fF lF walk filrms;
  end%for
  clear c3dnames_1 c3dnames_2 c3dnames_3 c3dnames_4 c3d_names;

  end%if
end%for
cd (data_ev_folder);
cd filter_try

save -mat-binary 'EMG_MVCc.mat' 'MVC_end' 'MVC_end2' 'MVC_normalizer' 'MVC_normalizer2' 'EMG_baseline';
save -mat-binary 'EMG_rawc.mat' 'EMG_raw'
save -mat-binary 'EMG_filteredc.mat' 'EMG_filtered' 'EMG_filtered2' 'EMG_baseclear' 'EMG_baseclear2';
save -mat-binary 'EMG_walkfilesc.mat' 'walkfiles';
save -mat-binary 'EMG_normc.mat' 'EMG_timenorm' 'EMG_amplnorm';
save -mat-binary 'eventsc.mat' 'cycle_1_an';
##abprobands = fieldnames(walkfiles);
##for prb = 1:length(abprobands);
##
