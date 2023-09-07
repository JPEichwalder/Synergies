history -c  %clears command history
clear         %clears workspace
clc           %clears command window
close all     %close all figures

pkg load signal %signal package is needed for filter functions-->has to be downloaded (pkg install -forge signal)

% add functions path
functionsPath = addpath(genpath('D:\octave\')); % replace with your folderpath with scripts and functions, can contain subfolders with scripts ('\' is important at the end of the path to read subfolders); have all your octave scripts and functions in this folder

% directory of c3d files and where to save the data --> replace with your path (always '\'at the ending)
data_folder = ['D:\synergies_constraints_july\c3d\'];
save_folder = ['D:\synergies_constraints_july\results'];

%conditions that should be analyzed
conditions  = {'Norm'; 'Met_indiv'; 'Bod_indiv'; 'Kombi_indiv'; 'Einh_kombi'}; % changing this names = change names in the script
nr_con      = length(conditions);

% muscles via sEMG that should be analyzed in the order they should be shown (I recommend the way its now) --> name of sensors in vicon --> send me list with differneces in namings and muscles so I can change the script for these
muscles = {'L_tib_ant'; 'L_sol'; 'L_gast_med'; 'L_vast_lat'; 'L_semitend'; 'L_bic_fem'; 'L_glut_med'; ...
           'R_tib_ant'; 'R_sol'; 'R_gast_med'; 'R_vast_lat'; 'R_semitend'; 'R_bic_fem'; 'R_glut_med'};
nr_mus  = length(muscles);

sensors_p_07_session_1 = {'4_tib_ant';  '16_unused'; '14_unused'; '9_unused'; 'Trigger';      '8_rect_fem'; '2_gast_med'; ...
                          '3_soleus';   '15_unused'; '13_unused'; '5_per_long'; '11_semitend';  '7_vast_med'; 'G_gast_lat'};
                          
sensors_p_03_session_2 = {'4_R_fingFlex'; 'Fy_right'; 'Fy_left'; '6_L_tricep'; 'not_used'; '8_R_bicep'; '2_L_fingFlex'; ...
                          '3_R_fingEx'; 'Fz_right'; 'Fz_left'; '5_L_bicep'; '11_R_lat'; '7_L_lat'; '1_L_fingEx';};
                          
% joints that should be analyzed
jointsc3d   = {'LAnkleAngles_M';                              'LKneeAngles_M';                            'LHipAngles_M';                           'LPelvisAngle';...
              'RAnkleAngles_M';                               'RKneeAngles_M';                            'RHipAngles_M';                           'RPelvisAngle'};
jointnames  = {{'L_ankle_flex'; 'L_ankle_ad'; 'L_ankle_rot'}; {'L_knee_flex'; 'L_knee_ad'; 'L_knee_rot'}; {'L_hip_flex'; 'L_hip_ad'; 'L_hip_rot'};  {'L_pelvis_tilt'; 'L_pelvis_obli'; 'L_pelvis_rot'};...
               {'R_ankle_flex'; 'R_ankle_ad'; 'R_ankle_rot'}; {'R_knee_flex'; 'R_knee_ad'; 'R_knee_rot'}; {'R_hip_flex'; 'R_hip_ad'; 'R_hip_rot'};  {'R_pelvis_tilt'; 'R_pelvis_obli'; 'R_pelvis_rot'}};
nr_jnt = length(jointsc3d);
h_co = 20;  %cutoff frequency for high pass filter (I recommend 20 - 35)
l_co = 9;   %cutoff frequency for low pass filter (I recommend 7 - 14 or duration dependent)

%gets probandnames
probands_1  = readdir([data_folder]);
prb = 1;
for file = 1:length(probands_1)
  switch probands_1{file}(1)
    case 'p' % first letter of the probands folders (no other folder in 'data_folder' should start with this)
      probands{prb,1} = probands_1{file};
      prb = prb +1;
  end%switch
end%for -->file
nr_prb      = length(probands);
clear probands_1 prb file;

for prb = 1:nr_prb
  proband = probands{prb}
  
  clear session;
  session_1  = readdir([data_folder, proband]);
  ses = 1;
  for file = 1:length(session_1)
    switch session_1{file}(1)
      case 'S' % first letter of the probands folders (no other folder in 'data_folder' should start with this) --> change something in the session names = have to change 'Session_...' in script
        sessions{ses,1} = session_1{file};
        ses = ses +1;
    end%switch
  end%for --> file
  nr_ses      = length(sessions);
  clear session_1 ses file;
  
  for ses = 1:nr_ses
    session = sessions{ses}
        
    cd ([data_folder, proband, '\', session]); % cd = change to specific folder
    c3dnames_1  = readdir([data_folder, proband, '\', session]);%gets names of c3d c3d files
    c3dnames_1  = sort(c3dnames_1);
    nr_c3d      = length(c3dnames_1);
    
    %rest_emg
    restc3d1  = strfind(c3dnames_1,'Ruhe');
    restc3d2  = find(~cellfun(@isempty,restc3d1'));
    restc3d   = c3dnames_1(restc3d2){1}
    file      = ezc3dRead(restc3d); %loads c3d into variable 'file'
##    restfiles.(proband).(session).(restc3d) = file;
    restfiles.(proband).(session) = restc3d;
    [an_fR, an_fF, an_lF] = analogsInfo(file);%gets basic information for the analogs (emg)
    for mus = 1:nr_mus
      muscle = muscles{mus};
      switch muscle%switches sensors - has to be similar to EMG_rest
        case 'L_vast_lat'
          if (strfind(proband, 'p_01')&&strfind(session, 'Session_2'))>0 || (strfind(proband, 'p_04')&&strfind(session, 'Session_1'))>0 || (strfind(proband, 'p_07'))>0 || (strfind(proband, 'p_09'))>0 || (strfind(proband, 'p_10')&&strfind(session, 'Session_1'))>0 || (strfind(proband, 'p_11')&&strfind(session, 'Session_1'))>0
            muscle = '9_unused';
          end%if
        case 'L_bic_fem'
          if (strfind(proband, 'p_02')&&strfind(session, 'Session_2'))>0 ||(strfind(proband, 'p_04')&&strfind(session, 'Session_2'))>0||(strfind(proband, 'p_05')&&strfind(session, 'Session_2'))>0||(strfind(proband, 'p_08')&&strfind(session, 'Session_2'))>0||(strfind(proband, 'p_09'))>0||(strfind(proband, 'p_10')&&strfind(session, 'Session_2'))>0
            muscle = '10_unused';
          end%if
      end%switch
      
      EMG_raw1          = collectANALOGSignals(file, muscle)';                                  %gets emg signal of muscle
      EMG_filtered1     = filterSignal_butter(EMG_raw1, 'high', an_fR, 'order', 4, 'cutoff', h_co);      %highpass filter
      EMG_demeaned      = EMG_filtered1 - mean(EMG_filtered1);                                        %demeaned
      EMG_rectified     = sqrt(EMG_demeaned.^2);                                                      %rectified
      EMG_filtered_low  = filterSignal_butter(EMG_rectified, 'low', an_fR, 'order', 4, 'cutoff', l_co);  % 4th order low-pass Butterworth filter 6 Hz
      EMG_filtered_low(EMG_filtered_low<0) = 0;                                                       %negative values resulting from low pass filter are set to zero
      
      EMG_rest_raw.(proband).(session)(mus,:) = EMG_raw1;
      EMG_rest.(proband).(session)(mus,:)     = EMG_filtered_low;
      
      EMG_rest_rms.(proband).(session)(mus,:) = rms(EMG_filtered_low(500:end-500));
      if EMG_rest_rms.(proband).(session)(mus,:) > 0.02
        display(['EMG rest too big for ', proband, ' ', session, ' ', muscles{mus}]); %if this happens control the signals  of rest emg 
      end%if
    end%for-->mus
    f = 0;
    %walk conditions
    for con = 1:nr_con
      condition = conditions{con}
      
      contrls1  = strfind(c3dnames_1,condition);
      contrls2  = find(~cellfun(@isempty,contrls1'));
      contrls   = c3dnames_1(contrls2);
      nr_ctr = length(contrls);
      for ctr = 1:nr_ctr % for each trial within one condition
        c3d = contrls{ctr};
        trial = c3d(1:end-4);
        clear file;
        file  = ezc3dRead(c3d); %loads c3d into variable 'file'
##        walkfiles.(proband).(session).(trial) = file; %saves all c3d files in variable 'walkfiles'
        f = f+1;
        walkfilenames.(proband).(session){f,1} = trial;
        
        [an_fR, an_fF, an_lF]     = analogsInfo(file);%gets basic information for the analogs (emg)
        [pts_fR, pts_fF, pts_lF]  = pointsInfo(file); %gets basic information for the points (f.example the markers, angles)
        anmulti                   = an_fR/pts_fR;     %ratio between framerates
        
        all_events = getEventTimes(file, 'Foot Strike', 'Right', 2);  %gets Foot Strikes from right leg in seconds
        all_events_pts = round(all_events*pts_fR)+1-pts_fF+1;         %in frames of points channels
        nr_stp = length(all_events)-1;                                %nr of steps within this trial
        
        for stp = 1:nr_stp
          stepframes_an.(proband).(session).(trial)(stp,:) = all_events_pts(stp:stp+1)*anmulti;
          stepframes_pt.(proband).(session).(trial)(stp,:) = all_events_pts(stp:stp+1);
          steptimes.(proband).(session).(trial)(stp,:)     = (stepframes_an.(proband).(session).(trial)(stp,2) -stepframes_an.(proband).(session).(trial)(stp,1))/an_fR;         
          
          for mus = 1:nr_mus
            muscle = muscles{mus};
            switch muscle%switches sensors - has to be similar to EMG_rest
              case 'L_vast_lat'
                if (strfind(proband, 'p_01')&&strfind(session, 'Session_2'))>0 || (strfind(proband, 'p_04')&&strfind(session, 'Session_1'))>0 || (strfind(proband, 'p_07'))>0 || (strfind(proband, 'p_09'))>0 || (strfind(proband, 'p_10')&&strfind(session, 'Session_1'))>0 || (strfind(proband, 'p_11')&&strfind(session, 'Session_1'))>0
                  muscle = '9_unused';
                end%if
              case 'L_bic_fem'
                if (strfind(proband, 'p_02')&&strfind(session, 'Session_2'))>0 ||(strfind(proband, 'p_04')&&strfind(session, 'Session_2'))>0||(strfind(proband, 'p_05')&&strfind(session, 'Session_2'))>0||(strfind(proband, 'p_08')&&strfind(session, 'Session_2'))>0||(strfind(proband, 'p_09'))>0||(strfind(proband, 'p_10')&&strfind(session, 'Session_2'))>0
                  muscle = '10_unused';
                end%if
            end%switch
            if (strfind(proband, 'p_07') && strfind(session, 'Session_1'))>0
              if strfind(condition, 'Met_indiv') > 0 || strfind(condition, 'Bod_indiv') > 0 || strfind(condition, 'Kombi_indiv') > 0
                muscle = sensors_p_07_session_1{mus};
              elseif strfind(trial, 'Einh_kombi_01') > 0 || strfind(trial, 'Einh_kombi_02') > 0
                muscle = sensors_p_07_session_1{mus};                  
              end%if
            end%if
            if (strfind(proband, 'p_03') && strfind(session, 'Session_2'))>0
              if strfind(condition, 'Met_indiv') > 0 || strfind(condition, 'Kombi_indiv') > 0 || strfind(trial, 'Bod_indiv_06') > 0 || strfind(trial, 'Bod_indiv_07') > 0 || strfind(trial, 'Bod_indiv_08') > 0 || strfind(trial, 'Bod_indiv_09') > 0 || strfind(trial, 'Bod_indiv_10') > 0 
                muscle = sensors_p_03_session_2{mus};
              end%if
            end%if
            
            EMG_raw1          = collectANALOGSignals(file, muscle)';                                          %gets emg signal of muscle
            EMG_filtered1     = filterSignal_butter(EMG_raw1, 'high', an_fR, 'order', 4, 'cutoff', h_co);     %highpass filter
            EMG_demeaned      = EMG_filtered1 - mean(EMG_filtered1);                                          %demeaned
            EMG_rectified     = sqrt(EMG_demeaned.^2);                                                        %rectified
##            EMG_filtered_low  = filterSignal_butter(EMG_rectified, 'low', an_fR, 'order', 4, 'cutoff', l_co);  % 4th order low-pass Butterworth filter l_co Hz --> choose either this line or the next (other in comments = strg+r; uncomment = strg+shift+r)
            EMG_filtered_low  = filterSignal_butter(EMG_rectified, 'low', an_fR, 'order', 4, 'cutoff', l_co/steptimes.(proband).(session).(trial)(stp,:));  % 4th order low-pass Butterworth filter l_co/duration Hz
            EMG_baseclear     = EMG_filtered_low - EMG_rest_rms.(proband).(session)(mus,:);
            EMG_baseclear(EMG_baseclear<0) = 0;                                                             %negative values resulting from low pass filter are set to zero
            timen             = timeNormalizeSig(EMG_baseclear(stepframes_an.(proband).(session).(trial)(stp,1):stepframes_an.(proband).(session).(trial)(stp,2)), 101);  %time normalize signal to 101 timepoints = 100%
            timen(timen<0)    = 0;
            
            
            EMG_raw.(proband).(session).(trial){stp,1}(mus,:)      = EMG_raw1(stepframes_an.(proband).(session).(trial)(stp,1):stepframes_an.(proband).(session).(trial)(stp,2));
            EMG_timenorm.(proband).(session).(trial){stp,1}(mus,:) = timen;
            
            if con == 1 && ctr == 1
              EMG_peak.(proband).(session)(mus,1) = max(timen);
              EMG_peaktrial.(proband).(session){mus,1} = trial;
            else
              if max(timen)> EMG_peak.(proband).(session)(mus,1);
                EMG_peak.(proband).(session)(mus,1) = max(timen);
                EMG_peaktrial.(proband).(session){mus,1} = trial;
              end%if
            end%if
          end%for-->mus
          %%%%joint angles
          jn3 = 0;
          for jnt = 1:nr_jnt
            ikraw = collectPointsSignals(file, jointsc3d{jnt});
            for jn2 = 1:length(jointnames{jnt});
              jn3 = jn3+1;
              ikraw2(1,:) = ikraw(jn2,1,:);
              ikfilt      = filterSignal_butter(ikraw2, 'low', pts_fR, 'order', 4, 'cutoff', 6); %filter with low pass filter (I recommend between 6 and 10)
              timen       = timeNormalizeSig(ikfilt(stepframes_pt.(proband).(session).(trial)(stp,1):stepframes_pt.(proband).(session).(trial)(stp,2)), 101);  %time normalize signal to 101 timepoints = 100%
              
              IK_raw.(proband).(session).(trial){stp,1}(jn3,:)      = ikraw2(stepframes_pt.(proband).(session).(trial)(stp,1):stepframes_pt.(proband).(session).(trial)(stp,2));
              IK_timenorm.(proband).(session).(trial){stp,1}(jn3,:) = timen;
              joints{jn3,1} = jointnames{jnt}{jn2};
              clear ikraw2;
            end%for-->jn2
          end%for-->jnt
        end%for-->stp    
      end%for --> ctr
    end%for-->con
    
    %amplnorm emg and reorder data
    trialnames = walkfilenames.(proband).(session);
    
    for con = 1:nr_con
      condition = conditions{con};
      stp2 = 0;
      contrls1  = strfind(trialnames,condition);
      contrls2  = find(~cellfun(@isempty,contrls1'));
      contrls   = trialnames(contrls2);
      nr_ctr = length(contrls);
      for ctr = 1:nr_ctr % for each trial within one condition
        trial = contrls{ctr};
        for stp = 1:length(EMG_timenorm.(proband).(session).(trial));        
          amplNorm  = 1./EMG_peak.(proband).(session);                            %multiplikator for amplnorm
          anorm     = EMG_timenorm.(proband).(session).(trial){stp,1}.*amplNorm;  %norm signal
          anorm(anorm<0) = 0;                                                     %set neg values (should not occour here) to 0
          stp2 = stp2+1;
          EMG_timenorm.(proband).(session).(trial){stp,1}     = anorm;
          EMG_contrl.(proband).(session).(condition){stp2,1}  = anorm;
          IK_contrl.(proband).(session).(condition){stp2,1}   = IK_timenorm.(proband).(session).(trial){stp,1};        
        end%for-->stp
      end%for-->ctr
      stepcount.(proband).(session).(condition) = length(EMG_contrl.(proband).(session).(condition));
      EMG_con.(proband).(session).(condition) = horzcat(EMG_contrl.(proband).(session).(condition){:});
      IK_con.(proband).(session).(condition)  = horzcat(IK_contrl.(proband).(session).(condition){:});
      
      Estorage{con} = horzcat(EMG_contrl.(proband).(session).(condition){:});
      Istorage{con} = horzcat(IK_contrl.(proband).(session).(condition){:});
    end%for-->con
    EMG_all.(proband).(session) = horzcat(Estorage{:});
    IK_all.(proband).(session)  = horzcat(Istorage{:});
  end%for --> ses
  clear c3dnames_1 nr_c3d restc3d1 restc3d2 restc3d file an_fR an_fF an_lF EMG_raw1 EMG_filtered1 EMG_demeaned EMG_rectified EMG_filtered_low;
  clear contrls1 contrls2 contrls nr_ctr c3d trial pts_fR pts_fF pts_lF anmulti all_events all_events_pts nr_stp EMG_baseclear timen;
  clear ikraw ikraw2 ikfilt amplNorm anorm Estorage Istorage;
end%for --> prb

cd (save_folder);

save -mat7-binary 'EMG_restpeak.mat' 'EMG_rest_raw' 'EMG_rest' 'EMG_rest_rms' 'EMG_peak' 'EMG_peaktrial';
##save -mat-binary 'c3dfiles.mat' 'restfiles' 'walkfiles';
save -mat-binary 'EMG_matrix.mat' 'EMG_raw' 'EMG_timenorm' 'EMG_contrl' 'EMG_con' 'EMG_all';
save -mat-binary 'IK_matrix.mat' 'IK_raw' 'IK_timenorm' 'IK_contrl' 'IK_con' 'IK_all';
save -mat7-binary 'events.mat' 'stepframes_an' 'stepframes_pt' 'steptimes' 'stepcount';
save -mat7-binary 'metadata.mat' 'probands' 'muscles' 'joints' 'conditions';
