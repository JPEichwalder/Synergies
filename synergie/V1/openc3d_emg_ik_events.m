clear         %clears workspace
clc           %clears command window
close all     %close all figures

%pkg load signal %signal package is needed for filter functions-->has to be downloaded (pkg install -forge signal)

% add functions path
functionsPath = addpath(genpath('C:\Users\jpeic\OneDrive\Documents\A_University\M_Rehatechnik\A_Masterarbeit\Code\Synergies\synergie\V1\fuctions\')); % replace with your folderpath with scripts and functions, can contain subfolders with scripts ('\' is important at the end of the path to read subfolders); have all your octave scripts and functions in this folder

% directory of c3d files and where to save the data --> replace with your path (always '\'at the ending)
data_folder = ['D:\VR_Synergy\PostProc\c3d\'];
save_folder = ['D:\VR_Synergy\PostProc\save\'];

%conditions that should be analyzed

% Conditions sind die Überbegriffe der Trails die in die gleiche Kategorie fallen (walking, balance, balance high) VR1,2,3 R1 1,2 
conditions  = {'R1_walking_'; 'R1_balance_'; 'VR1_walking_'; 'VR1_balance_'; 'VR1_balance_high_'}; % changing this names = change names in the script
nr_con      = length(conditions);

patternlist = {'^R1_walking_\d\d'; '^R1_balance_\d\d'; '^VR1_walking_\d\d'; '^VR1_balance_\d\d'; '^VR1_balance_high_\d\d'};


% muscles via sEMG that should be analyzed in the order they should be shown (I recommend the way its now) --> name of sensors in vicon --> send me list with differneces in namings and muscles so I can change the script for these
muscles = {'EMG03_m_gast_med'; 'EMG04_m_glut_max'; 'EMG05_m_tib_ant'; 'EMG06_m_pero_long';  'EMG07_m_soleus'; 'EMG09_m_qad_med'; 'EMG10_m_qad_lat'; ...
           'EMG11_m_qad_fem'; 'EMG12_m_gracilis'; 'EMG13_m_tens_fasc_lat';  'EMG14_m_glut_med'; 'EMG15_m_bic_fem'; 'EMG16_m_semiten'};
nr_mus  = length(muscles);

muscles_VAP3 = {'EMG01_m_tib_ant'; 'EMG02_m_pero_long'; 'EMG03_gast_med'; 'EMG04_m_glut_max';  'EMG07_m_soleus';  'EMG09_m_qad_med'; 'EMG10_m_qad_lat'; ...
                   'EMG11_m_qad_fem'; 'EMG12_m_gracilis'; 'EMG13_m_tens_fasc_lat'; 'EMG14_m_glut_med'; 'EMG15_m_bic_fem'; 'EMG16_m_semiten'};
    
muscles_VAP1_2 = {'EMG04_m_glut_max'; 'EMG05_m_tib_ant'; 'EMG06_m_pero_long'; 'EMG07_m_soleus'; 'EMG08_m_gast_med'; 'EMG09_m_qad_med'; 'EMG10_m_qad_lat'; ...
           'EMG11_m_qad_fem'; 'EMG12_m_gracilis'; 'EMG13_m_tens_fasc_lat'; 'EMG14_m_glut_med'; 'EMG15_m_bic_fem'; 'EMG16_m_semiten'};



                      
% joints that should be analyzed
jointsc3d   = {'LAnkleAngles_M';                              'LKneeAngles_M';                            'LHipAngles_M';                           'LPelvisAngle';...
              'RAnkleAngles_M';                               'RKneeAngles_M';                            'RHipAngles_M';                           'RPelvisAngle'};
jointnames  = {{'L_ankle_flex'; 'L_ankle_ad'; 'L_ankle_rot'}; {'L_knee_flex'; 'L_knee_ad'; 'L_knee_rot'}; {'L_hip_flex'; 'L_hip_ad'; 'L_hip_rot'};  {'L_pelvis_tilt'; 'L_pelvis_obli'; 'L_pelvis_rot'};...
               {'R_ankle_flex'; 'R_ankle_ad'; 'R_ankle_rot'}; {'R_knee_flex'; 'R_knee_ad'; 'R_knee_rot'}; {'R_hip_flex'; 'R_hip_ad'; 'R_hip_rot'};  {'R_pelvis_tilt'; 'R_pelvis_obli'; 'R_pelvis_rot'}};
nr_jnt = length(jointsc3d);
h_co = 20;  %cutoff frequency for high pass filter (I recommend 20 - 35)
l_co = 9;   %cutoff frequency for low pass filter (I recommend 7 - 14 or duration dependent)

% DA

dir_info = dir(data_folder);
probands = {dir_info.name};

file_names = {dir_info.name};

probands = cell(1, numel(file_names)); % Initialisiere probands

prb = 1;
for file = 1:numel(file_names)
    if file_names{file}(1) == 'v' || file_names{file}(1) == 'V'
        probands{prb} = file_names{file};
        disp(probands); % Hier hinzugefügt, um den Inhalt von probands zu überprüfen
        prb = prb + 1;
    end
end

probands = probands(~cellfun('isempty', probands)); % Entferne leere Einträge
nr_prb = length(probands);
clear probands_1 prb file;

for prb = 1:nr_prb
  proband = probands{prb};
    
    cd ([data_folder, proband, '\new\',]); % cd = change to specific folder
    c3dnames_1  = dir('*.c3d');%gets names of c3d c3d files
    c3dnames_1 = {c3dnames_1.name};
    c3dnames_1  = sort(c3dnames_1);

    nr_c3d      = length(c3dnames_1);
    
    f = 0;
    %walk conditions

    for con = 1:nr_con
      pattern = patternlist{con};
      condition = conditions{con};
      
        matchingFilenames = {};

            % Loop through each filename and check for a match
        for i = 1:length(c3dnames_1)
             if ~isempty(regexp(c3dnames_1{i}, pattern, 'once'))
                % If a match is found, add the filename to the matchingFilenames list
                 matchingFilenames{end+1} = c3dnames_1{i};
            end
        end

      contrls1 = matchingFilenames;

      contrls   = contrls1;

      nr_ctr = length(contrls);
      for ctr = 1:nr_ctr % for each trial within one condition
        c3d = contrls{ctr};
        trial = c3d(1:end-4);
        clear file;
        file  = ezc3dRead(c3d); %loads c3d into variable 'file'
        walkfiles.(proband).(trial) = file; %saves all c3d files in variable 'walkfiles'
        f = f+1;
        walkfilenames.(proband){f,1} = trial;
        
        [an_fR, an_fF, an_lF]     = analogsInfo(file);%gets basic information for the analogs (emg)
        [pts_fR, pts_fF, pts_lF]  = pointsInfo(file); %gets basic information for the points (f.example the markers, angles)
        anmulti                   = an_fR/pts_fR;     %ratio between framerates
        
        all_events = getEventTimes(file, 'Foot Strike', 'Right', 2);  %gets Foot Strikes from right leg in seconds
        all_events_pts = round(all_events*pts_fR)+1-pts_fF+1;         %in frames of points channels
        nr_stp = length(all_events)-1;                                %nr of steps within this trial
        
        for stp = 1:nr_stp
          stepframes_an.(proband).(trial)(stp,:) = all_events_pts(stp:stp+1)*anmulti;
          stepframes_pt.(proband).(trial)(stp,:) = all_events_pts(stp:stp+1);
          steptimes.(proband).(trial)(stp,:)     = (stepframes_an.(proband).(trial)(stp,2) -stepframes_an.(proband).(trial)(stp,1))/an_fR;         
          
        if strcmp(proband, 'VAP_01') || strcmp(proband, 'VAP_02')
            muscle_set = muscles_VAP1_2;
           elseif strcmp(proband, 'VAP_03')
            muscle_set = muscles_VAP3;
           else
            muscle_set = muscles;
        end


          
          for mus = 1:length(muscle_set)
            muscle = muscle_set{mus};
            
            EMG_raw1          = collectANALOGSignals(file, muscle)';                                          %gets emg signal of muscle
            EMG_filtered1     = filterSignal_butter(EMG_raw1, 'high', an_fR, 'order', 4, 'cutoff', h_co);     %highpass filter
            EMG_demeaned      = EMG_filtered1 - mean(EMG_filtered1);                                          %demeaned
            EMG_rectified     = sqrt(EMG_demeaned.^2);                                                        %rectified
            EMG_filtered_low  = filterSignal_butter(EMG_rectified, 'low', an_fR, 'order', 4, 'cutoff', l_co);  % 4th order low-pass Butterworth filter l_co Hz --> choose either this line or the next (other in comments = strg+r; uncomment = strg+shift+r)
            timen             = timeNormalizeSig(EMG_filtered_low(stepframes_an.(proband).(trial)(stp,1):stepframes_an.(proband).(trial)(stp,2)), 101);  %time normalize signal to 101 timepoints = 100%
            timen(timen<0)    = 0;
            
            
            EMG_raw.(proband).(trial){stp,1}(mus,:)      = EMG_raw1(stepframes_an.(proband).(trial)(stp,1):stepframes_an.(proband).(trial)(stp,2));
            EMG_timenorm.(proband).(trial){stp,1}(mus,:) = timen;
            
            if con == 1 && ctr == 1
              EMG_peak.(proband)(mus,1) = max(timen);
              EMG_peaktrial.(proband){mus,1} = trial;
            else
              if max(timen)> EMG_peak.(proband)(mus,1)
                EMG_peak.(proband)(mus,1) = max(timen);
                EMG_peaktrial.(proband){mus,1} = trial;
              end%if
            end%if
          end%for-->mus


          %%%%joint angles
          jn3 = 0;
          for jnt = 1:nr_jnt
            ikraw = collectPointsSignals(file, jointsc3d{jnt}); %gleich umschreiben wie analog signals
            for jn2 = 1:length(jointnames{jnt})
              jn3 = jn3+1;
              ikraw2(1,:) = ikraw(jn2,1,:); % Problem!!!!
              ikfilt      = filterSignal_butter(ikraw2, 'low', pts_fR, 'order', 4, 'cutoff', 8); %filter with low pass filter (I recommend between 6 and 10)
              timen       = timeNormalizeSig(ikfilt(stepframes_pt.(proband).(trial)(stp,1):stepframes_pt.(proband).(trial)(stp,2)), 101);  %time normalize signal to 101 timepoints = 100%
              
              IK_raw.(proband).(trial){stp,1}(jn3,:)      = ikraw2(stepframes_pt.(proband).(trial)(stp,1):stepframes_pt.(proband).(trial)(stp,2));
              IK_timenorm.(proband).(trial){stp,1}(jn3,:) = timen;
              joints{jn3,1} = jointnames{jnt}{jn2};
              clear ikraw2;
            end%for-->jn2
            %disp(jointsc3d{jnt}); 
          end%for-->jnt
          %disp([trial '... -> finished']); 

        end%for-->stp
        disp([trial '']); 
      end%for --> ctr
    end%for-->con
    
    %amplnorm emg and reorder data
    trialnames = walkfilenames.(proband);
    


    for con = 1:nr_con
      pattern = patternlist{con};
      condition = conditions{con};
      
        matchingFilenames = {};

            % Loop through each filename and check for a match
        for i = 1:length(c3dnames_1)
             if ~isempty(regexp(c3dnames_1{i}, pattern, 'once'))
                % If a match is found, add the filename to the matchingFilenames list
                 matchingFilenames{end+1} = c3dnames_1{i};
            end
        end

      contrls1 = matchingFilenames;

      contrls   = contrls1;

      nr_ctr = length(contrls);



      for ctr = 1:nr_ctr % for each trial within one condition
        trial = contrls{ctr};
        trial = trial(1:end-4);
        stp2 = 0;
        for stp = 1:length(EMG_timenorm.(proband).(trial))        
          amplNorm  = 1./EMG_peak.(proband);                            %multiplikator for amplnorm
          anorm     = EMG_timenorm.(proband).(trial){stp,1}.*amplNorm;  %norm signal
          anorm(anorm<0) = 0;                                                     %set neg values (should not occour here) to 0
          stp2 = stp2+1;
          EMG_timenorm.(proband).(trial){stp,1}     = anorm;
          EMG_contrl.(proband).(condition){stp2,1}  = anorm;
          IK_contrl.(proband).(condition){stp2,1}   = IK_timenorm.(proband).(trial){stp,1};        
        end%for-->stp
      end%for-->ctr
      disp([trial ' Angles finished']); 
      stepcount.(proband).(condition) = length(EMG_contrl.(proband).(condition));
      EMG_con.(proband).(condition) = horzcat(EMG_contrl.(proband).(condition){:});
      IK_con.(proband).(condition)  = horzcat(IK_contrl.(proband).(condition){:});
      
      Estorage{con} = horzcat(EMG_contrl.(proband).(condition){:});
      Istorage{con} = horzcat(IK_contrl.(proband).(condition){:});
    end%for-->con
    EMG_all.(proband) = horzcat(Estorage{:});
    IK_all.(proband)  = horzcat(Istorage{:});
    disp([proband ' -> finished']);

  clear c3dnames_1 nr_c3d restc3d1 restc3d2 restc3d file an_fR an_fF an_lF EMG_raw1 EMG_filtered1 EMG_demeaned EMG_rectified EMG_filtered_low;
  clear contrls1 contrls2 contrls nr_ctr c3d trial pts_fR pts_fF pts_lF anmulti all_events all_events_pts nr_stp EMG_baseclear timen;
  clear ikraw ikraw2 ikfilt amplNorm anorm Estorage Istorage;
end%for --> prb

cd (save_folder);

save(fullfile(save_folder, 'EMG_peak.mat'), 'EMG_raw', 'EMG_all', 'EMG_con', 'EMG_peak', 'EMG_peaktrial', 'EMG_contrl', 'EMG_timenorm');
save(fullfile(save_folder, 'c3dfiles.mat'), 'walkfiles');
save(fullfile(save_folder, 'EMG_matrix.mat'), 'EMG_raw', 'EMG_timenorm', 'EMG_contrl', 'EMG_con', 'EMG_all');
save(fullfile(save_folder, 'IK_matrix.mat'), 'IK_raw', 'IK_timenorm', 'IK_contrl', 'IK_con', 'IK_all');
save(fullfile(save_folder, 'events.mat'), 'stepframes_an', 'stepframes_pt', 'steptimes', 'stepcount');
save(fullfile(save_folder, 'metadata.mat'), 'probands', 'muscles', 'joints', 'conditions');


