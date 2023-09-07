history -c  %clears command history
clear       %clears workspace
clc         %clears command window
close all   %close all figures

% set breakpoints at lines with variables 'breakpointhere' --> F5 to every breakpoint

% add functions path
functionsPath = addpath(genpath('D:\octave\')); % replace with your folderpath with scripts and functions, can contain subfolders with scripts ('\' is important at the end of the path to read subfolders); have all your octave scripts and functions in this folder

% directory of c3d files and where to save the data --> replace with your path (always '\'at the ending)
data_folder = ['D:\synergies_constraints_july\c3d\'];
save_folder = ['D:\synergies_constraints_july\results'];

cd (save_folder);
load('EMG_restpeak.mat');
load('EMG_matrix.mat');
load('events.mat');
load('metadata.mat');

probands    = probands;
nr_prb      = length(probands);
conditions  = conditions;
nr_con      = length(conditions);
muscles     = muscles;
nr_mus      = length(muscles);

EMG_rest      = EMG_rest;
EMG_timenorm  = EMG_timenorm;
time          = steptimes;

%%%%%%%%%%% EMG_rest --> check if 500 to end-500 EMG has artefacts/high; if yes write down what parts we have to cut away
for prb = 1:nr_prb
  proband = probands{prb};
  
  sessions = fieldnames(EMG_rest.(proband));
  nr_ses   = length(sessions);
    
  for ses = 1:nr_ses
    session = sessions{ses};
    for mus = 1:nr_mus
      figure(1);
      s = subplot(2,nr_mus/2,mus, 'align');
      plot(EMG_rest_raw.(proband).(session)(mus,:));      
      hold all;
      plot(EMG_rest.(proband).(session)(mus,:), 'Linewidth', 4);
      title(muscles{mus});
    end%for-->mus
    figure(1); suptitle([proband, session]);
    breakpointhere=3;
  end%for --> ses
end%for --> prb

close all
%%%%%%%%%%% walkfiles --> write down which proband, session, trial, step (or all steps), muscle has artefacts and which might has artefacts --> so that we can later decide what we have to exclude/could exclude
for prb = 1:nr_prb
  proband = probands{prb};
  
  sessions = fieldnames(EMG_timenorm.(proband));
  nr_ses   = length(sessions);
    
  for ses = 1:nr_ses
    session = sessions{ses};
    trials  = fieldnames(EMG_timenorm.(proband).(session));
    nr_trl  = length(trials);
    for trl = 1:nr_trl
      trial = trials{trl};
      for stp = 1:length(EMG_timenorm.(proband).(session).(trial));
        for mus = 1:nr_mus
          figure(1);
          s = subplot(2,nr_mus/2,mus, 'align');
          plot([0:1:100], EMG_timenorm.(proband).(session).(trial){stp,1}(mus,:));
          title(muscles{mus});
          
          figure(2);
          s = subplot(2,nr_mus/2,mus, 'align');
          plot(EMG_raw.(proband).(session).(trial){stp,1}(mus,:));
          title(muscles{mus});
        end%for-->mus
        figure(1); suptitle([proband, session, '     ', trial, ' step', num2str(stp), ' time ', num2str(time.(proband).(session).(trial)(stp))]);
        figure(2); suptitle([proband, session, '     ', trial, ' step', num2str(stp), ' time ', num2str(time.(proband).(session).(trial)(stp))]);
        breakpointhere=3;
      end%for-->stp
    end%for --> trl
  end%for --> ses
end%for --> prb

