clear
clc
close all

functionsPath = addpath(genpath('C:\Users\jpeic\OneDrive\Documents\A_University\M_Rehatechnik\A_Masterarbeit\Code\Synergies\synergie\V1\fuctions\')); % replace with your folderpath with scripts and functions, can contain subfolders with scripts ('\' is important at the end of the path to read subfolders); have all your octave scripts and functions in this folder

save_folder = ['C:\Users\jpeic\OneDrive\Documents\A_University\M_Rehatechnik\A_Masterarbeit\Data\save\'];


cd (save_folder);
load('EMG_peak.mat');
load('EMG_matrix.mat');
load('events.mat');
load('metadata.mat');
load('muscles_muscleNumberation.mat');

probands    = probands;
nr_prb      = length(probands);
conditions  = conditions;
nr_con      = length(conditions);

muscles_set     = muscles;




nr_mus      = length(muscles_set);

EMG_timenorm  = EMG_timenorm;
time          = steptimes;

close all

for prb = 1:nr_prb
    proband = probands{prb};
    trails  = fieldnames(EMG_timenorm.(proband));
    nr_trl = length(trails);

    if strcmp(proband, 'VAP_01') || strcmp(proband, 'VAP_02')
        muscles_set = muscles_VAP1_2;
    elseif strcmp(proband, 'VAP_03')
        muscles_set = muscles_VAP3;
    else muscles_set = muscles;
    end
    
    for trl = 1:nr_trl
        trail = trails{trl};
        
        for stp = 1:length(EMG_timenorm.(proband).(trail))
            figure(1);
            sgtitle([proband, ' - ', trail, ' Step ', num2str(stp), ' Time ', num2str(time.(proband).(trail)(stp))]);
            set(gcf, 'Position', [1.6, 1.6, 638.6, 719.3]);

            for mus = 1:nr_mus
                subplot(2, nr_mus/2, mus, 'align');
                plot(0:1:100, EMG_timenorm.(proband).(trail){stp,1}(mus,:));
                title(muscles_set{mus});

                xlabel('Gait');
                if mus == 1
                    ylabel('Normalized EMG');
                end
                grid on;
            end

            figure(2);
            sgtitle([proband, ' - ', trail, ' Step ', num2str(stp), ' Time ', num2str(time.(proband).(trail)(stp))]);
            set(gcf, 'Position', [641.6,1.6,638.6,719.3]);

            for mus = 1:nr_mus
                subplot(2, nr_mus/2, mus, 'align');
                plot(EMG_raw.(proband).(trail){stp,1}(mus,:));
                title(muscles_set{mus});

                xlabel('Gait');
                if mus == 1 
                    ylabel('Raw EMG');
                end
                grid on;
            end
        
        end
    end
end


  
 %% nr_ses   = length(sessions);
    
%   for ses = 1:nr_ses
%     session = sessions{ses};
%     trials  = fieldnames(EMG_timenorm.(proband).(session));
%     nr_trl  = length(trials);
%     for trl = 1:nr_trl
%       trial = trials{trl};
%       for stp = 1:length(EMG_timenorm.(proband).(session).(trial));
%         for mus = 1:nr_mus
%           figure(1);
%           s = subplot(2,nr_mus/2,mus, 'align');
%           plot([0:1:100], EMG_timenorm.(proband).(session).(trial){stp,1}(mus,:));
%           title(muscles{mus});
%           
%           figure(2);
%           s = subplot(2,nr_mus/2,mus, 'align');
%           plot(EMG_raw.(proband).(session).(trial){stp,1}(mus,:));
%           title(muscles{mus});
%         end%for-->mus
%         figure(1); suptitle([proband, session, '     ', trial, ' step', num2str(stp), ' time ', num2str(time.(proband).(session).(trial)(stp))]);
%         figure(2); suptitle([proband, session, '     ', trial, ' step', num2str(stp), ' time ', num2str(time.(proband).(session).(trial)(stp))]);
%         breakpointhere=3;
%       end%for-->stp
%     end%for --> trl
%   end%for --> ses
% end%for --> prb
% 
