clear
clc
close all

% Add functions path
functionsPath = addpath(genpath('C:\Users\jpeic\OneDrive\Documents\A_University\M_Rehatechnik\A_Masterarbeit\Code\Synergies\loadsol\'));

% Add ezc3d path
ezc3dPath = addpath('C:\Users\jpeic\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\ezc3d');

% Directory of the folder with the raw c3d and loadsol data for each participant
data_folder = 'C:\Users\jpeic\OneDrive\Documents\A_University\M_Rehatechnik\A_Masterarbeit\Pilot_Study\Pilot_Alex\Pilot_all\';
save_folder = 'C:\Users\jpeic\OneDrive\Documents\A_University\M_Rehatechnik\A_Masterarbeit\Pilot_Study\Pilot_Alex\Pilot_test\';

cd(save_folder);
soles = load('soles.mat');
soles = soles.soles;

probands = fieldnames(soles);
nr_prb = length(probands);

for prb = 1:nr_prb
    clear trialnames
    proband = probands{prb};
    soletrials = fieldnames(soles.(proband)); % Trials of proband in soles struct
    nr_strl = length(soletrials);

    for strl = 1:nr_strl
        [solesevents.(proband).(soletrials{strl}).con, solesevents.(proband).(soletrials{strl}).off, nozeros] = det_contacts(soles.(proband).(soletrials{strl}), 'R_total', 'gen_th', 30, 'min_cycle', 200, 'max_cycle', 10000, 'transmission_error', 2, 'cycle_th', 100);
    end

    cd([data_folder, proband, '/c3d']); % Folder with raw c3d files

    c3dnames_1 = dir([data_folder, proband, '/c3d']); % Get names of c3d files
    c3dnames = c3dnames_1(3:end); % Skip '.' and '..' entries
    nr_c3d = length(c3dnames);
    a = 1;

    for c3 = 1:nr_c3d
        searchString = {'static', 'MVC', 'EMG'};
        found = false;

        % Get the current file name from the structure
        currentFileName = c3dnames(c3).name;

        % Loop through each search string and check if it exists in the file name
        for i = 1:numel(searchString)
            if ~isempty(strfind(currentFileName, searchString{i}))
                found = true;
                break; % Exit the loop early if a match is found
            end
        end

        % Check if none of the search strings were found
        if ~found
            trialnames{a, 1} = currentFileName(1:end-4);
            a = a + 1;
        end
    end

    nr_trl = length(trialnames);

    for trl = 1:nr_trl
        trial = trialnames{trl};

        con = solesevents.(proband).(trial).con;
        off = solesevents.(proband).(trial).off;

        % Open the C3D file for reading
        file = btkReadAcquisition([trial, '.c3d']);

        % Get the original events
        original_events = btkGetEvents(file);

        % Calculate condelay
        %condelay = con(1) - original_events.('Right_Foot_Strike').Right(1); %name the Events like in the c3d file event 
        condelay = con(1) - original_events.Right_Foot_Strike.Right(1);
        % Calculate contimes and offtimes
        contimes = con - condelay;
        offtimes = off - condelay;

        % Update events in the C3D file
        file = btkSetEvent(file, 'Foot Strike', 'Right', contimes);
        file = btkSetEvent(file, 'Foot Off', 'Right', offtimes);

        % Save the modified C3D file
        btkWrite([trial, '_modified.c3d'], file);

        % Move the modified file to the desired location
        movefile([trial, '_modified.c3d'], [save_folder, proband, '/new']);

        % Change the working directory back to the original C3D folder
        cd([data_folder, proband, '/c3d/new']); % Folder with raw c3d files
    end

    clear c3dnames_1 proband
end
