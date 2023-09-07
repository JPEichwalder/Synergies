history -c
clear
clc
close all

% add functions path
functionsPath = addpath(genpath('D:\octave\')); %folder with skript and functions

% directory of the folder with the raw c3d and loasol data for each participant
data_folder = ['D:/master_data/raw/'];

save_folder = ['D:/master_data/events/'];
cd (save_folder);
load('soles_1.mat');
%gets probandnames
probands_1  = readdir([data_folder]);
prb = 1;
for fil = 1:length(probands_1)
  switch probands_1{fil}(1)
    case 'P'
      probands{prb,1} = probands_1{fil};
      prb = prb +1;
  end%switch
end%for
nr_prb      = length(probands);
clear probands_1;

for prb = 10%11:nr_prb
  proband = probands{prb};
  cd ([data_folder, proband, '/loadsol']);
  
  asciinames_1  = readdir([data_folder, proband, '/loadsol']);%gets names of ascii loadsol files
  asciinames    = asciinames_1((3:length(asciinames_1)),1);
  nr_asc        = length(asciinames);
  
  for asc = 1:nr_asc
    ascii = asciinames{asc};
    
    if asc == 1 && prb == 1
    [soles, celll, raw] = loadsol_open(ascii,proband);
    else
    [soles, celll, raw] = loadsol_open(ascii, soles, proband);
    end%if
  
    clear ascii celll raw;
  end%for
  clear asciinames_1;
end%for

cd (save_folder);
save -mat-binary 'soles_correct.mat' 'soles';
