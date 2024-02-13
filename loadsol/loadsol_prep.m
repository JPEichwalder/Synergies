clear;
functionspath = addpath(genpath('C:\Users\jpeic\OneDrive\Documents\A_University\M_Rehatechnik\A_Masterarbeit\Code\loadsol')); %functionspath
cd ('D:\VR_Synergy\RecordingData\VAP_13\loadsole'); %datapath
a = dir(pwd);

j = {a(3).name};
disp(j{1});

[soles, cell2, raw2] = loadsol_open('VAP_13ASCII_23-10-23 18-06-21-899.txt','VAP_13');

for i = 4:numel(a)
    zelleInErsterSpalte = a(i).name;
    [soles, cell2, raw2] = loadsol_open(zelleInErsterSpalte, 'VAP_13', soles);
    fprintf('Zelle in der ersten Spalte: %s\n', zelleInErsterSpalte);
end

safepath = ('D:\VR_Synergy\RecordingData\VAP_13\loadsole\soles');
save(safepath, 'soles');

fprintf('Michi.ENDE');