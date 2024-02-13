clear
history -c
clc

functionspath = addpath(genpath('C:\Users\jpeic\OneDrive\Documents\A_University\M_Rehatechnik\Masterarbeit\Code\')); %functionspath

cd ('D:\masterarbeit'); %datapath

    
[soles, cell2, raw2] = loadsol_open('line_01.txt');
[soles, cell2, raw2] = loadsol_open('line_01.txt', soles);
[soles, cell2, raw2] = loadsol_open('line_02.txt', soles);
[soles, cell2, raw2] = loadsol_open('line_03.txt', soles);
[soles, cell2, raw2] = loadsol_open('line_03_notsame.txt', soles);
[soles, cell2, raw2] = loadsol_open('noname1.txt', soles);
[soles, cell2, raw2] = loadsol_open('noname1.txt', soles);
[soles, cell2, raw2] = loadsol_open('noname2.txt', soles);
[soles, cell2, raw2] = loadsol_open('noname3.txt', soles);
[soles, cell2, raw2] = loadsol_open('noname3.txt', soles);
[soles, cell2, raw2] = loadsol_open('noname4.txt', soles);

%struct soles schould contain of 8 fields: line_01; line_02; line_03_c001; line_03_c002; trial_c001; trial_c002; trial_c003; trial_c004;

[con, off] = det_contacts(soles.Elias1.line_01, 'R_total', 'gen_th', 500, 'max_cycle', 7000); 
line_01.con = con; line_01.off = off;
[con, off] = det_contacts(soles.Elias1.line_02, 'R_total', 'con_th', 30, 'off_th', 26, 'max_cycle', 100000); 
line_02.con = con; line_02.off = off;
[con, off] = det_contacts(soles.Elias1.line_03_c001, 'L_total', 'gen_th', 20, 'min_cycle', 100); 
line_03_c001.con = con; line_03_c001.off = off;
%try some manipulation of input parameters

