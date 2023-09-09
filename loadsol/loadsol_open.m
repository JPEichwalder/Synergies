function [Lstruct, Lcell, raw_file] = loadsol_open(filename, varargin) % loads ASCII file from loadsol data.

%Lstruct ... structure where data is stored; Lcell ... data stored as cell; raw_file...imported txt file with text in rows and just on collumn
%filename ... name of the .txt file to open; varargin ... if struct already exists, datas are stored in this --> must have the same name as Lstruct

Lstruct = struct();
proband{1,1} = 'proband';

while (numel(varargin)>0) %changes values from default to specified input values: *'inputname', value*
    switch class(varargin{end});
      case 'struct'
        Lstruct = varargin{end};
        varargin(end) = [];
      case 'char'
        proband{1,1} = varargin{end};
        varargin(end) = [];
      otherwise
        error ("check input arguments");
    end%switch
  end%while

filestr = importdata(filename); %loads file as stings in cells per rows
raw_file = filestr;

nr_rows = length(filestr)-4; %number of rows = nr of datas + headings - empty appendix in ascii file

for he = 1:4  % 4 rows of headings in cells with strings
    head = textscan(filestr{he}, '%s%s%s%s%s%s%s%s%s%s%s%s%s%s', 'delimiter', '\t', 'EndOfLine', ''); %takes strings as strings seperated by tabs
    for col = 1:14
      if length(head{col}) == 1
        storage = head{col};
        Lcell{he,col} = storage{1,1};
      else
        Lcell{he,col} = head{col};
      end
    end
    clear head;
end

switch proband{1,1}
  case 'proband'
    a = strfind(Lcell{1,1}, '_23');
    if a<1
      a = 18;
    end%if
    proband{1,1} = Lcell{1,1}(7:a-1); %proband name
end%switch


if length(Lcell{2,1}) == 8 %no comment in the file
  trial{1,1} = 'trial'; %trialname = trial
else
  trial{1,1} = Lcell{2,1}(9:end); %trialname = comment
end

namesprb = fieldnames(Lstruct); %names of probands in struct
if length(namesprb)>=1 %if at least 1 proband exists
  if sum(strcmp(proband, namesprb)==1)>=1 %if one proband is the exact same as proband{1,1}
    namestrl = fieldnames(Lstruct.(proband{1,1})); %already existing trials of proband{1,1}
    same = strfind(namestrl, trial);
    nr_same = sum(cell2mat(same)==1); %nr of trials that are the same or longer than trial{1,1}
    if nr_same == 1 %if only 1 trial has the same name
      trial{2,1} = [trial{1,1}, '_c002']; %trial{2,1} = name of trial if datas in existing and this trial is not the same
      trial{3,1} = [trial{1,1}, '_c001']; %if data is not the same, then existing data has to be renamed
    elseif nr_same > 1 && nr_same < 9
      trial{2,1} = [trial{1,1}, '_c00', num2str(nr_same+1)];
      trial{3,1} = [trial{1,1}, '_c001'];
    elseif nr_same > 1 && nr_same > 8 && nr_same < 99
      trial{2,1} = [trial{1,1}, '_c0', num2str(nr_same+1)];
      trial{3,1} = [trial{1,1}, '_c001'];
    elseif nr_same > 1 && nr_same > 98
      trial{2,1} = [trial{1,1}, '_c', num2str(nr_same+1)];
      trial{3,1} = [trial{1,1}, '_c001'];
    elseif nr_same == 0
      clear nr_same;
    end%if
  end%if
end%if
parameters = {'L_time', 'L_heel', 'L_medial', 'L_lateral', 'L_total', 'R_time', 'R_heel', 'R_medial', 'R_lateral', 'R_total', 'Acc_time', 'AccX', 'AccY', 'AccZ'}; %headers of struct 


for dat = 5:nr_rows % for data rows
  datrow = filestr{dat}; %selects specific row
  datrow(datrow == ',') = '.'; %changes all commas in the row with dots
  data = textscan(datrow, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'delimiter', '\t', 'EndOfLine', '', 'EmptyValue', NaN); %takes strings as numbers
  for col = 1:14
    if length(data{col}) == 0 %if no number exists (cell 0x1) --> NaN
      data{col} = NaN;
    end
    Lcell(dat,col) = data(col); %creates cell with numbers
    datas(dat-4,col) = data{col}; %creates matrix with numbers (doubles)
  end
  clear datarow data 
end

if exist('nr_same') == 1 %if same trialname exists
  for col = 1:14 %creates Lstructfield with trialname trial{2,1}
  datastore = datas(:,col); %gets data per parameter
  nans = find(isnan(datastore));  % finds indexes of NaN's
  if length(nans) == 0 % no Nan values
    Lstruct.(proband{1,1}).(trial{2,1}).(parameters{col}) = datastore;
  else %NaN values exist
    firstnan = nans(1);
    Lstruct.(proband{1,1}).(trial{2,1}).(parameters{col}) = datastore(1:firstnan-1); % just not NaN's in the struct
  end%if
  clear datastore nans firstnan;
  end%for
  
  if nr_same == 1 %if only 1 trial same as trial[1,1} already exists
    equ = isequal(Lstruct.(proband{1,1}).(trial{2,1}), Lstruct.(proband{1,1}).(trial{1,1})); %checks if existing is the same as current trial
    if equ == 1 %is the same
      Lstruct.(proband{1,1}) = rmfield(Lstruct.(proband{1,1}), trial{2,1}); %remove current data from struct, cause it is already in
      disp([proband{1,1}, ' trialname: ', trial{1,1}, ' already exists with same data']);
    else %is not the same
      Lstruct.(proband{1,1}).(trial{3,1}) = Lstruct.(proband{1,1}).(trial{1,1}); %rename existing to _c001
      Lstruct.(proband{1,1}) = rmfield(Lstruct.(proband{1,1}), trial{1,1}); %remove named trial{1,1}
    end%if
  elseif nr_same > 1 %more than one already exist
    for trl = 1:length(namestrl)
    equ = isequal(Lstruct.(proband{1,1}).(trial{2,1}), Lstruct.(proband{1,1}).(namestrl{trl,1}));
    if equ == 1 %if it is eaqual remove current, else, let all in
      Lstruct.(proband{1,1}) = rmfield(Lstruct.(proband{1,1}), trial{2,1});
      disp([proband{1,1}, ' trialname: ', trial{1,1}, ' already exists with same data as trial: ' namestrl{trl,1}]);
    end%for
    end%if  
  end%if
  
  

else %no field with same trial name exists already
for col = 1:14
  datastore = datas(:,col); %gets data per parameter
  nans = find(isnan(datastore));  % finds indexes of NaN's
  if length(nans) == 0 % no Nan values
    Lstruct.(proband{1,1}).(trial{1,1}).(parameters{col}) = datastore;
  else %NaN values exist
    firstnan = nans(1);
    Lstruct.(proband{1,1}).(trial{1,1}).(parameters{col}) = datastore(1:firstnan-1); % just not NaN's in the struct
  end
 
  clear datastore nans firstnan;
end
end%if
Lstruct.(proband{1,1}) = orderfields(Lstruct.(proband{1,1}));%orders struct
end%function