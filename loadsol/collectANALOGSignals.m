function [analog]= collectANALOGSignals(dataSet, analogLabel) %dataSet...c3dfile; analogLabel...muscle

ind = index (dataSet.parameters.ANALOG.LABELS.DATA, [analogLabel]);
if sum(ind) == 0
  error('label is not existing in the file')
else
  [i,channel] = max(ind); %finds index of label
  if length(channel)>1
    error('label is not unique in the file')
  end%if
end%if

analog = dataSet.data.analogs(:,channel); %gets data from the analogs folder with index (channel)

end%function

%gets the EMG data from the c3d file for the input muscle
%Author: Paul Kaufmann