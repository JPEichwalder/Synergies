function [analog]= collectANALOGSignals(dataSet, analogLabel) %dataSet...c3dfile; analogLabel...muscle

[i,channel] = max(index (dataSet.parameters.ANALOG.LABELS.DATA, [analogLabel])); %finds index of muscle

analog = dataSet.data.analogs(:,channel); %gets data from the analogs folder with index (channel)

end%function

%gets the EMG data from the c3d file for the input muscle
%Author: Paul Kaufmann