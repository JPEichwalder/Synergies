function [EMG_raw]= collectEMGSignals(dataSet, electrode) %dataSet...c3dfile; electrode...muscle

[i,channel] = max(index (dataSet.parameters.ANALOG.LABELS.DATA, ['Voltage.',electrode])); %finds index of muscle

EMG_raw = dataSet.data.analogs(:,channel); %gets data from the analogs folder with index (channel)

end%function

%gets the EMG data from the c3d file for the input muscle
%Author: Paul Kaufmann