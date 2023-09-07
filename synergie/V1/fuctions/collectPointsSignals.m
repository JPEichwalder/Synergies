function [points_raw]= collectPointsSignals(dataSet, pointslabel) %dataSet...c3dfile; pointslabel...label of the point

[i,row] = max(index (dataSet.parameters.POINT.LABELS.DATA, pointslabel)); %finds index of label

points_raw = dataSet.data.points(:,row,:); %gets data from the points folder with index (:,row,:); 3deminsional so ist smth like (3,row,datapoints)

end%function

%gets the points data from the c3d file for the input label
%Author: Paul Kaufmann