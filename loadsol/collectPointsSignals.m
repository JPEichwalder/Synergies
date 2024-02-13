function [points_raw]= collectPointsSignals(dataSet, pointslabel) %dataSet...c3dfile; pointslabel...label of the point

ind = index (dataSet.parameters.POINT.LABELS.DATA, pointslabel);
if sum(ind) == 0
  error('label is not existing in the file')
else
  [i,row] = max(ind); %finds index of label
  if length(row)>1
    error('label is not unique in the file')
  end%if
end%if

points_raw = dataSet.data.points(:,row,:); %gets data from the points folder with index (:,row,:); 3deminsional so ist smth like (3,row,datapoints)

end%function

%gets the points data from the c3d file for the input label
%Author: Paul Kaufmann