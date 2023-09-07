function [points_frameRate, points_firstFrame, points_lastFrame] = pointsInfo(dataSet)
  points_frameRate=dataSet.header.points.frameRate;
  points_firstFrame=dataSet.header.points.firstFrame;
  points_lastFrame=dataSet.header.points.lastFrame;
end%function

%Author: Kaufmann Paul
% gets frame Rate, first and last Frame of the points, from a c3d file obend with ezc3d read