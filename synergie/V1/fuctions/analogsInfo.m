function [analogs_frameRate, analogs_firstFrame, analogs_lastFrame] = analogsInfo(dataSet)
  analogs_frameRate=dataSet.header.analogs.frameRate;
  analogs_firstFrame=dataSet.header.analogs.firstFrame;
  analogs_lastFrame=dataSet.header.analogs.lastFrame;
end%function
  
%Author: Kaufmann Paul
% gets frame Rate, first and last Frame of the analogs, from a c3d file obend with ezc3d read