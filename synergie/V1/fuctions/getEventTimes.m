function [event_times]= getEventTimes(dataSet, term, context, row) %dataSet...c3d data file; term...name of the event--> has to be written in ''; row... row, in which the event is located (2 for right foot strike) located in folder dataSet.parameters.EVENT.TIMES.DATA(row,:)
 event_labels= dataSet.parameters.EVENT.LABELS.DATA';  %gets the labels file
 event_contexts = dataSet.parameters.EVENT.CONTEXTS.DATA'; %gets the labels contexts (f.ex.: general, right, left)
 event_times_raw= dataSet.parameters.EVENT.TIMES.DATA(row,:); %gets the file with the times (correct row)
 
if size (event_labels,1) > 1 %if new
    event_labels = event_labels';
end

 for i=1:length(event_labels)
    switch (event_contexts{i})
      case context
        event_labels_context{1,i} = event_labels{1,i};
    end%switch
 end%for
 for i=1:length(event_labels_context)
    switch (event_labels_context{i})
      case term
      event_term_times_raw(1,i)= event_times_raw(1,i);
    end%switch
    %if the event label is the same as the term, it is copied to the event_term_times_raw matrix to the same position, if not, it will bee zero, at this position
  end%for
  event_times_with_zeros = sort(event_term_times_raw); %sorts the times in correct order, attention: zeros are at the beginning
  event_times= nonzeros(event_times_with_zeros)'; %clears the zeros, so only the times of the Foot Strike are in the matrix, with the first starting at (1,1), row style
end%function

%Author: Kaufmann Paul