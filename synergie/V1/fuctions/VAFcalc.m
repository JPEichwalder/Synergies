function VAF=VAFcalc(muscleMatrix_all_Norm, synergyMatrix)
  
  [x,y]=size(synergyMatrix);
  
  for fa=1:y
    w = synergyMatrix(fa).W;
    h = synergyMatrix(fa).H;
    if length(h)<length(muscleMatrix_all_Norm)
      nr_trials = length(muscleMatrix_all_Norm)/length(h);
      for trl = 1:nr_trials
        mM = muscleMatrix_all_Norm(:,trl*length(h)-length(h)+1:trl*length(h));
        VAF(fa,trl) = (1 - (norm(mM - w*h)^2 / norm(mM)^2))*100;
        clear mM
      endfor
    else
      VAF(fa,1) = (1 - (norm(muscleMatrix_all_Norm - w*h)^2 / norm(muscleMatrix_all_Norm)^2))*100;
    endif
  endfor
endfunction

%calculates total Variance-accounted-for for a synergy and the muscle Matrix

% @ Paul Kaufmann