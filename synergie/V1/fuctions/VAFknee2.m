function [NoS, angles] = VAFknee2(vaf, minVAF) %NoS = Nr of Synergies, at which the knee point is; angles = all angles between synergies; vaf = VAF matrix
  
  if exist('minVAF')
    miV = numel(vaf(vaf<minVAF));
  end%if
  for syn = 2:length(vaf)-1 %calc angles from 2 to lenght(vaf) - 1 synergies
    
    P1 = [syn-1,vaf(syn-1)];
    P2 = [syn,vaf(syn)];
    P3 = [syn+1,vaf(syn+1)];
    angl = rad2deg(angle3Points(P1,P2,P3));

    angles(syn,1)= angl;% all angles
  end%for
mins = sort(nonzeros(angles));
mi = 1;
mindex = find(angles == mins(mi));
if exist('minVAF')
  while mindex<=miV
    mi = mi+1;
    mindex = find(angles == mins(mi));
  end%while
end%if


NoS = mindex;
end%function

##  [minval, idx] = min(angles(2:length(angles)));%index of minimum value without idx 1 couse this is 0
##  NoS = idx+1; %idx+1 to get the synergie with the minimum value
  
%@Paul Kaufmann