function [signal_timenorm] = timeNormalizeSig(signal, dp) %signal= signal das zeitnormalisiert werden soll; dp...deserved number of points(auf soviele datenpunkte normalisieren)
  onp = length(signal); %original number of points; soviele sinds gerade
      r=onp/dp; % ratio ... verh�ltnis zwischen dp und onp
      for li= 1:dp
        s= 1 + r * (li - 1);
        
        if li == dp
          j= fix(s - 0.5); % rundet zum n�heste kleinsten integer (runden nach null)
        else
          j= round(s - 0.5); % rundet zum n�hesten integer
        end%if
        
        f = s - j;
        
        signal_timenorm(li) = (1 - f) * signal(j) + f * signal(j + 1);
       end%for
end%function

%Author: Kaufmann Paul
