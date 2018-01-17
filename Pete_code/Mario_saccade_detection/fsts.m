% fsts finds triggered ensembles of the signal.
%
% Written by Sungho Hong, CNS unit, OIST, 2013
function r = fsts(stim, triggers, Tpre, Tpost)
  too_early = (triggers<=Tpre);
  too_late = (triggers>=(length(stim)-Tpost));
  
  if sum(too_early)>0
    fprintf(1,'%d samples before Tpre -- ignored.', sum(too_early));
  end
  
  if sum(too_late)>0
    fprintf(1,'%d samples too late -- ignored.', sum(too_late));
  end
  
  t = triggers(~too_late & ~too_late);
  r = stim(uint32(bsxfun(@plus,-Tpre:Tpost,t)));

end
