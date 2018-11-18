function [p,stats ] =  powerttest(dist1,dist2,randdraw)
% performs 2sample t-test or Wilcoxon ranksum if data not normally
% distributed.
%
% dist1 and dist2 two are the vectors you want to compare by ttest. 
%
% randomdraw takes 1 or 0 as imputs. If 1, ttest/ranksum is performed on a
% random draw of each sample, the size of the randomly drawn sample is set
% to achieve a power of .80
% If randdraw=0, the sample will be drawn from the beginning and end of the
% distirbution. This is because i used this to compare beginning and end of
% saccade amplitude in saccadic adaption. 

if isempty(randdraw)
    randdraw=1;
end

  pwr = .80;
  mean1 = mean(dist1);
  mean2 = mean(dist2);
  
  std1 = std(dist1);
  std2 = std(dist2);
  
  % power to detected 10% difference in means
  difference = 1.1;
 
  % determine number of samples neede to determine 10% difference 
  % with .80 power
  nout = sampsizepwr('t',[mean1 std1 ],mean1*difference,pwr);
  
  % see if there are enough samples, if not show how how large the power is
  if nout > numel(dist1) || nout> numel(dist2) 
      
      % smallest of the two distributions determines 
      minsize=  min([length(dist1), length(dist2)]);

      pwr = sampsizepwr('t',[mean1 std1 ],mean1*difference,[],minsize);
      nout= minsize;
      disp('sample size too small for power=.8, power in this smaple:')
      disp(pwr)
  end
   
   % if not normal distributed Wilcoxon ranksum, if normal distributed ttest
   if kstest(datasample(dist1,nout,'Replace',false))
        
          if randdraw
           x =  datasample(dist1,nout,'Replace',false) ;
           y = datasample(dist2,nout,'Replace',false) ;
           [p, h] = ranksum( x , y);
            testtype='ranksum';
          else
            [p, h] = ranksum( (dist1(1:nout)), (dist2(end-nout+1:end)));
            testtype='ranksum';
          end
   else
       x =  datasample(dist1,nout,'Replace',false) ;
       y = datasample(dist2,nout,'Replace',false) ;
           [h, p] = ttest2( x , y);
           
      testtype='ttest';

   end
   
stats.h=h; % hypothesis (0/1)
stats.p=p; % p-value
stats.pwr =pwr; % power of test
stats.n = nout; % number of sample used in ttesst/wilcoxon
stats.randdraw = randdraw; % samples randomly drawn from population (0/1)
stats.testtype=testtype; % testtype used. Wilcoxon/ttest