function [hF] = createConfLimSdf(alignedTrials,timeWindow,varargin)
%funtion that will create a plot with shaded confidence limits around the
%mean sdf for pro and anti overlayed on each other

p = inputParser;
p.addParamValue('hA',[],@(x) ishandle(x));
p.parse(varargin{:});
hA = p.Results.hA;
if isempty(hA)
    hF = figure('color','w');
    hA = axes('parent',hF);
else
    hF = get(hA,'parent');
end
%dirLabels = {'Pro','Anti'};
dirCc = [0 0.6 0;0.8 0 0]; %green fro pro red for anti
patchCc = [0 1 0;1 0 0];
%stick all sdfs on top of each other
proSdfs = vertcat(alignedTrials.Pro.sdf);
antiSdfs = vertcat(alignedTrials.Anti.sdf);
if ~isempty(alignedTrials.Pro)
sdfT = alignedTrials.Pro(1).sdfT; %careful will only work if all the same whih they should be but should have an error check
else
    sdfT = alignedTrials.Anti(1).sdfT; 
end
numProTrials = size(proSdfs,1);
numAntiTrials = size(antiSdfs,1);


if numAntiTrials == 1;
    antiMu = antiSdfs;
    antiSem = zeros(size(antiMu));
    
elseif numAntiTrials == 0
    antiMu = zeros(size(sdfT));
    antiSem = zeros(size(sdfT));
    
    
else
    
    antiMu = mean(antiSdfs);
    %exactly the same for anti
    %now get upper and lower confidence limits
    antiSem = std(antiSdfs)./sqrt(numAntiTrials);
    
end
antiUpp = antiMu+antiSem;
antiLow = antiMu-antiSem;

if numProTrials == 1
    proMu = proSdfs;
    proSem = zeros(size(proSdfs));
elseif numProTrials == 0
    proMu = zeros(size(sdfT));
    proSem = zeros(size(sdfT));
    
else
    
    proMu = mean(proSdfs);
    %now get upper and lower confidence limits
    proSem = std(proSdfs)./sqrt(numProTrials);
end

proUpp = proMu+proSem;
proLow = proMu-proSem;
%draw the mean as a line
line(sdfT,proMu,'parent',hA,'color',dirCc(1,:))
%and the confidence limits as a patch around it
patch([sdfT fliplr(sdfT)],[proUpp fliplr(proLow)],patchCc(1,:),'edgecolor','none','facealpha',0.5)


%draw the mean as a line
line(sdfT,antiMu,'parent',hA,'color',dirCc(2,:))
%and the confidence limits as a patch around it
patch([sdfT fliplr(sdfT)],[antiUpp fliplr(antiLow)],patchCc(2,:),'edgecolor','none','facealpha',0.5)

set(get(hA,'ylabel'),'string','Firing Rate (Hz)')
set(get(hA,'xlabel'),'string','Time (s)')
set(hA,'xlim',[-timeWindow(1) timeWindow(2)],'ycolor','k','xcolor','k')

topYlim = ceil(max([antiUpp proUpp]));
set(hA,'ylim',[0 topYlim])