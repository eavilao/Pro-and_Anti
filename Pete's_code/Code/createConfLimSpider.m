function [hAx] = createConfLimSpider(datVal,datErr,anglesToTest,hAx,plotLims)
%function that should plot 1 or more data sets in a polar plot along with a
%patch representing confidence limits. Is designed to work with Eric's data
%with firing rates to sacccades made in 8 directions.

% hAx = [];
% datVal = {[10*rand(1,8)+1] [8*rand(1,8)+1]};
% datErr = {[rand(1,8)] [rand(1,8)]};
% anglesToTest = linspace(0,2*pi,9);
% anglesToTest = anglesToTest(1:8); %8 evenly spaced going counterclockwise

if isempty(hAx)
    hF3 = figure;
    hAx = axes('parent',hF3);
end
cc = [0 1 0; 1 0 0]; %TODO make input arg with parser


%find how many data sets we have and ensure all are column vectors
if iscell(datVal)
    
    numDataSets = numel(datVal);
    maxWidth = max(max([datVal{:}]));
    datVal = cellfun(@(x) x(:),datVal,'uniformoutput',false);
    datErr = cellfun(@(x) x(:),datErr,'uniformoutput',false);
else
    maxWidth = max(datVal);
    numDataSets = 1;
    datVal = datVal(:);
    datVal = {datVal};
    datErr = datErr(:);
    datErr = {datErr};
    
end
anglesToTest = anglesToTest(:);

%if user specifies limit then overwrite. 
if ~isempty(plotLims)
    
    maxWidth = plotLims(2);
end

% build the background of the polar plot
% draw a unit circle
zz = exp(1i*linspace(0, 2*pi, 101)) * maxWidth;
line(real(zz),imag(zz),'color','k','linestyle',':','parent',hAx)
%also add marker line at half max size
zz = exp(1i*linspace(0, 2*pi, 101)) * maxWidth/2;
line(real(zz),imag(zz),'color','k','linestyle',':','parent',hAx)
%add crosshairs
line([-maxWidth maxWidth], [0 0],'color','k','linestyle','-','parent',hAx)
line([0 0], [-maxWidth maxWidth],'color','k','linestyle','-','parent',hAx)
line([-maxWidth maxWidth]*cos(pi/4), [maxWidth -maxWidth]*cos(pi/4),'color','k','linestyle','-','parent',hAx)
line([maxWidth -maxWidth]*cos(pi/4), [maxWidth -maxWidth]*cos(pi/4),'color','k','linestyle','-','parent',hAx)
set(hAx,'xtick',[])
set(hAx,'ytick',[])

%now add some axis labels
text(0,0,num2str(0))
text(maxWidth, 0, num2str(maxWidth));
text(maxWidth/2, 0, num2str(maxWidth/2));

hDat = cell(4,numDataSets);
for dataSetNum = 1:numDataSets
    
    %plot curve
    [X Y] =pol2cart([anglesToTest; anglesToTest(1)], [datVal{dataSetNum}; datVal{dataSetNum}(1)]);
    hDat{1,dataSetNum} = line(X, Y,'color',cc(dataSetNum,:),'linewidth',2,'parent',hAx)
    
    %now for confidence limit patch
    
    lowLim = [datVal{dataSetNum}; datVal{dataSetNum}(1)] - [datErr{dataSetNum}; datErr{dataSetNum}(1)];
    uppLim = [datVal{dataSetNum}; datVal{dataSetNum}(1)] + [datErr{dataSetNum}; datErr{dataSetNum}(1)];
    
    [lowX,lowY] = pol2cart([anglesToTest; anglesToTest(1)], lowLim);
    [uppX,uppY] = pol2cart([anglesToTest; anglesToTest(1)], uppLim);
    %patch( ,  , ,,'parent',hAx
    
    %hDat{2,dataSetNum} = line(lowX, lowY,'color',cc(dataSetNum,:),'linewidth',2,'parent',hAx)
    %hDat{3,dataSetNum} = line(uppX, uppY,'color',cc(dataSetNum,:),'linewidth',2,'parent',hAx)
    
    hDat{4,dataSetNum} = patch([uppX; flipud(lowX)],[uppY; flipud(lowY)],cc(dataSetNum,:),...
        'facealpha',0.5 ,'parent',hAx,'edgecolor',cc(dataSetNum,:),'edgealpha',0.5);
    
end

%uistack([hDat{:}])