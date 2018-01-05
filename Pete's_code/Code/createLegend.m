function hL = createLegend(legColours,legStrings,legStyles,parentAxis,posInAxis,varargin)
%my function for creating a legend using easily defined inputs, for now it
%takes a list of colours (RGB) and strings. It then plots a simple box with
%those colours next to those strings and places it at posInAxis(normlised
%units) within parentAxis
% 
% legColours = [1 0 0;0 1 0;0 0 0];
% legStrings = {'An','Ca','Sh'};
% posInAxis = [0.75 0.75 0.2 0.2];
% parentAxis  =gca;
p = inputParser;
p.addParamValue('textFontSize',10,@(x) isnumeric(x));
p.addParamValue('legLineWidth',2,@(x) isnumeric(x));
p.parse(varargin{:});
textFontSize = p.Results.textFontSize;
legLineWidth = p.Results.legLineWidth;
set(parentAxis,'units','normalized')
parPos = get(parentAxis,'position')

legPos(1) = posInAxis(1)*parPos(3)+parPos(1); %start in x
legPos(2) = posInAxis(2)*parPos(4)+parPos(2); %start in y
legPos(3) = posInAxis(3)*parPos(3); %size in x
legPos(4) = posInAxis(4)*parPos(4); %size in y

hL = axes('parent',get(parentAxis,'parent'),'units','normalized','position',legPos);

%numEntries = size(legColours,1);
numEntries = size(legStyles,1); %switched to stop errors when tryign to use same colourmap on multiple graphs
set(hL,'xlim',[0 4],'ylim',[0 numEntries+1])
for entryNum = 1:numEntries
    line([0.5 1.5],[entryNum entryNum],'color',legColours(entryNum,:),...
        'parent',hL,'linewidth',legLineWidth,'linestyle',legStyles{entryNum})
    text(2.5,entryNum,legStrings{entryNum},'parent',hL,'fontSize',textFontSize)
    %TODO add more lineprops 
end


set(hL,'visible','off')

