% function [hscatter,hbar,ax]=scatterDiagHist(x,y,dedge,scattervarargin)
% scatterDiagHist(x,y)
% scatterDiagHist(x,y,dedge)
% scatterDiagHist(x,y,dedge,scattervarargin)
% [hscatter,hbar]=scatterDiagHist(...)
% INPUT
% x,y: vectors of x and y data. same size
% dedge: edge for x,y difference histogram, a vector or an integer(no. of bin).
%           use [] if dont want to set it.
% scattervargin: input variables same as 'scatter' function
% OUTPUT
% hscatter: handle for scatter points
% hbar: handle for histogram bar
% ax: axes
% Kefei Liu. 2017.5.5
function [h1,h2,a1]=scatterDiagHist(x,y,dedge,varargin)
d=x-y;
rlim=[min([x(:);y(:)]),max([x(:);y(:)])];rlim(2)=rlim(2)+rlim(2)-rlim(1);
diflim=sqrt(2)/2* ((rlim(2)-rlim(1))/2);
if nargin<3 || isempty(dedge)
dedge=linspace(-diflim,diflim,13);
end
if length(dedge)==1
    dedge=linspace(-diflim,diflim,dedge+1);
end
dx=dedge(1:end-1)+(dedge(2)-dedge(1))/2;
dn=histcounts(d,dedge);
% s=sqrt(2);% scaling factor
a1=gca;
h1=scatter(x,y,varargin{:});
hold on;
diagxy=[rlim(1),rlim(2)-(rlim(2)-rlim(1))/2];
plot(diagxy,diagxy,'k:');
a2=axes('position',get(a1,'position'));
h2=bar(dx,dn);
hold on;plot([0,0],[0,max(dn)*1.5],'k:');
pos1=get(a1,'Position');posfig=get(gcf,'Position');
tmp=min(pos1([3,4]).*posfig([3,4]));
pos1(3)=tmp/posfig(3);pos1(4)=tmp/posfig(4);
pos2=[pos1(1)+(pos1(3))/2-pos1(3)/6,pos1(2)+(pos1(4))/2-pos1(4)/6,pos1(3),pos1(4)];
% axis(a1, 'square');
set(a1,'ActivePositionProperty','position','position',pos1,'xlim',rlim,...
    'ylim',rlim,'color','none');
set(a2,'Position',pos2,'ActivePositionProperty','position','color','none');
    set(a2,'View',[45,90],'box','off','ycolor','w','xlim',[-1,1]*diflim,...
        'xtick',[],'Nextplot','add');
% axis(a2, 'square');
% % uistack(a2,'top');
% axes(a1);
% linkaxes([a1,a2],'xy');
set(gcf,'CurrentAxes',a1)