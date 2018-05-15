function rosePlotProAnti(data,lowerBound,upperBound )
% data in N x 8 array format for directions from 0 to 315 degrees
% CIl = lower conficene interval CI u is upper CI is confidence interval in
% N x 8
% mw specifies the radius of the plotted polar plot
            mw = ceil(max(max(upperBound))*10)/10;

% pervent undershooting zero
    lowerBound(find(lowerBound<=0))=0;

%  coordinates
ori = [0 22.5 45 67.5 90 112.5 135 157.5];
ori = circ_axial(circ_ang2rad(ori),2);

% spacing of bins
dori = diff(ori(1:2));

 clear i
       figure; hold on;

cat=size(data,1); % the amount of conditions to plot(pro/anti)
for g=1:cat
    clearvars x y x2 y2 x3  y3 w
    w = data(g,:);
    ci_min = lowerBound(g,:);
    ci_plus  =  upperBound(g,:);

 

    r = circ_r(ori,w(1,:),dori) * mw;
    r(find(isnan(r)==1))=0;
    phi = circ_mean(ori,w(1,:));
    zm = r*exp(i*phi');
    
    
    
    if g==1
        h1 =  polar([ori ori(1)], [w(1,:) w(1,1)]);
        set(h1,'color','r','linewidth',2);
        
        h1 = polar([ori ori(1)], [w(1,:) w(1,1)]);
        h1.Marker='none';
        h1.LineStyle = '-';
        set(h1,'color','r','linewidth',2);
        
        h11= polar([ori ori(1)], [w(1,:) w(1,1)]);
        h11.Marker='.';
        h11.LineStyle = 'none';
        set(h11,'color','r','linewidth',1);

        
        
        [x,y]=pol2cart(ori,ci_plus(1,:)); % polar to cartesian)
        [x2,y2]=pol2cart(ori,ci_min(1,:));
        
        % trow out nans - patch has to be edited in illustartor 
        if any(isnan(ci_plus))
            ind= isnan(ci_plus);
            ci_plus(ind)=[];
            ci_min(ind)=[];
            ori_reduced=ori(~ind);
            [x,y]=pol2cart(ori_reduced,ci_plus(1,:)); % polar to cartesian)
        	[x2,y2]=pol2cart(ori_reduced,ci_min(1,:));
             
        end
        
        %patch std
        
        x3=[x,x(1),x2(1),fliplr(x2)];
        y3=[y,y(1),y2(1),fliplr(y2)];
        
        
        
        
        patch(x3,y3,[.5 0 0 ],'EdgeColor','none','FaceAlpha',0.3) % blue patch for non lesion
        
    elseif g==2
         
        h2 = polar([ori ori(1)], [w(1,:) w(1,1)]);
            h2.Marker='none';
                h2.LineStyle = '-';
                set(h2,'color','g','linewidth',2);
                
                h3= polar([ori ori(1)], [w(1,:) w(1,1)]); 
                h3.Marker='.';
                h3.LineStyle = 'none';
                set(h3,'color','g','linewidth',1);

        [x,y]=pol2cart(ori,ci_plus(1,:)); % polar to cartesian)
        [x2,y2]=pol2cart(ori,ci_min(1,:));
        if any(isnan(ci_plus))
            ind= isnan(ci_plus);
            ci_plus(ind)=[];
            ci_min(ind)=[];
            ori_reduced=ori(~ind);
            [x,y]=pol2cart(ori_reduced,ci_plus(1,:)); % polar to cartesian)
        	[x2,y2]=pol2cart(ori_reduced,ci_min(1,:));
             
        end
        %patch std
        x3=[x,x(1),x2(1),fliplr(x2)];
        y3=[y,y(1),y2(1),fliplr(y2)];
        
     
      
        patch(x3,y3,[ 0 .5 0 ],'EdgeColor','none','FaceAlpha',0.3); % red patch for   lesion
        
   
end
% draw a unit circle
 
zz = exp(1i*linspace(0, 2*pi, 101)) * mw;
plot(real(zz),imag(zz),'k:')
plot([-mw mw], [0 0], 'k:', [0 0], [-mw mw], 'k:')


formatSubplot(gca,'ax','square','box','off','lim',[-mw mw -mw mw])
formatSubplot(gca, 'xTickNumber',3,'yTickNumber',3)
axis square


 end