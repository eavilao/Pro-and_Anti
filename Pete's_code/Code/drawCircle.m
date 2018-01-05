function hC = drawCircle(x,y,r,hA)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%hA is axis handle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
hC = line(x+xp,y+yp,'parent',hA);
end