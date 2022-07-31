% Calculate the epicenter distance (km, degree)

function [dis1,dis2]=caldisdeg(evla,evlo,stla,stlo)
r0=6371.0;
lo=evlo/180.0*pi;
la=evla/180.0*pi;
z=r0*sin(la);
xy=r0*cos(la);
x=xy*cos(lo);
y=xy*sin(lo);
evx=x;
evy=y;
evz=z;
% for i=1:lev
lo=stlo/180.0*pi;
la=stla/180.0*pi;
z=r0*sin(la);
xy=r0*cos(la);
x=xy*cos(lo);
y=xy*sin(lo);
stx=x;
sty=y;
stz=z;

r=sqrt((evx-stx)^2+(evy-sty)^2+(evz-stz)^2);
dis1=2.0*asin(r/2/r0)*r0; % km
dis2=2.0*asin(r/2/r0)/pi*180; % degree