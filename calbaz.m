% Calculate azimuth and inverse azimuth

function [az,baz]=calbaz(evla,evlo,stla,stlo)
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
lo=stlo/180.0*pi;
la=stla/180.0*pi;
z=r0*sin(la);
xy=r0*cos(la);
x=xy*cos(lo);
y=xy*sin(lo);
stx=x;
sty=y;
stz=z;
a=acos((evx*stx+evy*sty+evz*stz)/(6371.0^2));
b=(90.0-evla)/180.0*pi;
c=(90.0-stla)/180.0*pi;
p=(a+b+c)/2.0;
da=2*acos(sqrt(sin(p)*sin(p-a)/sin(b)/sin(c)));
db=2*acos(sqrt(sin(p)*sin(p-b)/sin(a)/sin(c)));
dc=2*acos(sqrt(sin(p)*sin(p-c)/sin(b)/sin(a)));

if (evlo-stlo>=0&evlo-stlo<=180)|(evlo+360-stlo>=0&evlo+360-stlo<=180)
    baz=db/pi*180;
    az=360-dc/pi*180;
else
    baz=360-db/pi*180;
    az=dc/pi*180;
end
