% Tilt-interface-converted PS-wave travel time difference
% Written by Yifan Lu, more details can be found in https://doi.org/10.1093/gji/ggac260

function tps = cal_tps(stla,stlo,evla,evlo,event_depth,vp1,vp2,vs1,vs2,alpha,xita,h)

[~,baz] = calbaz(evla,evlo,stla,stlo);
[~,deg] = caldisdeg(evla,evlo,stla,stlo);
P = taupTime([],event_depth,'P','deg',deg);
p = P.rayParam / 6371;

n = [sind(xita)*cosd(alpha),sind(xita)*sind(alpha),-cosd(xita)];
P1 = [-p*cosd(baz),-p*sind(baz),-sqrt(1/vp2/vp2-p*p)];

slowness_p = sqrt(1/vp2/vp2-p*p) - cosd(xita)*(dot(n,P1)) + cosd(xita) * sqrt(1/vp1/vp1 - 1/vp2/vp2 + (dot(n,P1))^2 );

slowness_s = sqrt(1/vp2/vp2-p*p) - cosd(xita)*(dot(n,P1)) + cosd(xita) * sqrt(1/vs1/vs1 - 1/vp2/vp2 + (dot(n,P1))^2 );

tps = h*(slowness_s - slowness_p);

