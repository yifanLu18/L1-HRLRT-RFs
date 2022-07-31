% Loop to calculate the coordinates of the hypocenter that satisfy the specified azimuth 
% and epicentral distance (currently mainly used for forward modeling data)
% Written by Yifan Lu, more details can be found in https://doi.org/10.1093/gji/ggac260

function  [event_la,event_lo] = cal_event_location(sta_la,sta_lo,deg,back_az,error,print_screen)

for i = 1:length(back_az)
    for eve_la = -180:.1:180
        for eve_lo = -90:.1:90
            
            [~,baz] = calbaz(eve_la,eve_lo,sta_la,sta_lo);
            [~,dis2] = caldisdeg(eve_la,eve_lo,sta_la,sta_lo);
            if ( abs(baz-back_az(i))<error && abs(deg-dis2)<error )
                if (print_screen == 1)
                    fprintf('eve_la=%.2f eve_lo=%.2f ------ baz=%.2f deg=%.2f\n',eve_la,eve_lo,baz,dis2);
                    fprintf("-------------------------------------------------------\n");
                end
                event_la = eve_la;
                event_lo = eve_lo;
                return;
            end
        end
    end
end

