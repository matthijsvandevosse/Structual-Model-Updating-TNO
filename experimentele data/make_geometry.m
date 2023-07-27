g_beam = [];
for x = 0:size(sensor_loc,2)
    if x == 0
        g = [ 3;4;
            0;
            sensor_loc(x+1);
            sensor_loc(x+1); 
            0;
            h-t;
            h-t;
            h;
            h];
    elseif x == size(sensor_loc,2)
        g = [ 3;4;
            sensor_loc(x);
            l;
            l;
            sensor_loc(x);
            h-t;
            h-t;
            h;
            h];
    else
        g = [ 3;4;
            sensor_loc(x);
            sensor_loc(x+1);
            sensor_loc(x+1); 
            sensor_loc(x);
            h-t;
            h-t;
            h;
            h];
    end
        g_beam = [g_beam, g];
end
g_sides = [ 3 3;
            4 4;
            0 l-t;
            t l;
            t l; 
            0 l-t;
            0 0;
            0 0;
            h-t h-t;
            h-t h-t];

g_sides = [ 3 3;
            4 4;
            0 l-t;
            t l;
            t l; 
            0 l-t;
            0 0;
            0 0;
            h-t h-t;
            h-t h-t];

g_all = [g_beam g_sides];
% g_all = g_beam;

gdescg = decsg(g_all);
geometry = geometryFromEdges(structuralmodelinit,gdescg);