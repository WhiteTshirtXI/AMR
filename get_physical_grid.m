function [x, y, parms] = get_physical_grid( parms )

%Build physical grid of x and y points from parameters in parms (e.g. # of
%desired points, length of domain, and desired offset)

m = parms.m; n = parms.n; mg = parms.mg; len = parms.len; 
offsetx = parms.offsetx; offsety = parms.offsety;

for lev = 1 : mg
    
    fac = 2^(lev-1); 
    
    % Grid spacing in both directions for current grid
    delta = len ./ m *fac;              
    
    % Offset in x direction for current grid
    offx = 2^(lev-1) * len/2 - len/2 + offsetx ;    
    
    % Offset in y direction for current grid
    offy = 2^(lev-1) * (n*len/m)/2 - (n*len/m)/2 + offsety ;
    
    %1st grid level is easy
    if lev == 1
    
        %x and y points on streamfunction / circulation grid (vertices)
        x_gamm = delta : delta : (m-1) * delta;
        y_gamm = delta : delta : (n-1) * delta;

        x_gamm = x_gamm - offx;
        y_gamm = y_gamm - offy;
        
        x.gamm = repmat( x_gamm, [1, length(y_gamm)] )';
        y.gamm = repelem( y_gamm, length(x_gamm) )';
        
        %x and y points on pressure grid (centers)
        x_press = delta/2 : delta : (2*m-1) * delta/2;
        y_press = delta/2 : delta : (2*n-1) * delta/2;

        x_press = x_press - offx;
        y_press = y_press - offy;
        
        x.press = repmat( x_press, [1, length(y_press)] )';
        y.press = repelem( y_press, length(x_press) )';
        
        %x and y points on xvel grid (vertical edges)
        x_velx = delta : delta : (m-1) * delta;
        y_velx = delta/2 : delta : (2*n-1) * delta/2;

        x_velx = x_velx - offx;
        y_velx = y_velx - offy;
        
        x.velx = repmat( x_velx, [1, length(y_velx)] )';
        y.velx = repelem( y_velx, length(x_velx) )';
        
        %x and y points on yvel grid (horizontal edges)
        x_vely = delta/2 : delta : (2*m-1) * delta/2;
        y_vely = delta : delta : (n-1) * delta;

        x_vely = x_vely - offx;
        y_vely = y_vely - offy;
        
        x.vely = repmat( x_vely, [1, length(y_vely)] )';
        y.vely = repelem( y_vely, length(x_vely) )';
                
    else
        
        %x and y points on streamfunction / circulation grid (vertices)
        x_gamm = delta : delta : (m-1) * delta;
        y_gamm = delta : delta : (n-1) * delta;

        x_gamm = x_gamm - offx;
        y_gamm = y_gamm - offy;
        
        indgammx = repmat(1 : m-1, [1,n-1]); 
        indgammy = repelem(1 : n-1, m-1 );
        
        gamm_ind = get_vort_ind( indgammx, indgammy, lev, parms );
        
        x_gammb = repmat( x_gamm, [1, length(y_gamm)] )';
        y_gammb = repelem( y_gamm, length(x_gamm) )';
        
        x.gamm = [x.gamm; x_gammb( gamm_ind ~= 0 )];
        y.gamm = [y.gamm; y_gammb( gamm_ind ~= 0 )];
        
        %x and y points on pressure grid (centers)
        x_press = delta/2 : delta : (2*m-1) * delta/2;
        y_press = delta/2 : delta : (2*n-1) * delta/2;

        x_press = x_press - offx;
        y_press = y_press - offy;
        
        x_pressb = repmat( x_press, [1, length(y_press)] )';
        y_pressb = repelem( y_press, length(x_press) )';
        
        indpressx = repmat(1 : m, [1,n]); 
        indpressy = repelem(1 : n, m );
        
        press_ind = get_press_ind( indpressx, indpressy, lev, parms );
        
        x.press = [x.press; x_pressb( press_ind ~= 0 )];
        y.press = [y.press; y_pressb( press_ind ~= 0 )];
        
        %x and y points on xvel grid (vertical edges)
        x_velx = delta : delta : (m-1) * delta;
        y_velx = delta/2 : delta : (2*n-1) * delta/2;

        x_velx = x_velx - offx;
        y_velx = y_velx - offy;
        
        x_velxb = repmat( x_velx, [1, length(y_velx)] )';
        y_velxb = repelem( y_velx, length(x_velx) )';
        
        indvelx = repmat(1 : m-1, [1,n]); 
        indvely = repelem(1 : n, m-1 );
        
        velx_ind = get_velx_ind( indvelx, indvely, lev, parms );
        
        x.velx = [x.velx; x_velxb( velx_ind ~= 0 )];
        y.velx = [y.velx; y_velxb( velx_ind ~= 0 )];
        
        %x and y points on yvel grid (horizontal edges)
        x_vely = delta/2 : delta : (2*m-1) * delta/2;
        y_vely = delta : delta : (n-1) * delta;

        x_vely = x_vely - offx;
        y_vely = y_vely - offy;
        
        x_velyb = repmat( x_vely, [1, length(y_vely)] )';
        y_velyb = repelem( y_vely, length(x_vely) )';
        
        indvelx = repmat(1 : m, [1,n-1]); 
        indvely = repelem(1 : n-1, m );
        
        vely_ind = get_vely_ind( indvelx, indvely, lev, parms );
        
        x.vely = [x.vely; x_velyb( vely_ind ~= 0 )];
        y.vely = [y.vely; y_velyb( vely_ind ~= 0 )];
        
                
    end
    
    parms.xmin = -offx;
    parms.ymin = -offy;
    
end

