function D = get_div_BCs( D, glev, parms )

%get BCs for divergence operator (requires info from coarser grid 
%    for current grid)

m = parms.m; n = parms.n; mg = parms.mg;
nrows = parms.nrowsD; ncols = parms.ncolsD;

scl = 1/2; %multiply by 1/2 because we are dealing with vel fluxes
% scl = 1;

%--Left edge (requires x velocities)

    %Pressure indices
        indpressx = repmat( 1, [1, n] );
        indpressy = repelem( 1 : n, 1 );
        press_ind = get_press_ind( indpressx, indpressy, glev, parms );
        
    %Velocity indices on coarser grid
        indvelx = repmat( m/4, [1, n/2 + 2] );
        indvely = repelem( n/4 : 3 * n/4 + 1, 1 );
        velx_ind = get_velx_ind( indvelx, indvely, glev + 1, parms );
        
    %Add to div matrix: closest vels
        D = D - 3/4 * scl * sparse( press_ind( 1 : 2 : n-1), velx_ind(2:end-1), ...
            ones(size(velx_ind(2:end-1)) ), nrows, ncols );
        D = D - 3/4 * scl * sparse( press_ind( 2 : 2 : n ), velx_ind(2:end-1), ...
            ones(size(velx_ind(2:end-1)) ), nrows, ncols );
        
    %Add to div matrix: further vels
        D = D - 1/4 * scl * sparse( press_ind( 1 : 2 : n-1), velx_ind(1:end-2), ...
            ones(size(velx_ind(1:end-2)) ), nrows, ncols );
        D = D - 1/4 * scl * sparse( press_ind( 2 : 2 : n ), velx_ind(3:end), ...
            ones(size(velx_ind(3:end)) ), nrows, ncols );
%--
        
%--Right edge (requires x-velocities)

    %Pressure indices
        indpressx = repmat( m, [1, n] );
        indpressy = repelem( 1 : n, 1 );
        press_ind = get_press_ind( indpressx, indpressy, glev, parms );
        
    %Velocity indices on coarser grid
        indvelx = repmat( 3*m/4, [1, n/2 + 2] );
        indvely = repelem( n/4 : 3 * n/4 + 1, 1 );
        velx_ind = get_velx_ind( indvelx, indvely, glev + 1, parms );
        
    %Add to div matrix: closest vels
        D = D + 3/4 * scl * sparse( press_ind( 1 : 2 : n-1), velx_ind(2:end-1), ...
            ones(size(velx_ind(2:end-1)) ), nrows, ncols );
        D = D + 3/4 * scl * sparse( press_ind( 2 : 2 : n ), velx_ind(2:end-1), ...
            ones(size(velx_ind(2:end-1)) ), nrows, ncols );
        
    %Add to div matrix: further vels
        D = D + 1/4 * scl * sparse( press_ind( 1 : 2 : n-1), velx_ind(1:end-2), ...
            ones(size(velx_ind(1:end-2)) ), nrows, ncols );
        D = D + 1/4 * scl * sparse( press_ind( 2 : 2 : n ), velx_ind(3:end), ...
            ones(size(velx_ind(3:end)) ), nrows, ncols );

%--

%--Bottom edge (requires y-velocities)

    %Pressure indices
        indpressx = repmat( 1 : m, 1 );
        indpressy = repelem( 1 , m );
        press_ind = get_press_ind( indpressx, indpressy, glev, parms );
        
    %Velocity indices on coarser grid
    
        %index for yvels starts after xvels
        nadd = get_velx_ind( m-1, n, mg, parms );
        
        indvelx = repmat( m/4 : 3*m/4 + 1, 1 );
        indvely = repelem( n/4 , m/2 + 2 );
        vely_ind = get_vely_ind( indvelx, indvely, glev + 1, ...
            parms ) + nadd;
        
    %Add to div matrix: closest vels
        D = D - scl * 3/4 * sparse( press_ind( 1 : 2 : m-1), vely_ind(2:end-1), ...
            ones(size(velx_ind(2:end-1)) ), nrows, ncols );
        D = D - scl * 3/4 * sparse( press_ind( 2 : 2 : m ), vely_ind(2:end-1), ...
            ones(size(velx_ind(2:end-1)) ), nrows, ncols );
        
    %Add to div matrix: further vels
        D = D - scl * 1/4 * sparse( press_ind( 1 : 2 : m-1), vely_ind(1:end-2), ...
            ones(size(velx_ind(2:end-1)) ), nrows, ncols );
        D = D - scl * 1/4 * sparse( press_ind( 2 : 2 : m ), vely_ind(3:end), ...
            ones(size(velx_ind(2:end-1)) ), nrows, ncols );

%--

%--Top edge (requires y-velocities)

    %Pressure indices
        indpressx = repmat( 1 : m, 1 );
        indpressy = repelem( n , m );
        press_ind = get_press_ind( indpressx, indpressy, glev, parms );
        
    %Velocity indices on coarser grid
    
        indvelx = repmat( m/4 : 3*m/4 + 1, 1 );
        indvely = repelem( 3*n/4 , m/2 + 2 );
        vely_ind = get_vely_ind( indvelx, indvely, glev + 1, ...
            parms ) + nadd;
        
    %Add to div matrix: closest vels
        D = D + scl * 3/4 * sparse( press_ind( 1 : 2 : m-1), vely_ind(2:end-1), ...
            ones(size(velx_ind(2:end-1)) ), nrows, ncols );
        D = D + scl * 3/4 * sparse( press_ind( 2 : 2 : m ), vely_ind(2:end-1), ...
            ones(size(velx_ind(2:end-1)) ), nrows, ncols );
        
   %Add to div matrix: further vels
        D = D + scl * 1/4 * sparse( press_ind( 1 : 2 : m-1), vely_ind(1:end-2), ...
            ones(size(velx_ind(2:end-1)) ), nrows, ncols );
        D = D + scl * 1/4 * sparse( press_ind( 2 : 2 : m ), vely_ind(3:end), ...
            ones(size(velx_ind(2:end-1)) ), nrows, ncols );
 

%--



