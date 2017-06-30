function G = get_grad_BCs( G, glev, parms )

%get BCs for gradient operator (requires info from fine grid 
%    for current grid)

m = parms.m; n = parms.n; mg = parms.mg;
nrows = parms.nrowsG; ncols = parms.ncolsG;

%--Block corresponding to x-velocities
%  (takes pressures to left and right)

    %--left edge (takes pressures on right side)

        %Velocity indices 
            indvelx = repmat( m/4, [1,n/2] );
            indvely = repelem( n/4 + 1 : 3*n/4, 1 );
            velx_ind = get_velx_ind(indvelx,indvely, glev, parms);    


        %Pressure indices on finer grid
            %lower ones
            indpressx = repmat( 1, [1,n/2] );
            indpressy = repelem( 1 : 2 : n-1 , 1 );
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G + 1/4 * sparse( velx_ind, press_ind, ...
                ones(size(press_ind)), nrows, ncols);

            %upper ones
            indpressx = repmat( 1, [1,n/2] );
            indpressy = repelem( 2 : 2 : n , 1 );
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G + 1/4 * sparse( velx_ind, press_ind, ...
                ones(size(press_ind)), nrows, ncols);
            
            %lower ones
            indpressx = repmat( 2, [1,n/2] );
            indpressy = repelem( 1 : 2 : n-1 , 1 );
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G + 1/4 * sparse( velx_ind, press_ind, ...
                ones(size(press_ind)), nrows, ncols);

            %upper ones
            indpressx = repmat( 2, [1,n/2] );
            indpressy = repelem( 2 : 2 : n , 1 );
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G + 1/4 * sparse( velx_ind, press_ind, ...
                ones(size(press_ind)), nrows, ncols);

    %--

    %--right edge (takes pressures on left side)

        %Velocity indices 
            indvelx = repmat( 3*m/4, [1,n/2] );
            indvely = repelem( n/4 + 1 : 3*n/4, 1 );
            velx_ind = get_velx_ind(indvelx,indvely, glev, parms);    


        %Pressure indices on finer grid
            %lower ones
            indpressx = repmat( m, [1,n/2] );
            indpressy = repelem( 1 : 2 : n-1 , 1 );
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G - 1/4 * sparse( velx_ind, press_ind, ...
                ones(size(press_ind)), nrows, ncols);

            %upper ones
            indpressx = repmat( m, [1,n/2] );
            indpressy = repelem( 2 : 2 : n , 1 );
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G - 1/4 * sparse( velx_ind, press_ind, ...
                ones(size(press_ind)), nrows, ncols);
            
            %lower ones
            indpressx = repmat( m-1, [1,n/2] );
            indpressy = repelem( 1 : 2 : n-1 , 1 );
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G - 1/4 * sparse( velx_ind, press_ind, ...
                ones(size(press_ind)), nrows, ncols);

            %upper ones
            indpressx = repmat( m-1, [1,n/2] );
            indpressy = repelem( 2 : 2 : n , 1 );
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G - 1/4 * sparse( velx_ind, press_ind, ...
                ones(size(press_ind)), nrows, ncols);
    %--
%--


%--Block corresponding to y-velocities
%  (takes pressures above and below)

    %--top edge (requires pressures below)

        %Velocity indices
        
            %index for yvels starts after xvels
            nadd = get_velx_ind( m-1, n, mg, parms );

            %indices
            indvelx = repmat(m/4 + 1 : 3*m/4, 1);
            indvely = repelem(3*n/4 , m/2);
            vely_ind = get_vely_ind(indvelx,indvely, glev, parms);

        %Pressure indices

            %left part
            indpressx = repmat(1 : 2 : m-1, 1);
            indpressy = repelem(n, m/2);
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G - 1/4 * sparse( vely_ind + nadd, press_ind, ...
                ones(size(press_ind)), nrows, ncols);
            
            %right part
            indpressx = repmat(2 : 2 : m, 1);
            indpressy = repelem(n, m/2);
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G - 1/4 * sparse( vely_ind + nadd, press_ind, ...
                ones(size(press_ind)), nrows, ncols);
            
            %left part
            indpressx = repmat(1 : 2 : m-1, 1);
            indpressy = repelem(n-1, m/2);
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G - 1/4 * sparse( vely_ind + nadd, press_ind, ...
                ones(size(press_ind)), nrows, ncols);
            
            %right part
            indpressx = repmat(2 : 2 : m, 1);
            indpressy = repelem(n-1, m/2);
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G - 1/4 * sparse( vely_ind + nadd, press_ind, ...
                ones(size(press_ind)), nrows, ncols);

    %--

    %--bottom edge (requires pressures above)

     %Velocity indices
        
            %indices
            indvelx = repmat(m/4 + 1 : 3*m/4, 1);
            indvely = repelem(n/4 , m/2);
            vely_ind = get_vely_ind(indvelx,indvely, glev, parms);

        %Pressure indices

            %left part
            indpressx = repmat(1 : 2 : m-1, 1);
            indpressy = repelem(1, m/2);
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G + 1/4 * sparse( vely_ind + nadd, press_ind, ...
                ones(size(press_ind)), nrows, ncols);
            
            %right part
            indpressx = repmat(2 : 2 : m, 1);
            indpressy = repelem(1, m/2);
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G + 1/4 * sparse( vely_ind + nadd, press_ind, ...
                ones(size(press_ind)), nrows, ncols);
            
            %left part
            indpressx = repmat(1 : 2 : m-1, 1);
            indpressy = repelem(2, m/2);
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G + 1/4 * sparse( vely_ind + nadd, press_ind, ...
                ones(size(press_ind)), nrows, ncols);
            
            %right part
            indpressx = repmat(2 : 2 : m, 1);
            indpressy = repelem(2, m/2);
            press_ind = get_press_ind(indpressx,indpressy, glev-1, parms);

            %Add to grad matrix    
            G = G + 1/4 * sparse( vely_ind + nadd, press_ind, ...
                ones(size(press_ind)), nrows, ncols);
    %--
%--



