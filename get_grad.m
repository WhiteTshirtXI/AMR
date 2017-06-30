function G = get_grad( G, glev, parms )


%Get contribution to grad operator on gridlevel glev

%parameters
m = parms.m; n = parms.n; mg = parms.mg;
nrows = parms.nrowsG; ncols = parms.ncolsG;

%1st gridlevel is different from others
if glev == 1

    %--Block corresponding to x-velocities
    %  (takes pressures to left and right)

        %pressures to left of x-velocity point
        indvelx = repmat(1 : m-1,[1,n]);
        indvely = repelem(1 : n, m-1);
        velx_ind = get_velx_ind(indvelx,indvely, 1, parms);

        indpressx = repmat(1 : m-1,[1,n]);
        indpressy = repelem(1 : n, m-1);
        press_ind = get_press_ind(indpressx,indpressy, 1, parms);

        G = G - sparse( velx_ind, press_ind, ...
            ones(size(press_ind)), nrows, ncols);

        %pressures to right of x-velocity point
        indpressx = repmat(2 : m,[1,n]);
        indpressy = repelem(1 : n, m-1);
        press_ind = get_press_ind(indpressx,indpressy, 1, parms);

        G = G + sparse( velx_ind, press_ind, ...
            ones(size(press_ind)), nrows, ncols);
        
    %--Block corresponding to y-velocities
    %  (takes pressures above and below)

        %index for yvels starts after xvels
        nadd = get_velx_ind( m-1, n, mg, parms );
    
        %pressures below y-velocity point
        indvelx = repmat(1 : m,[1,n-1]);
        indvely = repelem(1 : n-1, m);
        vely_ind = get_vely_ind(indvelx,indvely, 1, parms);

        indpressx = repmat(1 : m,[1,n-1]);
        indpressy = repelem(1 : n-1, m);
        press_ind = get_press_ind(indpressx,indpressy, 1, parms);

        G = G - sparse( vely_ind + nadd, press_ind, ...
            ones(size(press_ind)), nrows, ncols);

        %pressures above y-velocity point
        indpressx = repmat(1 : m,[1,n-1]);
        indpressy = repelem(2 : n, m);
        press_ind = get_press_ind(indpressx,indpressy, 1, parms);

        G = G + sparse( vely_ind + nadd, press_ind, ...
            ones(size(press_ind)), nrows, ncols);
else
    
    
    %--Block corresponding to x-velocities
    %  (takes pressures to left and right)

        %--pressures to left of x-velocity point
        
        %Velocity indices
        
           %lower part (no overlap) 
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely, glev, parms);

            %middle part (overlap)
            indvelx = repmat([1 : m/4, 3*m/4+1:m-1],[1,n/2]);
            indvely = repelem(n/4 + 1 : 3*n/4, m/4 + (m/4-1) );
            velx_ind = [velx_ind, ...
                get_velx_ind(indvelx,indvely, glev, parms)];    
            
            %top part (no overlap)
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(3*n/4 + 1 : n, m-1);
            velx_ind = [velx_ind, ...
                get_velx_ind(indvelx,indvely, glev, parms)];
            
        %Pressure indices
        
            %lower part (no overlap)
            indpressx = repmat(1 : m-1,[1,n/4]);
            indpressy = repelem(1 : n/4, m-1);
            press_ind = get_press_ind(indpressx,indpressy, glev, parms);

            %middle part (overlap)
            indpressx = repmat([1 : m/4, 3*m/4+1:m-1],[1,n/2]);
            indpressy = repelem(n/4 + 1 : 3*n/4, m/2 - 1);
            press_ind = [press_ind, ...
                get_press_ind(indpressx,indpressy, glev, parms)];
            
            %top part (no overlap)
            indpressx = repmat(1 : m-1,[1,n/4]);
            indpressy = repelem(3*n/4 + 1 : n, m-1);
            press_ind = [press_ind, ...
                get_press_ind(indpressx,indpressy, glev, parms)];
            
        %Add to grad matrix    
        G = G - sparse( velx_ind, press_ind, ...
            ones(size(press_ind)), nrows, ncols);

        %--
        
        %--pressures to right of x-velocity point
        
         %Velocity indices
        
           %lower part (no overlap) 
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(1 : n/4, m-1);
            velx_ind = get_velx_ind(indvelx,indvely, glev, parms);

            %middle part (overlap)
            indvelx = repmat([1:m/4-1, 3*m/4 : m-1],[1,n/2]);
            indvely = repelem(n/4 + 1 : 3*n/4, m/2 -1);
            velx_ind = [velx_ind, ...
                get_velx_ind(indvelx,indvely, glev, parms)];    
            
            %top part (no overlap)
            indvelx = repmat(1 : m-1,[1,n/4]);
            indvely = repelem(3*n/4 + 1 : n, m-1);
            velx_ind = [velx_ind, ...
                get_velx_ind(indvelx,indvely, glev, parms)];
            
        %Pressure indices
        
            %lower part (no overlap)
            indpressx = repmat(2 : m,[1,n/4]);
            indpressy = repelem(1 : n/4, m-1);
            press_ind = get_press_ind(indpressx,indpressy, glev, parms);

            %middle part (overlap)
            indpressx = repmat([2: m/4, 3*m/4+1 : m],[1,n/2]);
            indpressy = repelem(n/4 + 1 : 3*n/4, m/2 -1);
            press_ind = [press_ind, ...
                get_press_ind(indpressx,indpressy, glev, parms)];
            
            %top part (no overlap)
            indpressx = repmat(2 : m,[1,n/4]);
            indpressy = repelem(3*n/4 + 1 : n, m-1);
            press_ind = [press_ind, ...
                get_press_ind(indpressx,indpressy, glev, parms)];
            
        %Add to grad matrix    
        G = G + sparse( velx_ind, press_ind, ...
            ones(size(press_ind)), nrows, ncols);
    
        %--
    %--
    
    
    %--Block corresponding to y-velocities
    %  (takes pressures above and below)

        %--pressures below y-velocity point
        
        %Velocity indices
        
            %index for yvels starts after xvels
            nadd = get_velx_ind( m-1, n, mg, parms );
        
           %lower part (no overlap) 
            indvelx = repmat(1 : m,[1,n/4]);
            indvely = repelem(1 : n/4, m);
            vely_ind = get_vely_ind(indvelx,indvely, glev, parms);

            %middle part (overlap)
            indvelx = repmat([1 : m/4, 3*m/4+1 : m],[1,n/2]);
            indvely = repelem(n/4 + 1 : 3*n/4, m/2);
            vely_ind = [vely_ind, ...
                get_vely_ind(indvelx,indvely, glev, parms)];    
            
            %top part (no overlap)
            indvelx = repmat(1 : m,[1,n/4-1]);
            indvely = repelem(3*n/4+1 : n-1, m);
            vely_ind = [vely_ind, ...
                get_vely_ind(indvelx,indvely, glev, parms)];
            
        %Pressure indices
        
            %lower part (no overlap)
            indpressx = repmat(1 : m,[1,n/4]);
            indpressy = repelem(1 : n/4, m);
            press_ind = get_press_ind(indpressx,indpressy, glev, parms);

            %middle part (overlap)
            indpressx = repmat([1 : m/4,3*m/4+1 : m],[1,n/2]);
            indpressy = repelem(n/4 + 1 : 3*n/4, m/2);
            press_ind = [press_ind, ...
                get_press_ind(indpressx,indpressy, glev, parms)];
            
            %top part (no overlap)
            indpressx = repmat(1 : m,[1,n/4-1]);
            indpressy = repelem(3*n/4+1 : n-1, m);
            press_ind = [press_ind, ...
                get_press_ind(indpressx,indpressy, glev, parms)];
            
        %Add to grad matrix    
        G = G - sparse( vely_ind + nadd, press_ind, ...
            ones(size(press_ind)), nrows, ncols);

        %--
        
        %--pressures above y-velocity point
        
         %Velocity indices
        
           %lower part (no overlap) 
            indvelx = repmat(1 : m,[1,n/4-1]);
            indvely = repelem(1 : n/4-1, m);
            vely_ind = get_vely_ind(indvelx,indvely, glev, parms);

            %middle part (overlap)
            indvelx = repmat([1 : m/4,3*m/4+1 : m],[1,n/2]);
            indvely = repelem(n/4 : 3*n/4-1, m/2);
            vely_ind = [vely_ind, ...
                get_vely_ind(indvelx,indvely, glev, parms)];    
            
            %top part (no overlap)
            indvelx = repmat(1 : m,[1,n/4]);
            indvely = repelem(3*n/4 : n-1, m);
            vely_ind = [vely_ind, ...
                get_vely_ind(indvelx,indvely, glev, parms)];
            
        %Pressure indices
        
            %lower part (no overlap)
            indpressx = repmat(1 : m,[1,n/4-1]);
            indpressy = repelem(2 : n/4, m);
            press_ind = get_press_ind(indpressx,indpressy, glev, parms);

            %middle part (overlap)
            indpressx = repmat([1 : m/4,3*m/4+1 : m],[1,n/2]);
            indpressy = repelem(n/4 + 1 : 3*n/4, m/2);
            press_ind = [press_ind, ...
                get_press_ind(indpressx,indpressy, glev, parms)];
            
            %top part (no overlap)
            indpressx = repmat(1 : m,[1,n/4]);
            indpressy = repelem(3*n/4+1 : n, m);
            press_ind = [press_ind, ...
                get_press_ind(indpressx,indpressy, glev, parms)];
            
        %Add to grad matrix    
        G = G + sparse( vely_ind + nadd, press_ind, ...
            ones(size(press_ind)), nrows, ncols);
    
        %--
    %--
end



