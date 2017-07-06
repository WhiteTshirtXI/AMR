clear all, close all, clc

%Test the curl operators. 

mvect = [20 40 80 160];

errinfC = zeros( length( mvect), 1);
errinfR = errinfC;

for j = 1 : length( mvect )

    %# of points in x and y directions
    parms.m = mvect(j); 
    parms.n = mvect(j);

    %length of domain for finest grid
    parms.len = 2;

    %# of grid levels
    parms.mg = 2;

    %Offset -- bottom left corner of finest grid is (-offsetx,-offsety)
    parms.offsetx = 0;
    parms.offsety = 0;

    %get points on physical domain
    [x, y, parms] = get_physical_grid( parms );
    
    %Get gradient and divergence operators
    [C, R] = get_curls( parms );
    [G, D] = get_divgrad( parms );
    
    %Scale by grid spacing
    M_vel = get_hscale_vel( parms );
    
%     errDC =  max(max(abs( D * M_vel * C ))) 
%     errRG = max( max( abs( R * G ) ) )
%     
   

    %error for C
        
        %Scale by grid spacing
        
        C = M_vel * C;
        
        %function to take curl of
        
        %Get interval for making sins periodic:
        l_sc = parms.mg * parms.len;
        shift = parms.xmin;
        
        s = sin( 3* (x.gamm-shift) * (2*pi/ l_sc)  ) .* sin( 4 *...
            (y.gamm-shift) * (2*pi/ l_sc) );
        
        %analytical curl
        ux = 4 * (2*pi/ l_sc) * cos( 4 *(y.velx-shift) * (2*pi/ l_sc) ) ...
            .* sin( 3* (x.velx-shift) * (2*pi/ l_sc)  ) ;
        uy = -3 * (2*pi/ l_sc) * cos( 3 *(x.vely-shift) * (2*pi/ l_sc) ) ...
            .* sin( 4 * (y.vely-shift) * (2*pi/ l_sc) );
        anal_curl = [ ux; uy ];
    
        %error in gradient
%         ( C*s - anal_curl)'
        errinfC(j) = max(abs( C*s - anal_curl) );
        
    %error for R
        
        %Scale by grid spacing
        
        M_vort = get_hscale_vort( parms );
%         M_vort_sqrt = sqrt( M_vort );
        M_veli = M_vel;
        M_veli( M_veli ~= 0 ) = 1./M_veli( M_veli ~= 0 );
        
        R = M_vort * R * M_veli ;
        
        %function to take div of
        %Get interval for making sins periodic:
        l_sc = parms.mg * parms.len;
        shift = parms.xmin;

        uvect = [ sin( (x.velx-shift) * (2*pi/ l_sc) ) .* ...
            sin( (y.velx-shift) * (2*pi/ l_sc) ) ; ...
            sin( 2*(x.vely-shift) * (2*pi/ l_sc) ) .* ...
            sin( 3*(y.vely-shift) * (2*pi/ l_sc) ) ];


        %analytical rot
        anal_rot = 2*(2*pi/ l_sc) * cos( 2*(x.gamm-shift) * (2*pi/ l_sc)) ...
            .*  sin( 3*(y.gamm-shift) * (2*pi/ l_sc) ) - ...
            (2*pi/ l_sc) * cos( (y.gamm-shift) * (2*pi/ l_sc)) ...
            .*  sin( (x.gamm-shift) * (2*pi/ l_sc) );

        %error in gradient
%         ( R*uvect - anal_rot)'
        errinfR(j) = max(abs( R*uvect - anal_rot) );
        
    
end

loglog( 1./mvect, errinfC, 'bx' ), hold on
loglog( 1./mvect, errinfR, 'ro' ), hold on
% loglog( 1./mvect, 30 * (1./mvect).^2, '--')
loglog( 1./mvect, 6 * (1./mvect), '--')



