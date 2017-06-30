clear all, close all, clc

%Test the scalar Laplacian obtained from D * G. 

mvect = [20 40 80 160 320];

errinf = zeros( length( mvect), 1);

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
    [G, D] = get_divgrad( parms );

    %Scale by grid spacing
    
%         M_vel = get_hscale_vel( parms );
%         G = M_vel * G;
        
        M_press = get_hscale_press( parms );
        
        D = M_press.^2 * D ;
        L = D * G;
        
%         full(L)
%         
%         pause

    %Get interval for making sins periodic:
        l_sc = parms.mg * parms.len;
        shift = parms.xmin;
        
        p = -cos( 3 * (x.press - shift) * (2*pi/l_sc) )  - ...
            cos( 4 * (y.press - shift) * (2*pi/l_sc) ) ;

    %analytical lap

        anal_lap = 9 * (2*pi/ l_sc)^2 * cos( 3 * (x.press-shift) * ...
            (2*pi/ l_sc)) + 16 * (2*pi/ l_sc)^2 * ...
            cos( 4 * (y.press - shift) * (2*pi/ l_sc));

    %error in lap
%         (L*p - anal_lap)'/ max(abs( anal_lap ) )
    
        errinf(j) = max(abs( L*p - anal_lap) ) / max(abs( anal_lap ) );
        
    
end

loglog( 1./mvect, errinf, 'bx' ), hold on
% loglog( 1./mvect, errinfD, 'ro' ), hold on
% loglog( 1./mvect, 30 * (1./mvect).^2, '--')
loglog( 1./mvect, 6 * (1./mvect), '--')



