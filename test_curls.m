clear all, close all, clc

%Test the gradient and divergence operators. 

mvect = 20;%[20 40 80 160];

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
    
    errDC =  max(max(abs( D * M_vel * C ))) 
    errRG = max( max( abs( R * G ) ) )
    
    pause

    %error for C
        
        %Scale by grid spacing
        M_vel = get_hscale_vel( parms );
        G = M_vel * G;

        %function to take gradient of
        p = sin( 3* x.press ) .* y.press + ...
            ( x.press ).^2 .* cos( y.press );

        %analytical gradient
        dpdx = 3.*cos(3.*x.velx).*y.velx + 2*x.velx.*cos(y.velx);
        dpdy = sin(3.*x.vely) - ( x.vely ).^2 .* sin( y.vely );
        anal_grad = [ dpdx; dpdy ];
    
        %error in gradient
        errinfG(j) = max(abs( G*p - anal_grad) );
        
%     %error for R
%         
%         %Scale by grid spacing
%         M_press = get_hscale_press( parms );
%         M_veli = M_vel;
%         M_veli( M_veli ~= 0 ) = 1./M_veli( M_veli ~= 0 );
%         D = M_press.^2 * D * M_veli ;
%         
%         %function to take div of
%         %Get interval for making sins periodic:
%         l_sc = parms.mg * parms.len;
%         shift = parms.xmin;
% 
%         uvect = [ sin( (x.velx-shift) * (2*pi/ l_sc) ) .* y.velx.^3; ...
%             sin( (y.vely-shift) * (2*pi/ l_sc) ) .* x.vely.^3];
% 
% 
%         %analytical div
% 
%         anal_div = (2*pi/ l_sc) * cos( (x.press-shift) * (2*pi/ l_sc)) ...
%             .* y.press.^3 + (2*pi/ l_sc) * cos( (y.press - shift) ...
%             * (2*pi/ l_sc)) .* x.press.^3;
% 
%         %error in gradient
%         errinfD(j) = max(abs( D*uvect - anal_div) );
%         
    
end

loglog( 1./mvect, errinfC, 'bx' ), hold on
loglog( 1./mvect, errinfR, 'ro' ), hold on
% loglog( 1./mvect, 30 * (1./mvect).^2, '--')
loglog( 1./mvect, 6 * (1./mvect), '--')



