function press_ind = get_press_ind( i, j, glev, grid_parms )

m = grid_parms.m;
n = grid_parms.n;

%Return vector indices for pressure (press_ind)
%corresponding to 2d array indices.

% e.g. get_press_ind(2, 3, 1, grid_parms) returns 2*(m) + 2, which is the
% vector index of the pressure at the point with xindex = 2, yindex = 3,
% and gridlevel = 1.

%Note: pressure on grid 1 is of size m * n.
%   For all other grids it is of size m*n/2 + m/2*n/2.


%1st gridlevel is easy...
if glev == 1
    
    press_ind = m .* (j-1) + i ;
    
else
    
    n_add = (m*n) + (glev-2) * (m*n/2 + m/2 * n/2);
    
   %points that are below finer grid:
    press_ind = (j <= n/4) .* ( n_add + m .* (j-1) + i );    
    
    %points that contain finer grid and are to left of finer grid:
    n_bott = m*n/4;
    press_ind = press_ind + (j > n/4 & j <= 3*n/4 & i <= m/4 ) .* ...
        ( n_add + n_bott + m/2*(j-n/4-1) + i );
    
    %points that contain finer grid and are to right of finer grid:
    press_ind = press_ind + (j > n/4 & j <= 3*n/4 & i > 3*m/4) .* ...
        ( n_add + n_bott + m/2*(j-n/4-1) + (i-m/2) );
    
    %points that are above finer grid:
    n_mid = m/2* n/2 ;
    press_ind = press_ind + (j > 3*n/4 ).* ...
        ( n_add + n_bott + n_mid + (j-3*n/4-1)*m + i );
    
end