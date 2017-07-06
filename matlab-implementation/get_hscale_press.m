function M_press = get_hscale_press( parms )

%Get diagonal matrix that scales the gradient by the grid spacing h

m = parms.m; n = parms.n; mg = parms.mg; len = parms.len;

del = len / m;

%Get size of matrix
nrows = get_press_ind( m, n, mg, parms );
ncols = nrows ;

M_press = sparse( nrows, ncols );

for glev = 1 : mg
    
    %grid spacing on current grid
    delb = del * 2.d0^(glev-1);
    
    if glev == 1
        ind_s = 1;
    else
        ind_s = 1 + get_press_ind(m,n,glev-1,parms);
    end
    
    ind_e = get_press_ind(m,n,glev,parms);
    
    ind = ind_s : ind_e;
    
    M_press = M_press + sparse( ind, ind, ...
        (1/delb) * ones(size(ind)), nrows, ncols);
    
end



