function M_vel = get_hscale_vel( parms )



m = parms.m; n = parms.n; mg = parms.mg; len = parms.len;

del = len / m;

%Get size of matrix
nrows = get_velx_ind( m-1, n, mg, parms ) + ...
    get_vely_ind( m, n-1, mg, parms ) ;
ncols = nrows ;

M_vel = sparse( nrows, ncols );

for glev = 1 : mg
    
    %grid spacing on current grid
    delb = del * 2.d0^(glev-1);
    
    %-- x-vel block
        if glev == 1
            ind_s = 1;
        else
            ind_s = 1 + get_velx_ind(m-1,n,glev-1,parms);
        end

        ind_e = get_velx_ind(m-1,n,glev,parms);
        
        ind = ind_s : ind_e;
        
        M_vel = M_vel + sparse( ind, ind, (1/delb) ...
            * ones(size(ind)), nrows, ncols);
    %--
    
    %-- y-vel block
        if glev == 1
            ind_s = 1;
        else
            ind_s = 1 + get_vely_ind(m,n-1,glev-1,parms);
        end

        ind_e = get_vely_ind(m,n-1,glev,parms);
        
        ind = ind_s : ind_e;
        
        %y-vel index starts after x-vels
        n_add = get_velx_ind(m-1,n,mg, parms);
        ind = ind + n_add;
        
        M_vel = M_vel + sparse( ind, ind, (1/delb) ...
            * ones(size(ind)), nrows, ncols);
    %--
    
    
    
end







