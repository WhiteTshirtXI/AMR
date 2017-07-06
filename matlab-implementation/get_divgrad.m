function [G, D] = get_divgrad( parms )

%Return the gradient and divergence operators

%--parameters used:
    m = parms.m; n = parms.n; mg = parms.mg;

    %size of mats
    parms.nrowsG = get_velx_ind( m-1, n, mg, parms ) + ...
        get_vely_ind( m, n-1, mg, parms ); 
    parms.ncolsG = get_press_ind( m, n, mg, parms );
    parms.nrowsD = parms.ncolsG; 
    parms.ncolsD = parms.nrowsG;

    G = sparse( parms.nrowsG, parms.ncolsG );
    
%--

%--build main part of matrices 
%  (excluding BCs from inner / outer grid):
    for j = 1 : mg

        G = get_grad( G, j, parms );

    end

    %Divergence is negative of gradient transpose 
    %(except for bcs that will be incorporated below)
    D = -G';
    
%--


%-- Incorporate BCs
    if mg > 1

        for j = 1 : mg

            %divergence requires info from coarser grid
            if j < mg
                D = get_div_BCs( D, j, parms );
            end

            %gradient requires info from finer grid
            if j > 1

                G = get_grad_BCs( G, j, parms );

            end


        end
    end
%--



