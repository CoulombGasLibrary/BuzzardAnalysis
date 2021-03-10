% Analytic expression for the NNN-spacing of Ginibre variables in the bulk.

function pNNN = GinibreNNN(s,cutoff)

    pNNN = 0;
    
    for j = 1:cutoff
        for k = [1:(j-1),(j+1):cutoff]
            pNNN = pNNN + 2*s.^(2*k+1).*exp(-s.^2)./( gammainc( s.^2,1 + k,'upper') .* gamma(1+k) )...
                .* gammainc( s.^2,1 + j,'lower') ./gammainc( s.^2,1 + j,'upper') ;
        end
    end
    
    pNNN = pNNN * prod(gammainc( s.^2,1 + (1:cutoff),'upper'));
    
end