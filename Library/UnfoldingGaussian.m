% Calculating the unfolded distances between points based on Gaussian
% distributions at each point used to find the local metric of the points.

function NewDistances = UnfoldingGaussian(sigs,Pos)
    
    N = length(Pos(:,1));
    
    % Prelocate for the density at each point.
    rhoVec = zeros(N,1);
    
    PDF = @(x,y,x0,y0,sig) 1/(2*pi*sig^2*N) * sum( exp( - ( (x-x0).^2 + (y-y0).^2 )/(2*sig^2)) );
    
    % For each point find the density
    parfor it = 1:N    
        rhoVec(it) =  PDF( Pos(it,1) , Pos(it,2) , Pos(:,1) , Pos(:,2) , sigs(it) );
    end
    
    % The unfolded distances are the non-unfolded distances multiplied by
    % the square root of the density at the point from which the distance
    % is calculated.
    NewDistances = sqrt(abs((repmat(Pos(:,1),1,length(Pos(:,1)))-repmat(Pos(:,1)',length(Pos(:,1)),1)).^2 ...
              + (repmat(Pos(:,2),1,length(Pos(:,2)))-repmat(Pos(:,2)',length(Pos(:,2)),1)).^2)) * diag(sqrt(rhoVec));
end