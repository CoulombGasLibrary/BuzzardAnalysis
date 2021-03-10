% Finds the variance of points with a certain distance. That is, the first
% entry in CovVec is the covariance between points in DataVec of distance
% 1, the second entry is of distance 2, and so on. The zeroth entry would
% be std(DataVec), though that is not calculated in this function.
% Range gives the length of CovVec.
% Periodic determines whether DataVec should be considered periodic or not.

function CovVec = CovCons(DataVec,Range,Periodic)

    if nargin < 3
        Periodic = false;
    end

    CovVec = zeros(Range,1);
    N = length(DataVec);

    for iRang = 1:Range
        
        if Periodic == true
            GatherVec = ones(length(DataVec),2)*NaN;
        else
            GatherVec = ones(length(DataVec)-Range,2)*NaN;
        end
        
        for iPick = 1:length(GatherVec(:,1))
            GatherVec(iPick,:) = [DataVec(iPick),DataVec(mod(iPick + iRang - 1,N)+1)];
        end
        
        CovVec(iRang) = mean(prod(GatherVec,2)) - prod(mean(GatherVec));
        
    end
end