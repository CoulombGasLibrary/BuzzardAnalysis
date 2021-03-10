% This file generates a number of realisations of a 2D Coulomb gas. The
% input parameters are chosen in the beginning of the file. The output is
% the real and imaginary parts saved in two separate files named after the
% input parameters and saved in the relative folder /data/.
% These files contain comma and line separated matrices of dimensions
% #particles X #realisations.

clear;

% The number of particles.
N = 200;

% The number of configurations. 1e4 gives reasonably smooth curves, but
% most operations with it requires some time.
Nconf = 1e4;

% The number of iterations. Comparison of \beta = 2 to Ginibre suggests
% that 100 iterations is enough.
Nsteps = 1e2;

% Step size. In accordance with Chafaï and G. Ferré, ep = N^(-1/3) is
% chosen.
ep = N^(-1/3);

% Energy function for a single point. Takes the configuration X and a
% point x with the number id. Only meant for finding the energy difference when
% perturbing a single point. Hence all interactions between points not
% X(id) are ignored.
H = @(x , id , X, N , bet) N * x' * x ... % First term comes from Gaussian
    - sum( log( abs( X( [1:(id-1) , (id+1):N] ) - x )).* bet ); % Second term is logarithmic interaction

% To keep track of the run, print the starting time.
disp('Started:')
disp(datetime('now'))

% The repulsion between points is iterated over, so several can be made in
% the same run.
for betaNum = [0,0.1,0.2,0.3] % Alternatively, the form 0:0.1:2 may be used.
    
    % Show when the program started calculating this value of \beta.
    disp(['Calculating \beta=',num2str(betaNum)])
    disp(datetime('now'))

    % The name of the files that will be saved.
    TEXT = ['Coulomb_N',num2str(N),'_Nsteps',num2str(Nsteps),...
                '_Nconf',num2str(Nconf),'_beta',num2str(betaNum),'_ep',num2str(ep)];

    % The particles will collapse into a point if \beta=0, so this Poisson
    % case has to be treated separately. This is done by generating
    % radially symmetric points and reshaping the radial part.
    if betaNum == 0
        eigsavetemp = randn(N,Nconf) + 1i*randn(N,Nconf);
        eigsave = eigsavetemp ./ abs(eigsavetemp) .* sqrt(rand(N,Nconf));
    else
        
    % Prelocate for the saved eigenvalues.
    eigsave = zeros(N,Nconf);

        parfor itConf = 1:Nconf

            % Start initial configuration as Poisson.
            Config = sqrt(rand(N,1)).*exp(2i*pi*rand(N,1));

            % Perform the given amount of iterations
            for it = 1:Nsteps
                
                % Each iteration includes a loop over the amount of points.
                % Rather than going through each point, we pick a random
                % point N times. On average this is the same.
                for iPoint = 1:N
                    
                    % Pick the random point and name it "id"
                    id = randi(N,1);
                    
                    % Make a perturbed configuration by adding a Gaussian
                    % to the single point.
                    xtilde = Config(id) + ep * randn + ep * 1i*randn;
                    
                    % Calculate the probability of accepting the new
                    % configuration by comparing the energies. If the
                    % energy of the new configuration is lower than the
                    % old, the exponential will be larger than 1, and 1 is
                    % simply picked. The anonymous function H is called
                    % from the beginning of the script.
                    p = min(1 , exp(- H(xtilde , id , Config, N, betaNum)...
                                    + H(Config(id) , id , Config, N, betaNum) ) );

                    % Generate random point and compare to probability.
                    if p > rand
                        % If accepted, update the configuration.
                        Config(id) = xtilde;
                    end

                end

            end
            
            % To keep track of the simulation, print status regularly.
            % This does not work so well inside the parrallel for-loop,
            % because the iterations are not done consecutively, but it
            % still gives a sense of how far the program is.
            if mod(itConf,Nconf/100) == 0
                disp([num2str(itConf),' configurations calculated at'])
                disp(datetime('now'))
            end
            
            % Add the finished configuration
            eigsave(:,itConf) = Config;

        end
    
    end
    
    % Save the real and imaginary parts of the eigenvalues. They are save
    % separately, because numbers are easier to handle than strings, which
    % would be the format, complex numbers are saved in.
    csvwrite(strcat('data/',TEXT,'_eigR.txt'),real(eigsave));
    csvwrite(strcat('data/',TEXT,'_eigI.txt'),imag(eigsave));
    
    % The saving part may be moved inside the previous loop if one wishes
    % to save more often. (For instance if the computer has a tendency to
    % crash or be restarted.) This requires the "parfor" to be a "for"
    % because file access does not work well with parallel loops. (So there
    % will be some time increase.)
    

end