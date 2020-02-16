% =========================================================================
%           ECONOMETRICS IV - SPRING 2020 - Prof. Schorfehide 
%                        STRUCTURAL VARs                     
%
% Author: Juan Castellanos Silván 
% Date  : 24/04/2020
% =========================================================================
%
% =========================================================================
%                             INPUTS
% =========================================================================
%   - Phi::Array =  (np+1) x n matirx of coefficients 
%   - Sigma::Array = n x n variance-covariance matrix 
%   - h::Integer = horizon of the IRFs
%
% =========================================================================
%                             OUTPUTS
% =========================================================================
%   - D_wold::Array = h x n matrix of impulse responses
%
% =========================================================================
function [D_wold]=IRF(Phi, Sigma, h)
        
    % ------------
    % Housekeeping
    % ------------
    
    n = size(Phi,2);                % number of variables
    p = (size(Phi,1)-1)/n;          % number of lags 
    
    Sigma_tr = chol(Sigma,'lower'); % Choleski decomposition
    
    % --------------------------
    % Prior: Acceptance sampling
    % --------------------------
    iter = 0;
    maxIter = 1000;
    while iter<maxIter
        
        % proposal draws
        Z = randn(n,1);   
        Omega1 = Z/norm(Z);

        % impose restrictions
        R = Sigma_tr * Omega1; 
        if R(2) < 0 && R(3) > 0 && R(4) < 0
            break
        end
        
        % update
        iter = iter+1;
    end

    
    % ---------------------------
    % Impulse Response Functions
    % ---------------------------
    
    % 2.1. Companion form 
    Phi_aux=zeros(n, n*p);
    j=1;
    for i=1:p
        phi = Phi(j:i*p,:);
        Phi_aux(:,j:i*p)=phi;
        j = i*p+1;
    end
    BigA=[Phi_aux; eye(n*p-n) zeros(n*p-n,n)]; % np x np matrix
    
    
    % 2.2. Wold representation
    C = zeros(n,n,h);
    D = zeros(n,h);
    for ih = 1:h        
        BigC = BigA^(ih-1);
        C(:,:,ih) = BigC(1:n,1:n);
        D(:,ih) = C(:,:,ih) * Sigma_tr * Omega1;
    end

    D_wold = D';
end


