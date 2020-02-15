%==========================================================================
%           ECONOMETRICS IV - SPRING 2020 - Prof. Schorfehide 
%                        STRUCTURAL VARs                     
%
% Author: Juan Castellanos Silván 
% Date  : 24/04/2020
%==========================================================================
%
%=========================================================================
%                             INPUTS
%=========================================================================
%   - Phi::Array =  (np+1) x n matirx of coefficients 
%   - Sigma::Array = n x n variance-covariance matrix 
%   - h::Integer = horizon of the IRFs
%
%=========================================================================
%                             OUTPUTS
%=========================================================================
%
%
function []=IRF(Phi, Sigma, h)
    
    % ------------------------
    % Orthogonal matrix OMEGA
    % ------------------------
    Omega 
    
    % --------------------------
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
    
    
    % 2.2. Wold + Choleski decomposition
    Sigma_tr = chol(Sigma,'lower');
     
    C = zeros(n,n,h);
    for ih = 1:h
        BigC = BigA^(ih-1);
        C(:,:,ih) = BigC(1:n,1:n);
        D(:,:,ih) = C(:,:,ih) * Sigma_tr * Omega;
    end
    
    irf =(reshape(permute(D,[3 2 1]),h,n*n,[]));
    
end


