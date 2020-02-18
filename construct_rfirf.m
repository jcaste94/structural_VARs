function [ YYirf ] = construct_rfirf(PHI, SIGMA, nirf);

nvar        = size(PHI,2);
nlags       = floor(size(PHI,1)/nvar);
sigmatr_mle = chol(SIGMA,'lower'); 
% remove the constant term (last row)
phi_mle     = PHI(1:nvar*nlags,:); 

% for lag = 1:nlags
%     phi_irf(nvar*(lag-1)+1:nvar*lag,1:nvar) = phi_mle(nvar*(lag-1)+1:nvar*lag,:)';   %phi_mle is [A1;A2;A3;..]; phi_irf is [A1';A2';A3';..] 
% end
% Compute Impulse Responses to Orthogonalized shock
% Format:
% h=1: (sh1 -> v1) (sh1 -> v2) (sh2 -> v1) (sh2 -> v2) 
% h=2: (sh1 -> v1) (sh1 -> v2) (sh2 -> v1) (sh2 -> v2) 

YYirf  = zeros(nirf,nvar^2);

for i =1: nvar
    %response of all variables to the i-th orthogonalized shock 
    imp    = sigmatr_mle(:,i)'; 
    XXpred = zeros(1,nvar*nlags);
    %loop over horizon
    for j=1:nirf 
          if j == 1
             YYirf(j,1+(i-1)*nvar:i*nvar) = XXpred*phi_mle + imp;
          else
             YYirf(j,1+(i-1)*nvar:i*nvar) = XXpred*phi_mle;
          end
          if nlags > 1
             XXpred(nvar+1:nvar*nlags) = XXpred(1:(nlags-1)*nvar);
          end
          XXpred(1:nvar) = YYirf(j,1+(i-1)*nvar:i*nvar);
    end
end

MM             = zeros(nvar^2,nvar^2);
IdentityMatrix = eye(nvar);

for i =1: nvar
   MM(1+(i-1)*nvar:i*nvar,:) = kron(IdentityMatrix,IdentityMatrix(i,:));
   i = i+1;
end

% Reorder Impulse Responses
% h=1: (sh1 -> v1) (sh2 -> v1) (sh1 -> v2) (sh2 -> v2) 
% h=2: (sh1 -> v1) (sh2 -> v1) (sh1 -> v2) (sh2 -> v2) 

YYirf = YYirf*MM;
end

