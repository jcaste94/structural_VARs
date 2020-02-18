function [psi_restricted, index_phi] = construct_psi_restricted( YYirfselect )
% column vector that contains restricted IRFs
% without zero elements that arise from Cholesky Decomp
psi_restricted = reshape(YYirfselect',size(YYirfselect,1)*size(YYirfselect,2),1); 

% construct S_phi
S_phi = eye(size(psi_restricted,1));
index_phi = [];
for i = 1:size(psi_restricted,1)
    if psi_restricted(i,1)~=0
         % this index keeps track of the rows in phi_hat that are zero 
         % (arising because some IRF are zero on impact)
        index_phi = [index_phi,i];
    end
end
S_phi = S_phi(:,index_phi);

% psi_restricted once zero elements have been dropped
psi_restricted = S_phi'*psi_restricted; 

end

