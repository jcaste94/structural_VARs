function [ varirf_id_l, varirf_id_h, q_id, q_indx, GTq] = construct_idset( q_grid, ...
          signrestrictionsindx, signrestrictions, nirf, YYirf)
% This procedure generates the identified set
% Step 1: id set for q
% Step 2: compute th for each q in id set 

% --------------------------------
% STEP 1 : construct id set for q
% --------------------------------
q_dim          = size(q_grid,2);
nvar           = size(q_grid,1);
mu_dim         = size(signrestrictionsindx,1);

% select the responses that are restricted
[YYirfselect]  = select_irf(YYirf,signrestrictionsindx);

% create column vector, omitting zero responsed due to Cholesky
[psi_restricted, index_phi] = construct_psi_restricted(YYirfselect);

q_id       = []; %initialize id set

q_indx = zeros(1,q_dim);

for i = 1 : q_dim
    
    q   = q_grid(:,i); 
    Sq  = kron(eye(mu_dim),q');
    Sq  = Sq(:,index_phi);

    % compute objective function
    csi_hat  = (Sq*psi_restricted).*signrestrictionsindx(:,3);
    mu_hat   = csi_hat .*(csi_hat >= 0 );
    delta    = csi_hat - mu_hat;
    GT_val   = delta'*delta;
    
    if GT_val <= 1E-10
        % include q into IDset
        q_id= [q_id , q];
        q_indx(1,i) = 1;
    end
    
    GTq(1,i) = GT_val;
    
    clear q Sq 
    
end %end of loop for grid of q

% --------------------------------------------------------
% STEP 2 : construct identified set for impulse responses
% --------------------------------------------------------
nq          = size(q_id,2);
varirf_id_l = zeros(nirf,nvar);
varirf_id_h = zeros(nirf,nvar);

% Iterate over variables/horizons 
for varindx = 1: nvar
    
    horizonindx = 1;
    for horizonindx = 1:nirf
        
        varirf_temp = zeros(nq,1);
        
        for q_ind = 1: nq
            q = q_id(:,q_ind);
            varirf_temp(q_ind,1) = YYirf(horizonindx,(varindx-1)*nvar+1:varindx*nvar)*q; %phi'*q
        end %end of loop for q_ind
               
        varirf_id_l(horizonindx,varindx) = min(varirf_temp);
        varirf_id_h(horizonindx,varindx) = max(varirf_temp);
        
        % truncate id set in view of sign restr
       
            if horizonindx <= size(signrestrictions,1) 
             if signrestrictions(horizonindx ,varindx)  == -1
                %  response is restricted to be <= 0
                varirf_id_l(horizonindx,varindx) = min( 0 , varirf_id_l(horizonindx,varindx)  );
                varirf_id_h(horizonindx,varindx) = min( 0 , varirf_id_h(horizonindx,varindx)  );
             elseif signrestrictions(horizonindx ,varindx) == 1
                %  response is restricted to be >= 0
                varirf_id_l(horizonindx,varindx)  = max( 0 , varirf_id_l(horizonindx,varindx)  );
                varirf_id_h(horizonindx,varindx)  = max( 0 , varirf_id_h(horizonindx,varindx)  );
            end
         end
       
    end % end of loop for horizonindx
    
end % end of loop for variable index

