%==========================================================================
%           ECONOMETRICS IV - SPRING 2020 - Prof. Schorfehide 
%                        STRUCTURAL VARs                     
%
% Author: Juan Castellanos Silván using replication files of Granziera, 
% Moon and Schorfheide (2017)
%
% Date: 24/04/2020
%==========================================================================

% PROBLEM 1. Impulse Responses, Bayesian CS, Identified Sets. 


%% HOUSEKEEPING

tic 
clear all
close all
clc


%% REDUCED-FORM VAR RESULTS

file = "/Users/Castesil/Documents/GitHub/Econ 722 - Schorfheide/PS2/reduced_form_VARs/resultsVAR.mat";
load(file,'Phip', 'Sigmap', 'X', 'Y','Phi_tilde');

T = size(X,1);                          % number of observations
ndraws = size(Phip,1);                  % number of simulated draws
nvar = size(Phip,3);                    % number of variables
nlags = floor(size(Phip,2)-1)/nvar;     % number of lags 

if size(Phip,2) == nvar*nlags
    itercept = 0;
else
    intercept = 1;
end

%% SIGN RESTRICTIONS

nirf = 40;     % maximum horizon of IRFs
nirfsign = 2;  % horizon for IRF restrictions ==> restrict IRFs on impact and one quarter after the shock

% restricted responses: inflation <= 0; interest rates >=0; money <=0;
signrestrictions      = zeros(nirfsign,nvar);
signrestrictions(:,2) = -ones(nirfsign,1);
signrestrictions(:,3) =  ones(nirfsign,1);
signrestrictions(:,4) = -ones(nirfsign,1);

% Construct summary of sign restrictions
NumberSignRestr      = sum(sum(abs(signrestrictions)));
signrestrictionsindx = zeros(NumberSignRestr,3);

count = 0;
for i = 1: nirfsign
    for j =1: nvar
        if signrestrictions(i,j) ~= 0
            
            count=count+1;
            
            % first column says which variable is restricted for each lag
            signrestrictionsindx(count,1) = j;
            
            % second column tells you which lag
            signrestrictionsindx(count,2) = i;
            
            % third column tells you which restriction is imposed
            signrestrictionsindx(count,3) = signrestrictions(i,j);
        
        end
    end
end


%% IMPULSE RESPONSES AND BAYESIAN CREDIBLE SETS 

disp('************************************');
disp('Computing IFR & Bayesian Cred. Bands');
disp('                                    ');

%=========================================================================
%             IMPULSE RESPONSES W/ ACCEPTANCE SAMPLER
%=========================================================================

nvalid = 0; % number of draws that satisfy the resrictions
ndrawsout = 2*1E3;
iDraws = 1; 

while nvalid < ndraws
    
    % recover results or draw new values from posterior distribution
    if iDraws < ndraws
        
        sigma = squeeze(Sigmap(iDraws,:,:));
        phi = squeeze(Phip(iDraws,:,:));
        
        iDraws = iDraws + 1;
    else
        
        % OLS estimator 
        Phi_ols = inv(X'*X)*X'*Y;
        SSR = (Y-X*Phi_tilde)'*(Y-X*Phi_tilde);

        % Draws from the density Sigma | Y
        sigma   = iwishrnd(SSR,T-nvar*nlags-1);
       
        % Draws from the density vec(Phi) |Sigma(j), Y
        phi_new = mvnrnd(reshape(Phi_tilde,nvar*(nvar*nlags+intercept),1),kron(sigma,inv(X'*X)));
        
        % Rearrange vec(Phi) into Phi
        phi = reshape(phi_new,nvar*nlags+intercept,nvar);
        
    end
    
    [ YYirf_new, valid ] = construct_rfirfsign(phi, sigma, nirf, signrestrictionsindx);

    if valid == 1
        nvalid = nvalid+1;
        YYirf_valid(:,nvalid) =  YYirf_new;
        %stack the responses: YYirf_valid is 
        %[var1h1;...;var1h40;var2h40;...var2h40;....var4h40]
        
       if mod(nvalid,ndrawsout) == 0
           disp(['Number of Draws Processed ', num2str(nvalid), ' out of ', num2str(ndraws)]);
       end
        
    end
        
end


%=========================================================================
%                    MEANS AND CREDIBLE SETS
%=========================================================================

alpha = 0.10;                   % desired significant level for theta set
bcslow = zeros(nirf, nvar);
bcshigh = zeros(nirf, nvar);
bcsmean = zeros(nirf,nvar);

for i = 1:nvar
    
    bcsmean(:,i) = mean(YYirf_valid((i-1)*nirf+1:i*nirf,:),2);
    
    bcslow(:,i) = prctile(YYirf_valid((i-1)*nirf+1:i*nirf,:)',(alpha/2)*100);
    
    bcshigh(:,i) = prctile(YYirf_valid((i-1)*nirf+1:i*nirf,:)',(1-(alpha/2))*100);
    
end

%=========================================================================
%                            FIGURES
%=========================================================================

figure_list = {'Output', 'Inflation', 'InterestRate', 'RealMoney'};

for i=1:nvar
    
    figure('Name',figure_list{i});
    plot(0:nirf-1, bcsmean(:,i), 'color','k', 'LineWidth', 2);
    hold on
    plot(0:nirf-1, bcshigh(:,i), 'LineStyle', '-.', 'color','k');
    hold on
    plot(0:nirf-1, bcslow(:,i), 'LineStyle', '-.', 'color','k');
    hold on
    yline(0, 'LineStyle','--', 'color','k');
    hold off
    grid on
    xlabel("Time horizon", 'FontSize', 14);
    
    x = 29.7;                  % A4 paper size
    y = 21.0;                  % A4 paper size
    xMargin = 1;               % left/right margins from page borders
    yMargin = 1;               % bottom/top margins from page borders
    xSize = x - 2*xMargin;     % figure size on paper (widht & hieght)
    ySize = y - 2*yMargin;     % figure size on paper (widht & hieght)

    set(gcf, 'Units','centimeters', 'Position',[0 0 xSize ySize]/2)

    set(gcf, 'PaperUnits','centimeters')
    set(gcf, 'PaperSize',[x y])
    set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
    set(gcf, 'PaperOrientation','portrait')
    
    cd('/Users/Castesil/Documents/EUI/Year II - PENN/Spring 2020/Econometrics IV/PS/PS3/LaTeX/')
    saveas(gcf, strcat('pIRF_',figure_list{i},'.pdf'));

end

cd('/Users/Castesil/Documents/GitHub/Econ 722 - Schorfheide/PS3/structural_VARs')

return

%% IDENTIFIED SETS 

%=========================================================================
%                          GENERATE GRID FOR q
%=========================================================================

rng(1) % replicability of results
nq = 2*1E4;
Z = randn(nvar,nq);
for n=1:nvar
    q_grid(n,:) = Z(n,:) ./ norm(Z(n,:));
end


%=========================================================================
%                   COMPUTE IDENTIFIED SETS FOR q AND theta
%=========================================================================

disp('************************');
disp('Computing Identified Set');
disp('                        ');

[ YYirf_hat ] = construct_rfirf(Phi_ols, SSR, nirf);

[ varirf_id_l, varirf_id_u, q_id, q_id_index , GTqpop] = construct_idset(q_grid, ...
    signrestrictionsindx, signrestrictions, nirf, YYirf_hat);

length_id = varirf_id_u - varirf_id_l;



toc