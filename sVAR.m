%==========================================================================
%           ECONOMETRICS IV - SPRING 2020 - Prof. Schorfehide 
%                        STRUCTURAL VARs                     
%
% Author: Juan Castellanos Silván 
% Date  : 24/04/2020
%==========================================================================

% Exercises 1,2,3, & 4: Read VAR.m, calculate IRFs, means and CS

%==========================================================================
%                             HOUSEKEEPING
%==========================================================================

tic 
clear all
close all
clc

%==========================================================================
%                       REDUCED-FORM VAR RESULTS
%==========================================================================

file = "/Users/Castesil/Documents/GitHub/Econ 722 - Schorfheide/PS2/reduced_form_VARs/resultsVAR.mat";
load(file,'Phip', 'Sigmap');

nsim = size(Phip,1);        % number of simulations
n = size(Phip,3);           % number of variables
p = (size(Phip,2)-1)/n;     % number of lags 


%==========================================================================
%                       IFRs + SIGN RESTRICTIONS
%==========================================================================

h = 40;  % time horizon
irf = zeros(h,n,nsim);

for s = 1:nsim
    
    irf(:,:,s) = IRF(squeeze(Phip(s,:,:)), squeeze(Sigmap(s,:,:)), h);

end


%=========================================================================
%                        MEANS AND CREDIBLE SETS
%=========================================================================

% -----
% Means
% ------
irf_mean = mean(irf,3);

% -------------
% Credible set
% -------------
irf_lb = zeros(h,n);
irf_ub = zeros(h,n);

for ih=1:h
    for i=1:n
        
        % lower bound
        irf_lb(ih,i) = prctile(irf(ih,i,:),5);
        
        % upper bound
        irf_ub(ih,i) = prctile(irf(ih,i,:),95);

    end
end


%=========================================================================
%                            FIGURES
%=========================================================================

figure_list = {'Output', 'Inflation', 'InterestRate', 'RealMoney'};

for i=1:n
    
    figure('Name',figure_list{i})
    plot(0:h-1, irf_mean(:,i), 'color','k', 'LineWidth', 2)
    hold on
    plot(0:h-1,irf_lb(:,i), 'LineStyle', '-.', 'color','k')
    hold on
    plot(0:h-1,irf_ub(:,i), 'LineStyle', '-.', 'color','k')
    hold on
    yline(0, 'LineStyle','--', 'color','k')
    hold off
    grid on
    xlabel("Time horizon", 'FontSize', 14)
    
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
    
   saveas(gcf, strcat('pIRF_',figure_list{i},'.pdf'))
    
end

toc