function [YYirfselect]  = select_irfs(YYirf,signrestrictionsindx);
% select the rfIRF that are sign restricted 
nvar        = sqrt(size(YYirf,2));
YYirfselect = zeros(size(signrestrictionsindx,1),nvar);
% signrestrictionsindx = zeros(NumberSignRestr,3); 
% first column: restricted variable
% second column: restricted lag
% third column: sign
for i = 1: size(signrestrictionsindx,1);
   vindx = signrestrictionsindx(i,1);
   hindx = signrestrictionsindx(i,2);
   YYirfselect(i,:) = YYirf(hindx,1+(vindx-1)*nvar:vindx*nvar);
end;


