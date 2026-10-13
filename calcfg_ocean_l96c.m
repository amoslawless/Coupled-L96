function [f,g] = calcfg_ocean_l96c(dX0_o,z_b,innov,z_lin,H_ocean,invB,invR,nsteps,h,na,no,Fx,Fy,alph,gamma,ob_ix,l_fgat_o)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% compute the cost function for the atmospheric part
% x=z_lin(1:na,:);
% y=z_lin(na+1:end,:);

%% Calculate cost function
% Background term
dXb_o = z_b(na+1:end,1) - z_lin(na+1:end,1);
fa = 0.5 * (dX0_o-dXb_o)' * invB * (dX0_o-dXb_o);

z = z_lin;
dz0 = [zeros(na,1);dX0_o];
if l_fgat_o == 1
%    dz(:,1:nsteps+1) = dz0;
    dz = repmat(dz0,1,nsteps+1);
else
    [dz,dfdz,dfdYi] = l96c_rk2_tl(z,dz0,h,nsteps,na,no,Fx,Fy,alph,gamma);
end
for i=1:nsteps
    fa = fa + 0.5 * ob_ix(i)*(innov(:,i)-H_ocean*dz(na+1:end,i+1))'...
        *invR*(innov(:,i)-H_ocean*dz(na+1:end,i+1));
end
f = fa;

% initialising the observation gradient, compute observation gradient at the last step:
g_hat = zeros(na+no,nsteps+1);
g_hat(na+1:na+no,nsteps+1) = g_hat(na+1:na+no,nsteps+1) - ob_ix(nsteps)*H_ocean'*invR*...
    (innov(:,nsteps)-H_ocean*dz(na+1:na+no,nsteps+1));
% adjoint interation:
for i=nsteps:-1:2
%   One step of adjoint model  
    if l_fgat_o == 1
        g_hat(na+1:end,i) = g_hat(na+1:end,i+1);
    else
        dfdz_adj = dfdz(:,:,i);
        dfdYi_adj = dfdYi(:,:,i);
        [g_hat_model] = l96c_rk2_adj(dfdz_adj,dfdYi_adj,g_hat(:,i+1),h,1,na,no);    
        g_hat(na+1:end,i) = g_hat_model(na+1:end,1) - ob_ix(i-1)*H_ocean'*invR*...
            (innov(:,i-1)-H_ocean*dz(na+1:end,i));
    end
end
if l_fgat_o == 1
    g_hat_model(:,1) = g_hat(:,2);
else
    [g_hat_model] = l96c_rk2_adj(dfdz(:,:,1),dfdYi(:,:,1),g_hat(:,2),h,1,na,no);
end
g0 = g_hat_model(:,1);
% l96c_rk2_adjoint(Jacob_out1, Jacob_out2, Yi, z_final, dz_final, h, nsteps)
gJb = invB * (dX0_o-dXb_o);
g = gJb + g0(na+1:end);

end

