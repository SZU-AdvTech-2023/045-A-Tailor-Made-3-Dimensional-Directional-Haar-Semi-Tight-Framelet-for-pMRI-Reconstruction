function [ y ] = L1_Spirit_ADMM(g,Mask,Lev,lambda,Rho,Iter_Out,Iter_In,Ker, Ker_Tra)
% initialization
Mask_Un = ~Mask;
W_IFFT_g = TightFrameDecomposition(ifft2_3D(g), Lev);
v = zeros(size(W_IFFT_g));
z = v;
A_g = imfilter(g,Ker)-g; %A_g=(C-I)g
Qo_At_A_g = Mask_Un.*(imfilter(A_g,Ker_Tra)-A_g);%Qo_At_A_g=Q(C-I)^T(C-I)g
u = g; 
Param = zeros(size(g,1),size(g,2),size(g,3),Lev);
for k = 1:Iter_Out
    % u-subquestion
    FFT_Wt_z_rho_Wt_y = fft2_3D(TightFrameReconstruction(Rho*v+z, Lev));
    target_b = FFT_Wt_z_rho_Wt_y.*Mask_Un - Qo_At_A_g - Rho*(Mask_Un.*g);
    u = CG(target_b,u,Iter_In,Mask_Un,Rho,Ker,Ker_Tra);
    % v-subquestion
    if(Itr<10&&mod(Itr,2)==1)
        Param = UpdateRegularizationParameter(Param,lambda);
    end
    t = TightFrameDecomposition(ifft2_3D(Mask_Un.*u+g), Lev) - z./Rho;
    v = ell1NormProximalOperator(t,Param);  
    % alpha
    z = z + Rho.*(v-TightFrameDecomposition(ifft2_3D(Mask_Un.*u+g), Lev));
end
y = sos(ifft2_3D(Mask_Un.*u+g));
end

function [res] = Afun(x,Mask_Un,Rho,Ker,Ker_Tra)
Qo_u = x.*Mask_Un;
A_Qo_u = imfilter(Qo_u,Ker)-Qy_z;
Qo_At_A_Qo_u = Mask_Un.*(imfilter(A_Qo_u,Ker_Tra)-A_Qo_u);
res = Qo_At_A_Qo_u + Rho.*Qo_u;
end

function [ u ] = CG(target_b, u, Iter, Mask_Un, Rho, Ker, Ker_Tra)
r = target_b - Afun(u, Mask_Un, Rho, Ker, Ker_Tra);
d = r;
r_t_r = r(:)'*r(:);
for k = 1:Iter
    Ad = Afun(d, Mask_Un, Rho, Ker, Ker_Tra);
    a = r_t_r/(d(:)'*Ad(:));
    u = u + a*d;
    r = r - a*Ad;
    r1_t_r1 = r(:)'*r(:);
    b = r1_t_r1/r_t_r;
    d = r + b*d;
    r_t_r = r1_t_r1;
end
end
