function [ y ] = L1_Spirit_PD3O(g,Mask,Lev,lambda,Itr_Max,Ker,Ker_Tra,Gamma,Delta)
% initialization
Mask_Un = ~Mask;
z = cat(4,real(g),imag(g));
y = z;
s = TightFrameDecomposition(ifft2_3D(Mask_Un .* y), Lev);
c = TightFrameDecomposition(ifft2_3D(z), Lev);
Param = zeros(size(g,1),size(g,2),size(g,3),Lev);
for Itr = 1:Itr_Max
    Qy_z = Mask_Un.*y+z;
    C_I_Qy_z = imfilter(Qy_z,Ker)-Qy_z; %B_y=(C-I)(Qy+z)
    gran_f = Mask_Un.*(imfilter(C_I_Qy_z,Ker_Tra)-C_I_Qy_z);%Qo_Bt_B_g=Q(C-I)^T(C-I)(Qy+z)
    y_Gamma_gran_f = y - Gamma .* gran_f;
    Gamma_Bt_s_y_Gamma_gran_f = Gamma.*Mask_Un.*fft2_3D(TightFrameReconstruction(s, Lev)) - y_Gamma_gran_f;
    t = s - Delta .* TightFrameDecomposition(ifft2_3D(Mask_Un .* Gamma_Bt_s_y_Gamma_gran_f), Lev);
    t_Delta_c = t + Delta .* c;
    if(Itr<10&&mod(Itr,2)==1)
        Param = UpdateRegularizationParameter(Param,lambda);
    end
    s = t_Delta_c - ell2NormProximalOperator(t_Delta_c,Param);
    y = y_Gamma_gran_f - Gamma .* Mask_Un .* fft2_3D(TightFrameReconstruction(s, Lev));
end
y = sos(ifft2_3D(Mask_Un.*y+z));
end