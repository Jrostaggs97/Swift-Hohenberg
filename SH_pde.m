function [phi_hat_a] = SH_pde(t,phi_hat_vec,alpha,beta,k_lap,k_biharm,s)

    %reshape to matrix (really only need for ifft to calculate NL part)
    phi_hat = reshape(phi_hat_vec,s,s);
    
    %let's go ahead and get phi in the time domain for computation of the
    %nonlinear part
    phi_time = real(ifft2(phi_hat)); %could try not taking real, but I guess real is what will be plotted. idk yet.
   
    %Linear part au - (1 + Laplacian)^2 phi 
    %(calculated in the frequency domain) 
    Lin_hat =  -1*(2*k_lap.*phi_hat_vec + k_biharm.*phi_hat_vec);
    
    %Non-linear (calculated in the time domain)
    N_Lin = (alpha - 1)*phi_time +  beta*phi_time.^2 - phi_time.^3;
    
    %Take our non-linear piece into the frequency domain
    N_Lin_hat = reshape(fft2(N_Lin),s^2,1);
    
    %The wonderous Swift Hohenberg Equation
    phi_hat_a = Lin_hat + N_Lin_hat;

end