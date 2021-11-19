function [X,Y,time,phi_plot,run_time] = SH_fft_main(alpha,beta,s,L,delta_t,t_end)
%------ Set up parameters and discretization -------%

%parameters
s=32; %spatial discretization; to be input parameter
L=2.5; %length of UNIFORM domain; 

alpha = .1; %How far the temp is above min temp difference to cause convection (so very small alpha are probably of interest) "temperature knob. negative there is no convection (test and show), fo
beta = 2; %Controls the quadratic non-linearity 
delta_t = .2;
t_end = 20;
t_span = 0:delta_t:t_end; %time span

x_space = linspace(0,2*L,s+1); %Periodic bc so add one
x_span = x_space(1:s); %then chop it off

%uniform discretization
delta = L/s;
y_space = x_space;
y_span = x_span;

%Use meshgrid to put X,Y into correctly oriented grid 
[X,Y] = meshgrid(x_span);

%initial conditions %try sines and cosines

phi_0 = .2*rand(s);

%Perturbed zero condition
% triv_pert_x = zeros(64,64);
% pert_x = .1*rand(64,1);
% triv_pert_x(1,:) = pert_x;

%cool functions and parameters to try out
% cos(X) + sin(X) + sin(Y) + cos(Y);
% try with beta = 1 alpha = .7
% (1/20)*(cos(X) + sin(2*X) + sin(Y) + cos(2*Y)) + exp(-((X-5*pi).^2 + (Y-5*pi).^2)) + exp(-((X-5*pi).^2 + (Y-15*pi).^2))+ exp(-((X-15*pi).^2 + (Y-15*pi).^2))+ exp(-((X-15*pi).^2 + (Y-5*pi).^2))+exp(-((X-10*pi).^2 + (Y-10*pi).^2));
% (1/20)*(cos(X) + sin(2*X) + sin(Y) + cos(2*Y)) + exp(-((X-5*pi).^2 + (Y-5*pi).^2)) + exp(-((X-5*pi).^2 + (Y-15*pi).^2))+ exp(-((X-15*pi).^2 + (Y-15*pi).^2))+ exp(-((X-15*pi).^2 + (Y-5*pi).^2))+exp(-((X-10*pi).^2 + (Y-10*pi).^2));
%(sqrt((sqrt((X.*sin(Y) + X.*cos(Y)).^2 + (Y*cos(Y) - Y.*sin(X)).^2)).^2 + (Y.*sin(Y) + X.*cos(X)).^2) - 2).^2 + (X.*cos(X) - X.*sin(X)).^2;
%cosh(1i.*X)./(Y+.5); try with alpha = .3, beta = 1;
% X.^3 - Y.^2 + cos(X); alpha = ~.7 beta = -.5
% (1/20)*(cos(X) + sin(2*X) + sin(Y) + cos(2*Y)) + exp(-((X-5*pi).^2 + (Y-5*pi).^2)) + exp(-((X-5*pi).^2 + (Y-15*pi).^2))+ exp(-((X-15*pi).^2 + (Y-15*pi).^2))+ exp(-((X-15*pi).^2 + (Y-5*pi).^2))+exp(-((X-10*pi).^2 + (Y-10*pi).^2));
%cos(X) + sin(X) + sin(Y) + cos(Y); try with beta = 1 alpha = .7
%cos(X) -Y ;


colormap("jet");
pcolor(X,Y,real(phi_0));shading interp;colorbar; xlabel("x"); ylabel("y");


% 
%---------- Taking everything into Freq Domain for Spectral integration in
%Matlab ---------%

%take initial phi_0 into to freq domain
phi_0_hat_mat = fft2(phi_0);

%make fourier matrices into vectors
phi_0_hat_vec = reshape(phi_0_hat_mat,s^2,1);

%Set up Ks for fourier transform, where Ks are the wavenumber
k_x_v = (2*pi/L)*[0:(s/2 -1) (-s/2):-1]; %rescale to periodic domain
[k_x,k_y] = meshgrid(k_x_v); %matrix of wavenumbers in x,y
k_lap_mat = -(k_x.^2 +k_y.^2); %square each component (to account for FT of Laplace)
k_biharm_mat = k_lap_mat.^2; %This the biharmonic operator! "Laplacian squared" 
k_lap = reshape(k_lap_mat,s^2,1); %Laplace FT in vector form
k_biharm = reshape(k_biharm_mat,s^2,1); %reshape into vector 


tic;
%integrate in time with ode45 - output still in freq domain

[time,phi_hat_a] = ode45(@(t,phi_hat_vec) SH_pde(t,phi_hat_vec,alpha,beta,k_lap,k_biharm,s),t_span,phi_0_hat_vec);

run_time = toc;
%------ Plotting ----- %
% colormap("jet");colorbar;

phi_plot = real(ifft2(reshape(phi_hat_a(length(time),:),s,s)));

pcolor(X,Y,phi_plot);shading interp;
xlabel("X");
ylabel("Y");
title("beta =0, alpha =",alpha);

%in order to plot we first need to bring our ode45 output into 
% 

%save as gif
% Video = VideoWriter("SH_MOVIE TEST");
% Video.FrameRate=60;
% open(Video)
% 
% colormap("jet"); colorbar;
% hold on;
% 
% %pre allocate movie frame
% for t = 1:length(time)
% 
%     phi_plot = real(ifft2(reshape(phi_hat_a(t,:),s,s)));
%     pcolor(X,Y,phi_plot);shading interp;
% 
%    title("Phi(X,Y), t =",t*delta_t);
%    xlabel("X");
%    ylabel("Y");
% %    drawnow;
%    pause(.25);
%    frame = getframe(gcf);
%    writeVideo(Video,frame);
% end
% hold off;
% close(Video);
% 
% videoinfo = aviinfo('SH_MOVIE TEST.avi');
% gif_name = 'SH_MOVIE TEST.gif'; 
% avi_name = VideoReader('SH_MOVIE TEST.avi');
% Frames = read(avi_name);
% 
% for n = 1:videoinfo.NumFrames
%       [a,b] = rgb2ind(Frames(:,:,:,n),255);
%       if n == 1
%           imwrite(a,b,gif_name,'gif', 'Loopcount',inf);
%       else
%           imwrite(a,b,gif_name,'gif',"DelayTime",.1,'WriteMode','append');
%       end
% end

end