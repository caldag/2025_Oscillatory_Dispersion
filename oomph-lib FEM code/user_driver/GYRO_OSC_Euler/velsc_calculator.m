% Part of the supplementary material for

% Caldag, H.O. & Bees, M.A. Fine-tuning the dispersion of active suspensions using oscillatory flows.

% velsc_calculator computes the velocity scaling factor that ensures the amplitude of the velocity field
% is equal to the root-mean-square velocity-based Peclet number.

clearvars;
ii=0;    % Loop counter
Sc=16.8; % Schmidt number
Wos=[0.3 0.4 0.5]; % Enter the desired Womersley numbers
for Wo = Wos
	ii=ii+1; % Increase counter
	period=2*pi/(Wo*Wo*Sc); % Period for given Wo
	yveccal=-1:0.01:1;
	tveccal=linspace(0,period,202)'; % Time and lateral position vectors
	fy=1/(1i*Wo^2)*( cosh(sqrt(1i)*Wo.*yveccal)/cosh(sqrt(1i)*Wo) - 1); % y-dependent part of the velocity field
	Woflow=real(fy.*exp(1i.*tveccal.*Wo^2*Sc)); % Complete Womersley flow field
	Woflow_rms=sqrt(1/(4*pi)*trapz(yveccal,trapz(tveccal.*Wo^2*Sc,Woflow.^2))); % Root-mean-square calculation
	velsc(ii)=1/Woflow_rms % Velocity scale is 1 over this root-mean-square value
end

disp(Wos);
disp(velsc);