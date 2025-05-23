% Part of the supplementary material for

% Caldag, H.O. & Bees, M.A. Fine-tuning the dispersion of active suspensions using oscillatory flows.

% This code processes the raw simulation data and generates dispersion
% related parameters such as the variance etc.

load('GYRO_OSC_data.mat'); % This is the ouput from GYRO_OSC_sim.m, not supplemented due to large file size.
swcount=5000;
Sc=16.8;
pp=0; % Plot loop counter

%% Dispersion computation and t_mix collection

for Wo=Wos % Going through each Wo

    pp=pp+1;
    yall=datt(pp).yall; % Pick up the yall which contains both x and y data
    SNRt(pp)=datt(pp).SNRt;

    period=2*pi/(Wo*Wo*Sc);
    ncycles=ceil(timeabs/period);
    tfinal=ncycles*period;
    tstep=period/8; % Period computation

    timevec=[1:81 ncycles*8-79:ncycles*8].*tstep; % Time steps recorded

    xdat(:,:,pp)=yall(:,1:swcount); % Receive x position data
	% Rows: Time instances
	% Columns: Swimmers
    meanx(pp,:)=mean(squeeze(xdat(:,:,pp))'); % Compute mean at all time instances
    varx(pp,:)=var(squeeze(xdat(:,:,pp))'); % Variance

    varx_eop(pp,:)=var(squeeze(xdat(1:8:end,:,pp))'); % Variance, end-of-period based
    varref=datt(pp).midpointvar; % Variance halfway through the simulation
    tref=datt(pp).midpointtime; % The corresponding time
    timevec_new=timevec(1:8:end); % Time vector containing the end-of-period values only
    dispx_mid(pp,:)=(varx_eop(pp,end-8:end)-varref)./(2*(timevec_new(end-8:end)-tref));
	% Dispersion

    tmix(pp)=datt(pp).SNRt; % Mixing time
    clear yall;
end

%% Drift computation
ft = fittype( 'poly1' ); % Fitting a line
ii=0;
for Wo=Wos % Going through Wos
    ii=ii+1;
    period=2*pi/(Wo*Wo*Sc);
    ncycles=ceil(timeabs/period);
    tfinal=ncycles*period;
    tstep=period/8; % Period computation

    opts = fitoptions('Method', 'LinearLeastSquares');
    xfit=tfinal-5*period+tstep:period:tfinal; % Fitting done for the last 5 periods
    mean_x=meanx(ii,end-5*8+1:8:end); % Evaluate mean position for each period
    yfit=mean_x;
  
    driftfit=fit(xfit',yfit',ft,opts); % Carry out the fit
    driftcoef(ii)=driftfit.p1; % Save the drift
end

clear datt % Clear datt as it is a very large file
clear xdat ydat ydat_bottom ydat_top ii
save('PP_GYRO_OSC.mat') % Save the processed data