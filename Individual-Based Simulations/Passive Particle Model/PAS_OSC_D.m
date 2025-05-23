% Part of the supplementary material for

% Caldag, H.O. & Bees, M.A. Fine-tuning the dispersion of active suspensions using oscillatory flows.

% This code piece simulates passive particles under oscillating flow. The
% diffusion is represented in space and there is no orientation tracking.
% This allows direct comparison with the results from Lee et al., see Fig.
% 2a in the paper.

% The code structure is similar to the one used for gyrotactic particles,
% except for the diffusion modelling. Postprocessing to obtain dispersion
% is done within the code. See GYRO_OSC_sim.m for more comments.

SAVE=1;
swcount=5000;
Pe=12.8;
Sc=10;

timeabs=80; % Simulation duration (minimum)
Wos=logspace(-1.1549,0.6990,32);

parfor pp=1:length(Wos) % For all Wos (replace parfor with for if parfor is not usable)
    disp(num2str(Wos(pp)))

    x0=zeros(swcount,1);
    y0=-1+2*rand(swcount,1);

    Wo=Wos(pp);
    OHM=Wo.^2*Sc;
    period=2*pi/OHM;
    tstep=period/8;
    ncycles=ceil(timeabs/period);
    tfinal=ncycles*period; % Period and simulation duration computations

    velsc=velsccal(Wo,period);

    datt(pp)=simrun(tstep, timeabs, x0, y0,swcount,Pe,Wo,Sc,velsc);

end


%% Postprocessing
for ii=1:32 % Sweeping all Wos
    xd=datt(ii).yall(:,1:swcount); % Get axial coordinates
    Wo=Wos(ii)
    OHM=Wo.^2*Sc;
    period=2*pi/OHM;
    tstep=period/8;
    ncycles=ceil(timeabs/period);
    tfinal=ncycles*period; % Time-related definitions
    
    laststepcount=ceil(size(xd,1)/2); % Number of steps that could be recorded after the first 80 instances

    timevec=[[1:80]*tstep (tfinal-(laststepcount-1)*tstep):tstep:tfinal];

    varx=var(xd'); % Variance of coordinates

    disps(ii).dispx_start=0.5*(varx(2:end)-varx(1))./timevec(2:end); % Dispersion
    dispx_avg(ii)=mean(disps(ii).dispx_start(end-7:end)); % Period-averaging
    clear varx timevec xd

end
%%
Sc=10; % Schmidt number, the one used in Lee et al. and closest to our case
OHMs=Wos.^2*Sc; % Parameter relating time scales
loglog(OHMs./Sc, ((dispx_avg-1).*OHMs.^2./Pe.^2),'-o'); % The D_2D^{*} in Lee et al.

hold on;
lee_fig6_data
OHMlee=lee10(:,1).*Sc; % Plotted data in Lee et al.'s (2014) Fig. 6A
loglog(lee10(:,1),(lee10(:,2)),'-o')

if SAVE
    save("PAS_OSC_SPATIAL.mat",'-v7.3')
end

%%
function [RETURNDAT] = simrun(tstep, timeabs, x0, y0,swcount,Pe,Wo,Sc,velsc)

kk=0;
ss=1;
steps=round(timeabs/tstep)+1;

for tsolve=tstep:tstep:timeabs

    kk=kk+1;
    dt=(tstep/40).*(Wo>=1)+(tstep/200)*(Wo<1);

    if tsolve==tstep
        y=[x0;y0];
        yall(kk,:)=y';
    end

    for t=tsolve-tstep:dt:tsolve
        [y]=movefun(t,y,swcount,Pe,dt,Wo,Sc,velsc);
    end

    if kk<81 || kk>steps-80
        ss=ss+1;
        yall(ss+1,:)=y';
    end
    % yall(kk,:)=y';
end

RETURNDAT.yall=yall(:,1:2*swcount);

end

function [y]=movefun(t,y,swcount,Pe,dt,Wo,Sc,velsc)

D=1; % Diffusion coefficient set to 1

tmod=t*Wo.^2*Sc;
% velsc=1; % Do this to remove tidal volume fixing
Uflow=Pe*velsc*real(1/(1i*Wo^2)*( cosh(sqrt(1i)*Wo.*y(swcount+1:2*swcount))/cosh(sqrt(1i)*Wo) - 1  ) .*exp(1i.*(tmod) ));
Uflow=Uflow.';

randy=sqrt(2*D*dt)*randn(swcount,1); % Random diffusion in the y-direction
% Computed beforehand to evaluate if the particles will hit the wall

crossvec=(abs(randy+y(swcount+1:2*swcount))>=1); % Find the particles hitting the wall.

y=y+[sqrt(2*D/dt)*randn(swcount,1)+Uflow';...
        crossvec.*( sign(y(swcount+1:2*swcount))-(2*y(swcount+1:2*swcount)+randy-sign(y(swcount+1:2*swcount)) ) )./dt + ~crossvec.*(randy./dt)].*dt;

% Particle moves in the x- direction under the effect of diffusion and
% flow.

% In the y- direction, the particle moves under the effect of diffusion
% only. The particle is reflected back if it is hitting the wall.

end

function velsc=velsccal(Wo,period)
Sc=10;
yveccal=-1:0.01:1;tveccal=linspace(0,period,202)';
fycal=1/(1i*Wo^2)*( cosh(sqrt(1i)*Wo.*yveccal)/cosh(sqrt(1i)*Wo) - 1);
Woflowcal=real(fycal.*exp(1i.*tveccal.*Wo^2*Sc));
Woflow_avgcal=sqrt(1/(4*pi)*trapz(yveccal,trapz(tveccal.*Wo^2*Sc,Woflowcal.^2)));
velsc=1/Woflow_avgcal;
end

