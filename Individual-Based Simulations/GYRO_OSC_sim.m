% Part of the supplementary material for

% Caldag, H.O. & Bees, M.A. Fine-tuning the dispersion of active suspensions using oscillatory flows.

% This code piece is used to run individual-based simulations of gyrotactic particles.
% The particles are subject to oscillatory flow, characterized by the
% Womersley number.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output structure: Data is recorded in a struct called datt. For the i-th
% Womersley number in the array Wos, the output is stored in datt(i).
% datt(i).yall contains the x- and y- coordinates of all particles during
% the first and last 10 periods of oscillations. The rows denote the time
% steps (each 1/8th of a period apart within the first and last 10 periods)
% while the j-th column denotes swimmer j. x and y coordinate data are
% stacked on top of each other, meaning that columns 1:swcount (swcount is
% the number of swimmers simulated) columns contain x- coordinates while
% columns swcount+1:2*swcount contain y- coordinates. Each datt(i) also
% contain a SNRt value, the time when the signal-to-noise ratio falls
% below 0.1.

clearvars;

SAVE=1; % Data recording Boolean

swcount=5000; % number of particles simulated
Pe=12.8; % Peclet number
Sc=16.8; % Schmidt number

AUGUSTAE=1; % Booleans to adjust the cell type simulated. AUGUSTAE=1 for C. augustae,
SALINA=0;   % SALINA=1 for D. salina
STRONG=0;   % STRONG=1 for strongly gyrotactic particle (see text).

if AUGUSTAE && ~SALINA && ~STRONG
	Bdim=3.4; % Gyrotactic reorientation time, chlamydomonas augustae
	drdim=0.067; % Rotational diffusivity, chlamydomonas augustae
elseif ~AUGUSTAE && SALINA && ~STRONG
	Bdim=10.5; % Gyrotactic reorientation time, chlamydomonas augustae
	drdim=0.23; % Rotational diffusivity, chlamydomonas augustae
elseif ~AUGUSTAE && ~SALINA && STRONG
	Bdim=0.3; % Gyrotactic reorientation time, chlamydomonas augustae
	drdim=0.01; % Rotational diffusivity, chlamydomonas augustae
else
	error('Set only one of AUGUSTAE, SALINA or STRONG values equal to 1, the rest should be zero.')
end

timeabs=1500; % This is the absolute minimum of the simulation duration
% Exact duration depends on Wo but it is always larger than or equal to timeabs.
% We always simulate an integer number of periods
Wos=logspace(-1.1549,0.6990,32); % Range of Womersley numbers tested

parfor pp=1:length(Wos) % parfor loop is used to speed up the simulations. Replace with for if parfor is unavailable
    disp(num2str(Wos(pp))) % Display the Wo simulated

    x0=zeros(swcount,1);
    y0=-1+2*rand(swcount,1); % Initial distribution, particles at x=0, spread across the channel.

    y0(1:swcount/2)=-1+rand(swcount/2,1); % We set the first half of the particles exclusively to one side of the channel to track their 
    y0(swcount/2+1:swcount)=rand(swcount/2,1); % signal-to-noise ratio (SNR) easily.

    thet0=-pi+2*pi*rand(swcount,1); % Initial orientation
    pvec0=[sin(thet0)';cos(thet0)']; % Orientation turned into a vector

    Wo=Wos(pp); % The Wo used for the current iteration of the loop
    OHM=Wo.^2*Sc; % OHM relates the oscillation and cell time scales
    period=2*pi/OHM; % Period in cell-based scaling
    tstep=period/8; % Time stepping while recording data
    ncycles=ceil(timeabs/period); % Number of cycles
    tfinal=ncycles*period; % The actual duration of the simulation (a value >= 1500, complete number of periods)

    velsc=velsccal(Wo,period); % Velocity scaling ensures the root-mean-square velocity ends up equal to Pe.

    Vs=1; % Swimming velocity

    datt(pp)=simrun(tstep, tfinal, x0, y0,pvec0, Vs,swcount,Pe,Wo,Sc,period,velsc,Bdim,drdim); % Send it to running code

end

if SAVE % Saving data
    save("GYRO_OSC.mat",'-v7.3')
end


% simrun function runs through the steps of simulation
function [RETURNDAT] = simrun(tstep, tfinal, x0, y0,pvec0, Vs,swcount,Pe,Wo,Sc,period,velsc,Bdim,drdim)

kk=0; % Loop counter (time steps)
ss=1; % Loop counter (only for recorded instances)
% We record only the first and last 10 periods of motion to save space.

steps=round(tfinal/tstep)+1; % Number of time steps
ncycles=tfinal/period; % The number of cycles
tmid=floor(ncycles/2)*period; % Time corresponding to halfway through a simulation
% Will be used to record variance data at that time.

steps=round(tfinal/tstep)+1; % Number of time steps to be taken

    SNRFALL=0; % Boolean checking if the SNR value has fallen below the reference value of 0.1
    SNRt=0; % The time the SNR value falls below 0.1, set as 0 initially.

for tsolve=tstep:tstep:tfinal % Going through the time steps

    kk=kk+1;
    dt=(tstep/10).*(Wo>=1)+(tstep/50)*(Wo<1); % We take smaller steps depending on the Womersley number

    if tsolve==tstep % If this is the first step
        y=[x0;y0];
        yall(kk,:)=y'; % Record the initial condition
        pvec=pvec0; % Set the orientation to be the initial orientation
    end

    for t=tsolve-tstep:dt:tsolve % Take small time steps in between the larger time steps where we record data
        [y,pvec]=movefun(t,y,pvec,Vs,swcount,Pe,dt,Wo,Sc,period,velsc,Bdim,drdim); % Simulate dt amount of swimming
        SNRval=mean(y(1.5*swcount+1:2*swcount))./std(y(1.5*swcount+1:2*swcount)); % Evaluate the SNR
        if SNRval<0.1 && SNRFALL==0 % If SNR has fallen below the reference, save the time it happened.
            SNRt=t;
            SNRFALL=1;
        end
    end

    if kk<81 || kk>steps-80 % Recording only the first and last 10 periods
        ss=ss+1;
        yall(ss,:)=y';
    end
	
	if abs(t-tmid)/tmid<1e-3 % If the current time is very close to halfway point
        midpointvar=var(y(1:swcount)); % Evaluate the variance of x- data
        midpointtime=t; % Get the exact time
    end
	
end

RETURNDAT.yall=yall(:,1:2*swcount);
RETURNDAT.SNRt=SNRt; % When the simulation ends, return the position data and SNRt in a struct.
RETURNDAT.midpointvar=midpointvar; % Also return the variance and time halfway through the simulation.
RETURNDAT.midpointtime=midpointtime;

end


% movefun is the function where the particle motion happens. Returns the position and orientation at the next time step.
function [y,pwithnoise]=movefun(t,y,pvec,Vs,swcount,Pe,dt,Wo,Sc,period,velsc,Bdim,drdim)

dr=1;Vs=1; % Both set to 1.
B=Bdim*drdim/Vs^2; % Non-dimensional reorientation rate

tmod=t*Wo.^2*Sc; % Time in oscillation terms
% % velsc=1; % Uncomment this line to remove tidal volume fixing
Uflow=Pe*velsc*real(1/(1i*Wo^2)*( cosh(sqrt(1i)*Wo.*y(swcount+1:2*swcount))/cosh(sqrt(1i)*Wo) - 1  ) .*exp(1i.*(tmod) )); % Oscillating flow
delU=-Pe*velsc*imag(1*(2^(1/2)*sinh(2^(1/2)*   Wo*y(swcount+1:2*swcount)*(1/2 + 1i/2)).*(cos(tmod) + sin(tmod)*1i)*(1/2 + 1i/2))...
    /(Wo*cosh(2^(1/2)*Wo*(1/2 + 1i/2)))); % Derivative (for vorticity)
Uflow=Uflow.';delU=delU.';

kmat=repmat([-1;0], 1, swcount); % Vertical direction
dp_gyro=(1/(2*B)*(kmat-dot(kmat,pvec).*pvec)).*dt; % Gyrotactic contribution
vormat=[zeros(size(delU));zeros(size(delU));delU]; % Vorticity in matrix form
pvec3=[pvec;zeros(1,swcount)]; % 3-dimensional form of the orientation vectors
dp_vor=(1/2*cross(vormat,pvec3)).*dt; % Vorticity-related contribution to orientation
dp_vor(3,:)=[]; % Remove third component (system is 2-dimensional)
pnoise=randn(2,swcount); % Noise component to orientation
pwithnoise=pvec+Vs*sqrt(2*dt).*pnoise+dp_vor+dp_gyro; % Add them all up to find the new set of orientation vectors
pwithnoise=pwithnoise./(vecnorm(pwithnoise));

crossvec=(abs(Vs*pwithnoise(2,:)*dt+y(swcount+1:2*swcount)')>=1); % These are the particles that attempt to cross the boundary, will be reflected off

% x component changes with the flow and the swimming
% y component changes only with swimming. Reflection is addressed as well.
y=y+[Vs*pwithnoise(1,:)'+Uflow';...
        crossvec'.*( sign(y(swcount+1:2*swcount))-(2*y(swcount+1:2*swcount)+Vs*pwithnoise(2,:)'*dt-sign(y(swcount+1:2*swcount)) ) )./dt + ~crossvec'.*(Vs*pwithnoise(2,:)')].*dt;

pwithnoise(2, crossvec)=-pwithnoise(2, crossvec); % Reflecting off the particles in terms of their orientation

end

% Function velsc computes the root-mean-square velocity and returns a velsc value such that the root-mean-square of the
% field is unity. This way, when the field is multiplied with Pe, we know that the rms-based definition of Pe is equal to-noise
% the value we prescribe.
function velsc=velsccal(Wo,period)
Sc=16.8;
yveccal=-1:0.01:1;tveccal=linspace(0,period,202)';
fycal=1/(1i*Wo^2)*( cosh(sqrt(1i)*Wo.*yveccal)/cosh(sqrt(1i)*Wo) - 1);
Woflowcal=real(fycal.*exp(1i.*tveccal.*Wo^2*Sc));
Woflow_avgcal=sqrt(1/(4*pi)*trapz(yveccal,trapz(tveccal.*Wo^2*Sc,Woflowcal.^2)));
velsc=1/Woflow_avgcal;
end

