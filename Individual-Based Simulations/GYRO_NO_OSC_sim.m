% Part of the supplementary material for

% Caldag, H.O. & Bees, M.A. Fine-tuning the dispersion of active suspensions using oscillatory flows.

% This code piece is used to run individual-based simulations of gyrotactic particles.
% The difference from GYRO_OSC_sim.m is that we simulate purely downwelling
% or upwelling flows.

clearvars;

SAVE=1; % Data recording Boolean

swcount=5000; % number of particles simulated
Pe=12.8; % Peclet number
Sc=16.8; % Schmidt number
Bdim=3.4; % Gyrotactic reorientation time, chlamydomonas augustae
drdim=0.067; % Rotational diffusivity, chlamydomonas augustae

timeabs=1500; % This is the absolute minimum of the simulation duration
% Exact duration depends on Wo but it is always larger than or equal to timeabs.
% We always simulate an integer number of periods
dirs=[-1 1]; % Flow directions

parfor pp=1:length(dirs) % parfor loop is used to speed up the simulations. Replace with for if parfor is unavailable
    disp(num2str(dirs(pp))) % Display the direction of the flow

    x0=zeros(swcount,1);
    y0=-1+2*rand(swcount,1); % Initial distribution, particles at x=0, spread across the channel.

    y0(1:swcount/2)=-1+rand(swcount/2,1); % We set the first half of the particles exclusively to one side of the channel to track their 
    y0(swcount/2+1:swcount)=rand(swcount/2,1); % signal-to-noise ratio (SNR) easily.

    thet0=-pi+2*pi*rand(swcount,1); % Initial orientation
    pvec0=[sin(thet0)';cos(thet0)']; % Orientation turned into a vector

    dir=dirs(pp); % The Wo used for the current iteration of the loop
    tstep=0.1; % An arbitrary low value to determine data recording intervals
    tfinal=timeabs; % No oscillations mean that we can use timeabs as the duration of the simulation

    velsc=1; % No velocity scaling.

    Vs=1; % Swimming velocity

    datt(pp)=simrun(tstep, tfinal, x0, y0,pvec0, Vs,swcount,Pe,dir,Bdim,drdim); % Send it to running code

end

if SAVE % Saving data
    save("GYRO_NOOSC.mat",'-v7.3')
end


% simrun function runs through the steps of simulation
function [RETURNDAT] = simrun(tstep, tfinal, x0, y0,pvec0, Vs,swcount,Pe,dir,Bdim,drdim)

kk=0; % Loop counter (time steps)
ss=1; % Loop counter (only for recorded instances)
% We record only the first and last 8 seconds of motion to save space.

steps=round(tfinal/tstep)+1; % Number of time steps to be taken

for tsolve=tstep:tstep:tfinal % Going through the time steps

    kk=kk+1;
    dt=(tstep/50); % Fixed dt

    if tsolve==tstep % If this is the first step
        y=[x0;y0];
        yall(kk,:)=y'; % Record the initial condition
        pvec=pvec0; % Set the orientation to be the initial orientation
    end

    for t=tsolve-tstep:dt:tsolve % Take small time steps in between the larger time steps where we record data
        [y,pvec]=movefun(t,y,pvec,Vs,swcount,Pe,dt,dir,Bdim,drdim); % Simulate dt amount of swimming
        % SNR computations omitted
    end

    if kk<81 || kk>steps-80 % Recording only the first and last 10 periods
        ss=ss+1;
        yall(ss+1,:)=y';
    end
end

RETURNDAT.yall=yall(:,1:2*swcount);

end


% movefun is the function where the particle motion happens. Returns the position and orientation at the next time step.
function [y,pwithnoise]=movefun(t,y,pvec,Vs,swcount,Pe,dt,dir,Bdim,drdim)

dr=1;Vs=1; % Both set to 1.
B=Bdim*drdim/Vs^2; % Non-dimensional reorientation rate

Uflow=Pe*dir*(1-y(swcount+1:2*swcount).^2);
delU=2*Pe*dir*y(swcount+1:2*swcount);
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