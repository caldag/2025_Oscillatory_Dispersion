% Part of the supplementary material for

% Caldag, H.O. & Bees, M.A. Fine-tuning the dispersion of active suspensions using oscillatory flows.

% This script computes the drift and dispersion coefficients from the
% Eulerian simulation data. After running the simulation and running
% postprocess_euler.m, run this script to extract drift and
% dispersion coefficients. The script should be in the same folder with the
% data recorded by postprocess_euler.m

clearvars;
RESstruct=natsortfiles(dir('PP*Wo*.mat')); % Find files with the matching names
RESlength=length(RESstruct); % Number of .mat files.
Wos=[0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.15 0.25 0.35 0.45 0.55 0.95 1 1.1 1.2 1.05 1.5];
% These Wo numbers are the simulated cases in the article.
% Update as necessary, should match the available data.

Sc=16.8; % Schmidt number
timeabs=300; % This is the simulation duration in terms of cell time scale
for rr=1:RESlength % Go through each result file.
    Wo=Wos(rr); % The Wo number for this case
    OHM=Wo^2*Sc; % Wo^2*Sc group relates time scales
    period=2*pi/OHM; % Period in oscillatory time scale
    ncycles=ceil(timeabs/period); % Number of cycles rounded up
    timeact=period*ncycles; % The actual simulation duration is rounded up as well
    tstep=period/8; % Time step between the recorded instances

    if ncycles>=20
        timevec=[0:80 ncycles*8-80:ncycles*8]*tstep;
    else
        timevec=0:tstep:timeact;
    end % timevec contains instances simulation data was recorded
    % The first and last 10 periods are recorded for each simulation.
    % If by chance the number of periods is less than 20, the simulation
    % recorded all periods.

    load(RESstruct(rr).name); % Load file
    
    yyr=0:yy(2,1,end):yy(end,1,end); % Assign x and y meshes
    xxr=0:xx(1,2,end):xx(1,end,end);
    
    for ii=2:1:size(xx,3) % Start from the second time step to evaluate measures
        thiscon=zz(:,:,ii);
        thisx=xx(:,:,ii)-1;
        thisy=yy(:,:,ii); % Concentration fields, x and y meshes for the step
        maxcon=max(max(thiscon)); % Find maximum concentration
        [meanpl]=find(thiscon>0.01*maxcon); % Find places where concentration is above 1% of the maximum
        weightedmeanx(ii)=mean(sum((thiscon(meanpl)).*thisx(meanpl))./sum(thiscon(thiscon>0.01*maxcon)));
        weightedmeany(ii)=mean(sum((thiscon(meanpl)).*thisy(meanpl))./sum(thiscon(thiscon>0.01*maxcon))); 
        % Evaluate weighted mean positions

        % Note that the x- direction here corresponds to the cross-channel
        % direction (y- direction in the article) and y- direction 
        % corresponds to the axial direction (x- direction in the article) 

        out(rr).vars(ii)=sum(sum((zz(:,:,ii).*(yy(:,:,ii)-weightedmeany(ii)).^2 )))./sum(sum(thiscon(thiscon>0.01*maxcon)));
        % Evaluate variance
    end

    for ii=130:1:size(xx,3) % Evaluate drift and dispersion
        out(rr).drifts(ii)=(weightedmeany(end)-weightedmeany(2))./(timeact);
        out(rr).disps(ii)=0.5*(out(rr).vars(ii)-out(rr).vars(2))./(timeact);
    end
    % We use the whole span of data as the simulation duration is short

    drifts(rr)=out(rr).drifts(end);
    disps(rr)=out(rr).disps(end);
    clear xx yy zz zz_vel xxr yyr thisx thisy timevec weightedmeanx weightedmeany timeact
end
