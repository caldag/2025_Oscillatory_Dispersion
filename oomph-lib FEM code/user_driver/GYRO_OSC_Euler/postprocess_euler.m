% Part of the supplementary material for

% Caldag, H.O. & Bees, M.A. Fine-tuning the dispersion of active suspensions using oscillatory flows.

% This code processes the output files from the oomph-lib simulations and
% saves the concentration field and coordinates of the meshing in a grid
% form that can be interpreted by Matlab.

% Raw data coming out of oomph-lib simulations are not supplemented here due to
% their excessive size. Please run the simulations first.

% This script should be put into the same folder as the oomph-lib code.
% It is assumed the data for each simulation are saved into folders
% starting with "RESLT" and followed by a number corresponding to the index
% of the Wo number in the main simulation loop.

% natsortfiles is required to run the code, can be installed via the link

% https://uk.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort

clearvars;close all;
RESstruct=natsortfiles(dir('RESLT*')); % Find folders with the name RESLT
RESlength=length(RESstruct); % Get the number of folders

Wos=[0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.15 0.25 0.35 0.45 0.55 0.95 1 1.1 1.2 1.05 1.5];
% These Wo numbers are the simulated cases in the article.
% Update as necessary, should match the available data.

for rr=1:RESlength % Go through all folders
    resfolder=RESstruct(rr).name % Display folder name studied
    cd(resfolder); % Move to the result folder
    fstruct = dir('soln*.dat'); % Get the names of all solution files.
    data_length=length(fstruct)-1; % The exact length

    parfor ii=1:1:data_length+1 % parfor here to speed things up. If you cannot use parfor, replace with for
        filefind=dir(strcat(['soln' num2str(ii-1) '_t=*.dat']));
        filename=filefind.name;
        tstr=filename([strfind(filename,'=')+1:strfind(filename,'.dat')-1]); % Get the time stamp information
        timevec(ii+1)=str2double(tstr); % Record the time in the time vector
        
        readdata=fileread(filename); % Read the contents
        readdata=regexprep(readdata,'ZONE I=5, J=5',''); % Remove these lines with strings
        aa=str2num(readdata); % Convert text to numbers

        rangevec=1:size(aa,1); % Number of rows corresponds to the range 

        % The element and domain sizes are determined from the data
        % structure
        % First column lists the x- components, second column lists the y-
        % components. Each element is divided into a 5x5 grid, this is why
        % the element on the 25th row correspond to the element sizes in
        % the x- and y- directions.

        % Note that the x- direction here corresponds to the cross-channel
        % direction (y- direction in the article) and y- direction here 
        % corresponds to the axial direction (x- direction in the article) 

        xxr=0:aa(25,1):aa(end,1);
        yyr=0:aa(25,2):aa(end,2);

        [xx(:,:,ii+1),yy(:,:,ii+1)]=meshgrid(xxr,yyr);
        zz(:,:,ii+1)=griddata(aa(rangevec,1),aa(rangevec,2),aa(rangevec,3),squeeze(xx(:,:,ii+1)),squeeze(yy(:,:,ii+1)));

        % The operations above form a mesh grid and matches the
        % concentration values in the grid
        
        zz_vel(:,:,ii+1)=griddata(aa(rangevec,1),aa(rangevec,2),aa(rangevec,end),squeeze(xx(:,:,ii+1)),squeeze(yy(:,:,ii+1)));
        % Velocity data

    end
    save(strcat(['PP_Wo_' num2str(Wos(rr)) '.mat']), 'xx','yy','zz','zz_vel')
end