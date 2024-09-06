%{
-------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 8/29/2021    
    Last update date: 7/23/2023

Description: this script 'diffAnalysis'
1) calculates the MSD( tau) of all tracks, of which the mean gives the EA-MSD
2) calculates the TA-MSD( tau) value for each track
3) perform linear fitting of squared displacment of individual tracks to
get apparent D coefficient & localization error (sigma)  
        MSD = 4Dt + 4* sigma^2 - 4D*dt/3
4) perform fitting of TA-MSD of individual tracks to get Dalpha & Alpha
        MSD = 4* {Dalpha}* t^{Alpha}
5) appends the values of each track to tracksFinal as a new field
6) output EnsMSD & EnsTAMSD in a matrix form

======Input=======
tracksFinal: U-track output

======Output=======
in-place modified tracksFinal (new field: MSD, Diff, LocErr, alpha, Dalpha)
EnsMSD:  MSD (M×N matrix), for all tracks (M) with maxTau=50 (N)
EnsTAMSD: time-averaged MSD for all tracks (M) with maxTau=50 (N)

-------------------------------------------------------------
%}

%% Calculating the MSD and fit the D & alpha
tic

maxTau = 50;            % the longest time lag used in calcuating MSD

nTracks = length( tracksFinal); % number of tracks

% initialization of the quantities with the correct size
EnsMSD = nan( nTracks, maxTau); 
EnsTAMSD = nan( nTracks, maxTau);

tracksFinal( nTracks).MSD = []; % TA-MSD for each track
tracksFinal( nTracks).Diff = [];
tracksFinal( nTracks).LocErr = [];
tracksFinal( nTracks).Alpha = [];
tracksFinal( nTracks).Dalpha = [];

for i = 1: nTracks
    
    traj = tracksFinal(i).tracksCoordXY* pixelSize; % store x and y coordinates
    nFrames = length( traj); % number of frames in each track
    
    % ~~~~ Ensemble-Averaged MSD ~~~~~
    % calculate first maxTau steps, assign NaN for shorter tracks that
    % doesn't have these steps, use mean(MSD,'omitnan') when plotting EAMSD
    MSD = nan( 1, maxTau);
    cutT = min( nFrames, maxTau+1); % cutoff of time maxTau+1, faster computation
    dr = traj( 2:cutT,:) - traj( 1,:); % displacement from origin with respect to time
    MSD( 1: cutT-1) = sum( dr.^2, 2); % squared displacement (t-1 timelags)
    EnsMSD( i,:) = MSD; % store it in ensemble MSD matrix (fixed length) for EA-MSD plotting
    
    % ~~~~ Time-Averaged MSD ~~~~~
    % each track has a TA-MSD with 'nFrame-1' data points
    taMSD = nan( 1, maxTau);
    for tau = 1: nFrames-1 % time gap (step size) from 1 to nFrames-1
        dr = traj( tau+1: end, :) - traj( 1: end-tau, :); % displacment of pairs with this timelag
        dr2 = sum( dr.^2, 2);
        taMSD( tau) = mean( dr2, 'omitnan'); % sum over all pairs with this timelag tau
    end
    
    EnsTAMSD( i,:) = taMSD( 1: maxTau); % store it in ensemble time-averaged MSD (fixed length) for EATA_MSD plotting
    tracksFinal(i).MSD = taMSD( 1: nFrames-1); % stored time-averaged MSD data into tracksFinal.MSD
    
    % fit the TA-MSD of individual tracks MSD = 4Dt + 4*locErr^2 - 4D*dt/3
    fitRange = 1: 3; % points of MSD used for fitting
    p = polyfit( timeStep* fitRange, taMSD( fitRange), 1); % linear fit, y = ax + b
    tracksFinal(i).Diff = p(1)/ 4;  % may be negative
    tracksFinal(i).LocErr = sqrt( p(2)+ timeStep*p(1)/3)/ 2; % unit: m, may be imaginary (should exclude later)    
    
    % fit the log(MSD) to get Dalpha & alpha
    fitRange = 1: 3; % max(3, floor( nFrames/4)); % 1/4 of the track
    pLog = polyfit( log( timeStep* fitRange), log( taMSD( fitRange)), 1); % log(MSD) = alpha* log(t) + log(4D)
    tracksFinal(i).Alpha = pLog(1);
    tracksFinal(i).Dalpha = exp( pLog(2))/ 4;
end

EnsMSD = EnsMSD* 1e12;      % unit: um^2
EnsTAMSD = EnsTAMSD* 1e12;  % unit: um^2

toc
fprintf( '    Diffusion Analysis (MSD, Diff, LocErr, alpha, Dalpha) Done\n');

