
function [ tracksSteps, tracksStepsLongi, tracksRg, tracksStepsNorm, tracksTheta] = stepToVectors(tf, cellInfo, pixelSize)
%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 1/12/2022
    Last update date: 7/23/2023

this script compute the single-step displacement and related quantities
i.e. radius of gyration Rg, step orientation angle theta, normalized
position of steps, step length projected to the long axis of cells
%}

    nTracks = length( tf);
    cellNum = [tf.ModCellNum]';
    
    % initialization, unit: um
    tracksRg = nan( nTracks, 1); % radius of gyration
    tracksSteps = cell( nTracks, 1);    tracksStepsLongi = cell( nTracks, 1);
    tracksTheta = cell( nTracks, 1);    tracksStepsNorm = cell( nTracks, 1);    
    
    for i = 1: nTracks

        traj = tf(i).tracksCoordXY; % dimension: nSteps*2
        
        % calculate single-step displacement
        singleSteps = sqrt( sum(( traj( 2:end, :)- traj( 1:end-1, :)).^2, 2)); % unit: pixel        
        tracksSteps{i} = singleSteps* pixelSize* 1e6; % unit: um 
        
        % calculate the mean radius of gyration Rg (RMSD from traj center), unit: um
        tracksRg(i) = sqrt( mean( sum(( traj- mean( traj)).^2, 2)))* pixelSize* 1e6;
        
        % Calculate normalized position of all steps [LNorm, xNorm]
        posNorm = tf(i).spotPosNorm;
        tracksStepsNorm{i} = posNorm( 1:end-1, :); % assign posNorm from where it jumps
        
        % Calculate the step orietation by the angle of steps (dirty way)
        % theta = arctan( d(width*xNorm/2)/ d(length*LNorm)
        c = cellInfo( cellNum(i));
        dStepNorm = posNorm( 2:end, :)- posNorm( 1:end-1, :); % [delta(LNorm), delta(xNorm)]        
        tracksTheta{i} = atand( c.width*dStepNorm(:,2)/2 ./ ( c.length*dStepNorm(:,1)));
        
        % calculate the step length projected to the longitudinal direction, unit: um
        tracksStepsLongi{i} = c.length* dStepNorm(:,1)* pixelSize* 1e6; 
        
%         tracksStepLongi{i} = diff( tf(i).spotPosInCell(:,1))* pixelSize* 1e6; % step projected to long axis
%         tracksStepTrans{i} = diff( tf(i).spotPosInCell(:,2))* pixelSize* 1e6; % step projected to short axis
    end
      
% disp( '~~~ StepToVectors Done ~~~')

end
