
function [ tracksLNorm, tracksxNorm, tracksLNorm4, tracksxNorm4, tracksDiff] = dataToVectors( tf)
%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 11/20/2021
    Last update date: 7/15/2023

this script convert the LNorm & xNorm & Diffusion properties of
tracksFinal into  matrix/array, covenient for future plotting 
%}

    nTracks = length( tf);
    
    % initialization 
    tracksLNorm = nan( nTracks, 4);
    tracksxNorm = nan( nTracks, 2);    
    tracksLNorm4 = nan( nTracks, 4);
    tracksxNorm4 = nan( nTracks, 4);
    tracksDiff = nan( nTracks, 4);
    
    for i = 1: nTracks
        
        LNorm = tf(i).spotPosNorm(:,1);
        xNorm = tf(i).spotPosNorm(:,2);
        
        tracksLNorm(i,:) = [ LNorm(1), LNorm(end), min( LNorm), max( LNorm)];
        tracksxNorm(i,:) = [ xNorm(1), xNorm(end)]; % [first pt, end pt]
        
        tracksLNorm4(i,:) = LNorm(1:4);
        tracksxNorm4(i,:) = xNorm(1:4);
        
        % tracksDiff: diff (um^2/s), locErr (nm), alpha, dalpha (um^2/s)
        tracksDiff(i,:) = [ tf(i).Diff*1e12, tf(i).LocErr*1e9, tf(i).Alpha, tf(i).Dalpha*1e12];
    end
    
% disp( '~~~ DataToVectors Done ~~~')
end

