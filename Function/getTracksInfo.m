
function [nTracks, totalCells, cellNum, tracksLength, tracksOrigin, tracksFilled, badMovie] = getTracksInfo( tf)
%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 1/12/2022
    Last update date: 7/14/2023
%}

    nTracks = length( tf);
    tracksLength = cellfun( @length, {tf.tracksFeatIndxCG}');
    tracksOrigin = [ tf.origin]'; 
    cellNum = [ tf.ModCellNum]';
    totalCells = max( cellNum);
    tracksFilled = [ tf.filled]';
    badMovie = [ tf.badMovie]';
end