function [cellInfo, poleBounds] = getCellInfo(tf, cellMesh, pixelSize)
%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 11/21/2021
    Last update date: 7/14/2023

this function returns the info (width, length, etc) of all cells & a vector
poleBounds for the cap region
%}

    cellNum = [ tf.ModCellNum]'; % renumbered cell number
    totalCells = max( cellNum);
    
    cellInfo( totalCells, 1) = struct; % cellInfo structure for all cells in the tf
    poleBounds = nan( totalCells, 1); % bounds for cap (pole) region

    for i = 1: totalCells
        
        thisCell = find( cellNum == i);
        cellNumOri = tf( thisCell(1)).cellNumber; % original cell number in the image
        cellFileName = char( tf( thisCell(1)).origin); % filename of the image
        
        cellInfo(i).origin = cellFileName;
        cellInfo(i).cellNumOri = cellNumOri;
        
        if cellNumOri == 0 || isnan( cellNumOri)
            warning( '~~~ No cell found for ModCellNum = %d', i);
            continue
        end
        
        meshOut = cellMesh(i).meshOut;
        cellInfo(i).area = polyarea( meshOut(:,1), meshOut(:,2))* pixelSize^2;
        cellInfo(i).length = sum( cellMesh(i).gridLen)* pixelSize;
        cellInfo(i).width = cellInfo(i).area/ cellInfo(i).length; % this is an estimation
        cellInfo(i).mesh = cellMesh(i).mesh;
        
        % define cell pole: cell width/2
        poleBounds(i) = cellInfo(i).width/ 2/ cellInfo(i).length; 
    end
    
%     disp('~~~ Cell Info Got ~~~');
end
