function tf = reNumCells_yh(tf)
%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 7/28/2022
    Last update date: 7/14/2023

Description: this function re-order the cell number in the combined
tracksFinal file by counting from the cells, there should be no repetition
or gaps among cell numbers. It also returns the collection of cell mesh
information

Output:
tracksFinal implemented with new cell numbers stored in tf.ModCellNum
cellMeshAll: collection of all cellMesh (new order) for future trackAveNorm

Example: [1 3], [ 1 2 4] --> [1 2], [3 4 5]
%}

    dataOrigin = unique( string( {tf.origin}')); % string array of fileName
    R = size( dataOrigin, 1); % number of files       
    chosenCell = cell( R, 1); % chosed cells in each movie
    pathFlag = true;
    
    [~, fileNum] = ismember( string( {tf.origin}'), dataOrigin); % find fileNum for each track
    cellNum = [ tf.cellNumber]'; % original cellNumber for each track
    index = unique( [fileNum cellNum], 'rows'); % unique combination of fileNum & cellNum
    
    % count the appearance of certain movie & cell combination in a 2D array
    % i.e. which cells in which movies appeared
    cellNumCount2D = zeros( max( index)); 
    cellNumCount2D( sub2ind( size( cellNumCount2D), index(:,1), index(:,2))) = 1;
    
    % re-align the cell numbers so no gaps in between, make a continous list
    cellNumBefore = [0; cumsum( sum(cellNumCount2D, 2))]; % vector of total cell numbers appeared before this file
    cellOrder2D = cumsum( cellNumCount2D, 2); % re-align in each movie
    cellNumNew = cellOrder2D + cellNumBefore( 1:end-1); % 2D array of new cellNumber for each [fileNum, cellNum]
    
    % assign the modified cellNum to the tracks by 2D indexing: [fileNum, cellNum]
    modCellNum = num2cell( cellNumNew( sub2ind( size( cellOrder2D), fileNum, cellNum)));         
    [tf.ModCellNum] = modCellNum{:};    
    
end
