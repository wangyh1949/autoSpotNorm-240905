 %{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    creation date: 9/1/2021
    Last update date: 11/17/2023

Description: autoSpotNorm will automatically run spotNorm for all cell
meshes & tracking movies in the folder, and combine them into the same
tracksFinal variable. It will then calcuate quantities for analysis (e.g.
MSD, Diff, alpha). Eventually, these quantities are converted into matrices
(arrays) which are convenient for future plotting

10/24/2023: version for SPT experiments, only saving combined Matrix file
and tf_Variables files, also save combined tf 

======Input=======
    BF00xm.mat (output mesh file of Oufti)
    Tracking00x (the output folder of Utrack)

======Output=======
    cell & track overlay image (from spotNorm)
    _Variables.mat for each movie (from spotNorm)
    tracksFinal (combined tf, no analysis, Archive)
    tracksFinal All (with all analysis fields)
    Matrix & Steps file

=====Functions=====
    spotNorm_yh.m  (major)
    reNumCells_yh.m
    getCellInfo.m
    getTracksInfo.m
    diffAnalysis.m
    dataToVectors.m
    stepToVectors.m

*** check out spotNorm_yh.m for more details ***

%}

%% Step 0. Set up the parameters & flags

clear, clc, close all

dateAdd = ' autoSB 3 pix';    % used for: Date = [expDate dateAdd];

pixelSize = 160e-9;     % Andor camera, unit: m
timeStep = 21.742e-3;   % frame time, unit: s  [check the movie output], needed in MSD fitting

imgSaveFlag = true;     % flag for saving the mesh & track overlay images

% flag for oufti cells (are all cells picked by oufti isolated?)
%     true: assume all cells are good (isolated) and autorun; 
%     false: manually input cell number for each BF image 
goodCellFlag = true;

% flag for combining tracksFinal
%     true: automatically combine tracks in all good (isolated) cells
%     false: manually select cell and only combine tracks in those cells
autoCombineMovieFlag = true; 

% set up the storage path for analysis result, users should set up their own path
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';


%% Step 1. Automatically run spotNorm for all movies in the folder

disp('~~~ Please Choose BF Mesh File ~~~');
[meshName, meshPath] = uigetfile( '*.mat', 'Please Select a BF mesh .mat file');
fprintf( '\n~~~ We are looking at %s ~~~ \n\n', meshPath);     cd( meshPath)
    
if ~exist( [meshPath 'Mesh & Tracks Image'], 'dir')
    mkdir( [meshPath 'Mesh & Tracks Image']);
    mkdir( [meshPath 'Mesh & Tracks Image\auto shiftBack']);
end

if logical( imgSaveFlag)
    fprintf( '~~~ Img save flag is on, cell & track overlay image will be saved ~~~\n\n')
else
    fprintf( '~~~ Img save flag is off, cell & track overlay image will NOT be saved ~~~\n\n')
end

% criteria for mostly-filled cells, needed in spotNorm_yh (autoShiftBack module) 
%       0.6 for membrane bound strain;  0.5 for cytoplasmic strain 
%       lower when necessary
cytoFlag = logical( input( '       Is this strain Cytoplasmic (1:yes, 0:no)?  ')); fprintf( '\n')
if logical( cytoFlag)
    areaRatioThresh = 0.4; % Cytoplasmic: 0.4
else
    areaRatioThresh = 0.5; % Membrane: 0.5
end


% find all BF .mat files in this folder [oufti output for cell mesh]
list = dir( [ meshPath  '*BF*m.mat']);
tStart = tic;

% run spotNorm for all files in this folder
for i = 1: length(list) 
    
    meshName = list(i).name;
    
    % run spotNorm and save _Variables file under the strain folder
    savePath = spotNorm_yh( meshPath, meshName, varPath, areaRatioThresh, goodCellFlag, imgSaveFlag);
    % savePath: [ varPath 'Variables\' strain '\' expDate extraName '\']; 
end

tEnd = toc(tStart);
fprintf( '~~~~~~~ Hello World, spotNorm All Done, Total Time: %.1fs ~~~~~~~\n\n\n', tEnd);


%% Step 2. Combine the tracksFinal & cellMesh of multiple movies into one file

clearvars tfComb
cellMeshAll = [];   cellRecord = {};  extraName = ''; 

if ~exist( 'savePath', 'var')
    disp('~~~ Please choose the spotNorm _Variables file ~~~');
    [~, savePath] = uigetfile( '*.mat', 'Please Choose a spotNorm _Variables file', varPath);
end

path = savePath;
list = dir( [savePath '*_Variables.mat']); % list files that match name
load( [savePath list(1).name])  % load some basic info of the strain & expDate
fprintf( '~~~ Cell & analysis info loaded, start combining tracksFinal ~~~\n\n')
    
for i = 1: length(list)  % you can also choose a subset of the list

%     savePath = [ varPath 'Variables\' strain '\' expDate extraName '\'];
    savePath = path;
    
    % load _Variables.mat file of individual movie
    clearvars tracksFinal_InIsolatedCell
    load( [savePath list(i).name], 'tracksFinal_InIsolatedCell', 'cellNum',...
        'folderName', 'meshName', 'meshPath', 'cellList', 'cellMesh',...
        'fileName', 'expDate', 'cellFilled', 'goodCells', 'badMovie')    
    fprintf( '  - %s loaded ~~~\n',  fileName);    

    % we only take tracks in isolated cells
    tfGood = tracksFinal_InIsolatedCell;

    % store filename & badMovie flag in new fields of the tfGood structure
    nTracks = length( tfGood);          value = cell( nTracks, 1);    
    value(:) = { string( fileName)};    [tfGood.origin] = value{:};    
    value(:) = { badMovie};             [tfGood.badMovie] = value{:};    

    % store tracksFilled flag (tracks in filled cells) in tfGood
    tracksFilled = num2cell( ismember( cellNum, find( cellFilled)));
    [ tfGood.filled] = tracksFilled{:};

    if ~exist('tfComb', 'var') 
    % create tfComb with same fields as tfGood
        tfComb = tfGood;
        tfComb(:) = [];
    end

    if ~autoCombineMovieFlag % manually select the good cells based on overlay image
        fprintf( '     You can choose from %s\n',  join( string( goodCells), ' '))
        flag = input( '     Are all cells good? (0:Yes, -1:None, [x x]:certain cells): ');
        if flag == -1 % skip this movie
            continue
        elseif isnumeric(flag) % combine tracks from only certain cells
            if ismember( flag, goodCells) % is the input part of goodCells?
                goodCells = flag;
                tfGood = tfGood( ismember( cellNum, flag));
            else
                error('      You have to select from good cells with tracks in it')
            end
        end
    end

    % stack tfGood & cellMesh from different movies together
    tfComb = [tfComb; tfGood]; 
    cellMeshAll = [cellMeshAll; cellMesh( goodCells)]; 
    cellRecord = [cellRecord; {fileName goodCells folderName meshName meshPath}]; % {fileName, goodCells}

    fprintf( '    Included Good cells with tracks: %s\n', join( string( goodCells), ' '))
    if sum( ~cellFilled) > 0
        fprintf( '    Not mostly filled cells: %s\n\n', join( string( find( ~cellFilled')), ' '))    
    else
        fprintf( '    All cells are mostly filled ~~~\n\n')
    end
end

fprintf( '~~~~~~~ tracksFinal Combination Finished ~~~~~~~\n\n');

%% Step 3. reNumCell & Get cell & tracks info

dateFlag = logical( input( sprintf( '  ~~~ dataAdd is ''%s'', is it right (1:yes, 0:no)?  ', dateAdd)));
if ~dateFlag
    dateAdd = input( '      ~~~ What is the dateAdd (like '' autoSB 2 pix''): ', 's');
end

Date = [expDate dateAdd];
tracksFinal = reNumCells_yh( tfComb); % reorder the cell number, load Variables

pixelSize = 160e-9;     % Andor camera, unit: m
timeStep = 21.742e-3;   % frame time, unit: s  [check the movie output], needed in MSD fitting

[cellInfo, poleBounds] = getCellInfo( tracksFinal, cellMeshAll, pixelSize);
[nTracks, totalCells, cellNum, tracksLength, tracksOrigin, tracksFilled, badMovie] = getTracksInfo( tracksFinal);

% set the savePath for the tracksFinal file
tfPath = [varPath 'Track Analysis Archive\Single Day\'];
if ~exist( tfPath, 'dir')
    mkdir( tfPath)
end

tfName = [strain  ' tracksFinal' extraName ' ' Date];

save( [tfPath tfName], 'tracksFinal', 'cellRecord', 'cellMeshAll',...
    'cellInfo', 'poleBounds', 'nTracks', 'totalCells', 'pixelSize', 'timeStep', ...
    'varPath', 'savePath', 'folderName', 'fileName', 'expDate', 'meshPath', ...
    'tfPath', 'tfName', 'Date', 'strain', 'extraName')

fprintf( '~~~ Step 3: ''%s'' saved to **Track Analysis Archive\\Single Day** ~~~\n\n', tfName)


%% Step 4. Diffusion analysis (calculate quantities)

% remove original uTrack output fields in tracksFinal to save space
tracksFinal = rmfield( tracksFinal, {'tracksFeatIndxCG' 'tracksCoordAmpCG' 'seqOfEvents'});

% Calculate MSD & Diff & alpha & Fitting
diffAnalysis


%% Step 5. Convert the tracksFinal structure into Matrix, for plotting

matrixPath = [ varPath '\Matrix\Single Day\'];
if ~exist( matrixPath, 'dir')
    mkdir( matrixPath)
end

[ tracksLNorm, tracksxNorm, tracksLNorm4, tracksxNorm4, tracksDiff] = dataToVectors( tracksFinal);
% tracksDiff: diff (um^2/s), locErr (nm), alpha, dalpha (um^2/s)

% cell subregion constraints (cell pole or middle part)
bound = poleBounds( cellNum); % cell pole bound for each track
tracksMid = false( nTracks, 2); % based on [first spot, whole track]

LNorm = tracksLNorm; % [first spot, end spot, minLNorm, maxLNorm]

tracksMid(:, 1) = LNorm(:,1) > bound & LNorm(:,1) < 1-bound; % first spot not at pole
tracksMid(:, 2) = LNorm(:,3) > bound & LNorm(:,4) < 1-bound; % whole track not at pole

tracksMid4 = tracksLNorm4 > bound & tracksLNorm4 < 1-bound;

    % get xNorm 12 frame, added 12/26/23
    longT = tracksLength >= 12;
    tf = tracksFinal( longT);

    tracksxNorm12 = nan( length( tf), 12);
    tracksLNorm12 = nan( length( tf), 12);

    for i = 1: length( tf)
        spotNorm = tf(i).spotPosNorm; % [LNorm, xNorm]
        tracksLNorm12(i,:) = spotNorm(1:12, 1);
        tracksxNorm12(i,:) = spotNorm(1:12, 2);
    end

    LNorm = tracksLNorm12;    bound = poleBounds( cellNum( longT));
    inCell = tracksxNorm12 >= -1 & tracksxNorm12 <= 1;
    tracksMid12 = LNorm > bound & LNorm < 1-bound & inCell; % first spot subregion constraint
    
    
matrixName = ['Matrix ' strain extraName ' ' Date];
save( [matrixPath matrixName], 'cellRecord', 'cellMeshAll',...
    'cellInfo', 'poleBounds', 'nTracks', 'totalCells', 'cellNum', ...
    'tracksLength', 'tracksOrigin', 'tracksFilled', 'badMovie',...
    'pixelSize', 'timeStep', 'fitRange', 'maxTau', 'EnsMSD', 'EnsTAMSD',...
    'tracksLNorm', 'tracksxNorm', 'tracksLNorm4', 'tracksxNorm4', 'tracksDiff',...
    'tracksLNorm12', 'tracksxNorm12', 'tracksMid12',... % added 10/26/2023
    'tracksMid', 'tracksMid4', 'folderName', 'Date', 'strain', 'extraName',... 
    'varPath', 'savePath', 'matrixPath', 'matrixName');

fprintf( '\n~~~ Step 4: ''%s'' saved to ''varPath\\Matrix\\Single Day\\''~~~\n\n', matrixName);


%% Step 7. Display Analysis Result Summary

cellShifted = length( unique( cellNum( ~badMovie)));
cellFilled = length( unique( cellNum( tracksFilled)));

fprintf( ['%s Analysis Result: \n\n     %d total cells with %d tracks, %d shifted cells with %d tracks,\n'...
    '     There are %d filled cells (shifted) with %d tracks,\n     xNorm (1st spot) has %d tracks not in cap region, in total %d tracks \n'...
    '     xNorm (first 4 spots) has %d tracks not in cap region, in total %d spots\n\n'],...
    folderName, totalCells, nTracks, cellShifted, sum( ~badMovie),...
    cellFilled, sum( tracksFilled), sum( tracksMid(:,1)), sum( ~isnan( tracksxNorm(:,1))),...
    sum( min( tracksMid4, [], 2)), sum( sum( ~isnan( tracksxNorm4( min( tracksMid4, [], 2),:)))))
    
% disp( cellRecord)


