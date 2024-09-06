function savePath = spotNorm_yh( meshPath, meshName, varPath, areaRatioThresh, goodCellFlag, imgSaveFlag)
%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
Last edited date: 8/9/2023

Description: spotNorm_yh.m is a script that combine oufti mesh files &
uTrack tracking outputs and align them together. It will then calculated
the normalized position of spots inside cells as [xNorm, lNorm]

~~~ We denote the short and long axes as x & l, respectively ~~~

1) load oufti mesh & uTrack output data
2) assign each track to the corresponding cell
3) make overlay plot of cell image & outline & tracks
4) shift back the shifted image by aligning the cell & tracks
5) find the normalized position of each track signal within the cell
6) save data into tracksFinal_Variables.mat file

======Input=======
    meshPath & methName
    varPath, abRatioThresh, goodCellFlag, imgSaveFlag
    from Oufti: BF00xm.mat {'cellList, cellListN'}
    from Utrack: Channel_1_tracking_result.mat {'tracksFinal'}

======Output=======
tracksFinal structure with new fields: (nTracks*1)
    - tracksCoordXY
    - cellNumber
    - spotPosInCell (optional)
    - spotPosNorm
    - tracksFinal_InIsolatedCell

=====Functions=====
    getExpInfo
    setCellMesh
    markIsoCells
    getBound
    getNormPos
    findPerpFoot
    findIntersect

Compare to oufti functions:
    projectToMesh.m
    getextradata.m

%}

%% Step 1: Load Oufti mesh & uTrack output data
global checkFlag img

    shiftFlag = true; % apply autoShift by default; false: don't apply autoShift
    checkFlag = false; % this flag is for checking the norm pos by plotting
    
    % get strain & Date information from the folder
    [folderName, fileName, Num, expDate, strain, extraName, imgName] = getExpInfo( meshPath, meshName);
    
    % set the savepath for spotNorm results: _Variables, including updated trackFinal
    savePath = [ varPath 'Variables\' strain '\' expDate extraName '\']; 
    if ~exist( savePath, 'dir') % create one if doesn't exist
        mkdir( savePath)
    end   
        
    % Load oufti output file (cell mesh info)
    load([meshPath meshName], 'cellList', 'cellListN');

    % load the tracking result file from uTrack
    trackFolder = dir([meshPath '*rack*' Num]);    % List the folder 'tracking001'
    tracksPath = '\TrackingPackage\tracks\';    tracksName = 'Channel_1_tracking_result.mat';
    load( [ trackFolder.folder '\' trackFolder.name tracksPath tracksName], 'tracksFinal')
     
    
%% Step 2: Assign cell numbers to tracks (tf.cellNumber = tracksCell)
    
    % Set up the cell mesh information (contains everything about cell mesh)
    cellMesh = setCellMesh( cellList, cellListN);    
    
    % Mark the cells as isolated good (1), or not isolated bad (0)
    cellIsolated = markIsoCells( goodCellFlag, cellListN, cellMesh, meshPath, imgName);
    
    nTracks = size( tracksFinal,1);
    tracksCell = nan( nTracks, 1); % cell number for each tracks
    
    for i = 1: nTracks
        
        nFrames = size( tracksFinal(i).tracksFeatIndxCG, 2); % number of frames in this track
        
        % save track X & Y coordinates (unit: pixel)
        Coords = tracksFinal(i).tracksCoordAmpCG;
        traj = nan( nFrames,2);
        traj(:,1) = Coords( 1: 8: end); % x coordinates
        traj(:,2) = Coords( 2: 8: end); % y coordinates        
        tracksFinal(i).tracksCoordXY = traj; % x & y coordiate in the unit of "pixels"
        
        % find which cell this track belongs to, If most of the spots in
        % the track are in a cell, then assign the cellNumber to this track
        for Cell = 1 : cellListN            
            if ~isempty( cellMesh( Cell).mesh)
                
                meshOut = cellMesh( Cell).meshOut;                
                % compare the track position & mesh to find which cell it belongs
                inFlag = inpolygon( traj(:,1), traj(:,2), meshOut(:,1), meshOut(:,2)); % return 0 or 1 in a vector
                if mode( inFlag) == 1 % most frequent values in the array
                    tracksCell(i) = Cell; % assign this track to this cell number
                    break
                end                
            end
        end        
        tracksFinal(i).cellNumber = tracksCell(i); % if not in cell: NaN
    end    
    goodTracks = ismember( tracksCell, find( cellIsolated)); % tracks inside selected isolated cells
    
    
%% Step 3: Plot & save overlay of cell image & outline & original tracks 
    
    % load greyscale tiff images
    img = imread([meshPath imgName]); imshow(img, []); set( gcf, 'Position', [1100 200 700 700]);  hold on

    tracksGood = cell2mat( { tracksFinal( goodTracks).tracksCoordXY}'); % inside cells
    tracksBad = cell2mat( { tracksFinal( ~goodTracks).tracksCoordXY}'); % outside cells
    
    % plot tracks
    if ~isempty( tracksGood)
        plot( tracksGood(:,1), tracksGood(:,2) , 'c.'); % cyan for tracks inside cells
    end
    if ~isempty( tracksBad)
        plot( tracksBad(:,1), tracksBad(:,2) , 'r.'); % red for bad tracks inside cells
    end
    
    % plot cell mesh   
    for Cell = 1 : cellListN       
        meshOut = cellMesh( Cell).meshOut;             
        plot( meshOut(:,1), meshOut(:,2), 'b-', 'LineWidth', 0.6)
        text( ( meshOut(10, 1)- 5), ( meshOut(10, 2)+ 5), num2str(Cell),...
            'Color', 'r', 'FontSize', 16, 'FontWeight', 'bold')
    end    
    text(5,13, fileName, 'Color', 'k', 'FontSize', 16)   % as '230712-SK187 Suc'
    if imgSaveFlag
        exportgraphics( gca, [meshPath 'Mesh & Tracks Image\' imgName(1:end-4) '_TracksMesh.png'], 'Resolution', 150)
        disp( ['       ~~~ ' imgName(1:end-4) ' Overlay Image Saved ~~~' ]);
    end
    
%% Step 4: Auto Shift Back (find shift value from filled cells)
    
    badMovie = false;   shift = [nan nan];
    
    if ~shiftFlag % don't run auto shift back
        cellFilled = false( cellListN, 1); % won't include in xNorm analysis
        badMovie = true;
        saveName = [fileName '_Variables'];
        fprintf( '~~~~ Shift Back NOT selected, Bad Movie! ~~~~ \n ')
    else
    % apply shift to all tracks & re-assign cell numbers to the tracks
    
        saveName = [fileName '_SB_Variables'];
        areaRatio = nan( cellListN, 1); % area ratio of tracks & cell outline
        cellShift = nan( cellListN, 2); % shift value in [x, y] for each cell

        tfGood = tracksFinal( goodTracks); % tracks inside cells (good)
        cellNum = [ tfGood.cellNumber]'; % cell number of each good track
        goodCells = unique( cellNum)'; % all goodCells are not empty and have tracks inside
        
        % Find shift value for each cell
        for Cell = goodCells

            thisCell = cellNum == Cell; % index of tracks in this cell
                                  
            % stack spot coordinates of all tracks in this cell
            traj = cell2mat( { tfGood( thisCell).tracksCoordXY}'); 

            % find the area covered by the tracks in this cell
            
            bound = getBound( traj(:,1), traj(:,2)); % get the boundary of the tracks
            
            spotsArea = polyarea( bound(:,1), bound(:,2));
            % find the area of the cell outline
            meshOut = cellMesh( Cell).meshOut;
            meshArea = polyarea( meshOut(:,1), meshOut(:,2)); 

            areaRatio( Cell) = spotsArea/ meshArea;

            % Compare the location of the mask & the cell mesh (edge + center) 
            cellShift( Cell, :) = mean( [mean(meshOut)-mean(bound); min(meshOut)-min(bound); max(meshOut)-max(bound)]); 
            
            % plot the boundary of tracks in Cell
%             plot( bound(:,1), bound(:,2), 'w', 'LineWidth', 2)
%             fprintf( '    Cell #%d:  spotsMask/meshSize = %.2f \n', Cell, spotsArea/meshArea );
        end
                
        % find 'shift' value from mostly-filled cells
        cellFilled = areaRatio > areaRatioThresh;  % mostly-filled cells        
        if sum( cellFilled) < 3
            fprintf( '~~~~ %s Lacks enough good cells, therefore lower the areaRatio threshold ~~~~ \n ', fileName)
            cellFilled = areaRatio > areaRatioThresh - 0.1;
        end
        if sum( cellFilled) > 2 % if <3 mostly filled cells, shift = NaN
            shift = mean( cellShift( cellFilled, :), 1); 
        end
        
        % shiftBack and examine the result
        cellColor = { "#EDB120", 'r'}; % two color for bad (yellow) & good (red) cells

        if isnan( sum( shift)) % not enough mostly-filled cells
            cellFilled = false( cellListN, 1);
            badMovie = true;
            fprintf( '~~~~ %s Lacks mostly-filled cells, Bad Movie! ~~~~ \n ', fileName)
        else
        % shift all tracks & replot the overlay
        
            fprintf( '    Mean Shift = [%.2f, %.2f],  %d/%d filled cells\n',...
                shift(1), shift(2), sum( cellFilled), cellListN)

            % plot shiftBack tracks & cell outline
            imshow( img, []); hold on

            for Cell = 1 : cellListN

                % plot shifted tracks
                thisCell = cellNum == Cell;
                if sum( thisCell) > 0 % make sure this cell is good            
                    traj = cell2mat( { tfGood( thisCell).tracksCoordXY}') + shift; % coordinates of all spots in tracks in that cell
                    plot( traj(:,1), traj(:,2) , 'c.') 
                    hold on
                end

                % plot cell outline & number (red for filled cells, yellow for non-filled cells)
                meshOut = cellMesh( Cell).meshOut;
                plot( meshOut(:,1), meshOut(:,2), 'b-', 'LineWidth', 0.6)
                text( ( meshOut(10,1)-4), ( meshOut(10,2)+4), num2str(Cell),...
                    'Color', cellColor{ cellFilled( Cell)+1}, 'FontSize', 16, 'FontWeight', 'bold')
            end
            hold off  

            text(5,13, [fileName '    auto ShiftBack'], 'Color', 'k', 'FontSize', 15)
            if imgSaveFlag
                exportgraphics( gca, [meshPath 'Mesh & Tracks Image\auto shiftBack\' imgName(1:end-4) ' autoShiftBack.png'], 'Resolution', 200)
            end
            close
            
            % shift back the coordiate to all tracks
            for j = 1: nTracks
                tracksFinal(j).tracksCoordXY = tracksFinal(j).tracksCoordXY + shift;
            end
        end
    end
    
    % re-assign the tracks to cells after the shiftBack
    tracksCell = nan( nTracks, 1); % cell number for each tracks
    for i = 1: nTracks
        traj = tracksFinal(i).tracksCoordXY; % x & y coordiate in the unit of "pixels"
        % find which cell this track belongs to
        for Cell = goodCells
            meshOut = cellMesh( Cell).meshOut;
            % compare the track position & mesh to find which cell it belongs
            inFlag = inpolygon( traj(:,1), traj(:,2), meshOut(:,1), meshOut(:,2)); % return 0 or 1 in a vector
            if mode( inFlag) == 1 % most frequent values in the array
                tracksCell(i) = Cell; % assign this track to this cell number
                break
            end
        end
        tracksFinal(i).cellNumber = tracksCell(i); % if not in cell: NaN
    end
    
    goodTracks = ismember( tracksCell, find( cellIsolated)); % tracks inside selected isolated cells
    tfGood = tracksFinal( goodTracks);  % tracks inside cells (good)
    cellNum = [ tfGood.cellNumber]';    % cell number of each good track
    goodCells = unique( cellNum)';      % all goodCells are not empty and have tracks inside
    
%% Step 5: Find the normalized positions of spots within the cell and save to tfGood  

    % find track position inside cells, add field 'spotPosNorm'
    %      & field 'spotPosInCell' (optional)
    tfGood = getNormPos( tfGood, cellMesh); 
    
    tracksFinal_InIsolatedCell = tfGood;
    
%% Step 6: Save Cells & Tracks Information of the Movie

    save( [ savePath saveName],...
         'tracksFinal', 'tracksFinal_InIsolatedCell', 'cellListN', 'cellList', 'cellMesh', ...
         'cellIsolated', 'areaRatio', 'areaRatioThresh', 'cellShift', 'cellFilled', 'shift', ...
         'goodTracks', 'cellNum', 'goodCells', 'img', 'badMovie',  ...
         'meshPath', 'meshName', 'goodCellFlag', 'varPath', 'savePath', 'saveName', ...
         'folderName', 'fileName', 'expDate', 'strain', 'extraName', 'imgName')
     
    fprintf( '    ''%s'' saved under ''varPath\\%s\\'' \n\n', saveName, [strain '\' expDate extraName]);
    
end


%% Functions

    function [folderName, fileName, Num, expDate, strain, extraName, imgName] = getExpInfo( meshPath, meshName)
    % Extract experimental info from meshPath & meshName
    % Naming matters! please follow the format rule example for the folder:
    % '230712-SK187' or '230712-SK187 Suc'  

        Num = meshName(3:5);                                % follow rules of naming: BF001, BF021
        tmp = split(meshPath, '\');
        folderName = tmp{end-1};                            % FolderName: '230712-SK187 Suc'
        fileName = [folderName ' ' meshName(1:end-5)];      % folderName + imgName as '230712-SK187 Suc BF001'

        tmp2 = split(folderName, ["-" " "]);                % split '230712-SK187 Suc'
        expDate = tmp2{1};                                  % date of the experiment '230712'
        strain = tmp2{2};                                   % strain number 'SK187'
        extraName = folderName( length( [tmp2{1:2}])+2: end); % anything after ' Suc'

        imgName = [meshName(1:end-5) '.tif'];                % m.mat --> .tif
        disp([' -- We are looking at ' meshName ' now ~~~']);
    end

    function cellMesh = setCellMesh( cellList, cellListN)
    % store mesh info of all cells into a structure called cellMesh
    % mesh, meshOut, meshMid, gridLen, gridLenCum

        cellMesh( cellListN, 1) = struct;
        for Cell = 1: cellListN
            
            % cell mesh: (n*4 single) - (x1,y1,x2,y2) - (left,right)
            mesh = double( cellList.meshData{ 1}{ Cell}.mesh); % convert to double

            if length( mesh) > 1 % sometimes the cell mesh is empty 

                meshMid = [ mean( mesh(:,[1 3]), 2), mean( mesh(:,[2 4]), 2)]; % midline along the long axis
                % save the cell mesh info for later use
                cellMesh( Cell).mesh = mesh;
                cellMesh( Cell).meshOut = [ mesh(:, 1:2); flipud( mesh(:, 3:4))]; % reshape the mesh matrix to form a circle [2n, 2]
                cellMesh( Cell).meshMid = meshMid;
                cellMesh( Cell).gridLen = vecnorm( diff( meshMid), 2, 2); % length of each grid (L direction)
                cellMesh( Cell).gridLenCum = cumsum( cellMesh( Cell).gridLen); % cumulative length of each grid (L direction)
            else
                % Sometimes the Oufti went wrong and some cells has no information
                fprintf( '~~~~~~ Cell #%d is having problem with its mesh! ~~~~~~', Cell)
            end
        end
    end

    function cellIsolated = markIsoCells( goodCellFlag, cellListN, cellMesh, meshPath, imgName)
    % Mark the cells as isolated good (1), or non-isolated bad (0)

        cellIsolated = true( cellListN, 1); % logical vector to store good isolated cells
        if logical( goodCellFlag) % flag set in autoSportNorm: autorun or manually choose good cells
            disp( '    You have assumed all cells are good, so I will run by myself');
        else 
        % plot cell image & mesh        
            img = imread( [ meshPath imgName]); imshow( img, []); 
            set( gcf, 'Position', [1100 200 700 700]);  hold on  
            for Cell = 1: cellListN
                meshOut = cellMesh( Cell).meshOut;
                plot( meshOut(:,1), meshOut(:,2), 'b-', 'LineWidth', 1)
                text( ( meshOut(10, 1)- 5), ( meshOut(10, 2)+ 5), num2str( Cell),...
                    'Color', 'r', 'FontSize', 16, 'FontWeight', 'bold')
            end
            text( 5,13, imgName, 'Color', 'k', 'FontSize', 20)
            
            flag = logical( input('    Hi there, are all cells isolated? Yes(1), No(0): '));
            if ~flag
                badCells = input('    Tell me the bad cell numbers (input as [x,x,x]): ');
                cellIsolated( badCells) = false;
            end
            close
        end
    end
    
    function bound = getBound( x, y)
    % returns a smoothed boundary of the input 2D coordinates (x,y),
    % created by Yu-Huan Wang @7/12/2023 
    
        k = boundary( x, y, 0.2); % Boundary of a set of points in 2-D        
        boundRaw = [ x(k), y(k)];

        % smooth the boundary
        boundN = [ boundRaw( end-1,:); boundRaw; boundRaw( 2,:)];
        tmp = movmean( boundN, 3);
        bound = tmp( 2: end-1, :);

        plot( x(k), y(k), 'r', 'LineWidth', 2)    
        plot( bound(:,1), bound(:,2), 'w', 'LineWidth', 1.5)
    end
        
    function tfGood = getNormPos( tfGood, cellMesh)
    % for each point in a track, calculate the distance to all semgment
    % edges (minor axis) with side info +/- to find which segment it is
    % located at. To find the normalized position inside cell, a parallel
    % line (with mean slope of the neighboring two segment edges) passing
    % through the point is created to find the intersection point to the
    % cell mid line & outline. The length of the cell width at that point
    % is used for normalization
    % created by Yu-Huan Wang @7/23/2023
        global checkFlag img
        
        for i = 1: length( tfGood) % only apply spotNorm to tracks inside selected cells

            thisCell = tfGood(i).cellNumber;
            traj = tfGood(i).tracksCoordXY;
            nFrames = length( traj);

            mesh = cellMesh( thisCell).mesh; % x1, y1 (left side), x2, y2 (right side)
            meshOut = cellMesh( thisCell).meshOut;
            meshMid = cellMesh( thisCell).meshMid;
            lenCum = [0; cellMesh( thisCell).gridLenCum]; % to match the mesh index
            
%             spotPosInCell = nan( nFrames, 2); % absolute position inside cell
            spotPosNorm = nan( nFrames, 2); % normalized position inside cell
            
            for frame = 1: nFrames
                
                pt = traj( frame, :); % point in a trajectory
                if ~inpolygon( pt(1), pt(2), meshOut(:,1), meshOut(:,2)) % pt outside cell
                    continue % position remains nan
                end
                
                % find the distance of this point to cell segments edges (nan for 1st & end element)
                [~, dist] = findPerpFoot( pt, mesh(:, 3:4), mesh(:, 1:2)); % -: counter-clock, +: clock, should be - to +
                bra = find( abs( diff( sign( dist))) == 2); % find the segment idx where distance changes sign
                % dist(bra) & dist(ket) are the distance of the point to its neighboring segment edges
                
                % Parallel Methods
                if isempty( bra) % not sandwiched between 2 segment edges, must be in two caps then
                    if sum( dist< 0) == 0       % all positive, it's in the 1st segment
                        bra = 1;
                    elseif sum( dist> 0) == 0   % all negative, it's in the last segment
                        bra = length( dist)- 1;
                    end
                end
                ket = bra + 1;
                
                % use parallel line of segment edges to find intersection
                % with midline (LNorm) & outline (xNorm), modified @7/23/2023 by YHW
                a = mesh( [bra ket],:);
                
                % find the mean slope of the neighboring two segments
                slopes = (a(:,4)-a(:,2))./ (a(:,3)-a(:,1));
                if sum( sign( slopes)) == 0 % two slopes have opposite sign
                    k = abs( diff( slopes))* sign( sum( slopes));
                else
                    k = mean( slopes, 'omitnan');
                end
                
                [D, xPos] = findIntersect( pt, k, meshMid(ket,:), meshMid(bra,:));
                % find the intersection point of pt-D & the cell outline (cell width at lPos for this pt)
                if xPos > 0 % on the right side of the mid line
                    [IntPt, ~] = findIntersect( pt, k, mesh( bra, 3:4), mesh( ket, 3:4));
                else        % on the left side of the mid line
                    [IntPt, ~] = findIntersect( pt, k, mesh( bra, 1:2), mesh( ket, 1:2));
                end
                xNorm = xPos/ norm( IntPt-D); % normalized by the width at that point
                lPos = lenCum( bra) + norm( D- meshMid( bra,:)); % real L value by portion
                lNorm = lPos/ lenCum(end); % normalized by the total length of the major axis
                
                if xNorm > 1 || xNorm < -1
                    warning( 'xNorm outside the [-1 1] region')
                end
                if lNorm < 0 || lNorm > 1
                    % if slope of the neighboring two segments are almost vertical,
                    % then the mean slope would be ~0, which is not right
                    warning( 'lNorm outside the [0 1] region')
                end
                
%                 xNorm = max( -1, min( 1, xNorm)); % set range at [-1 1]
%                 lNorm = max( 0, min( 1, lNorm)); % set range at [0 1]
                
%                 spotPosInCell( frame,:) = [lPos, xPos];
                spotPosNorm( frame,:) = [lNorm, xNorm];
                
                %% Check by plotting points & lines & intersections
                
                condCheck = cellMesh( thisCell).gridLenCum(end) < 20 && abs(lNorm-0.5)<0.2;
                % frame == 1; % && abs(lNorm-0.5)>0.45 % cell pole
                if condCheck && checkFlag
                    imshow( img, []); set( gcf, 'Position', [1100 200 700 700]); hold on
                    plot( meshOut(:,1), meshOut(:,2), 'b', 'LineWidth', 1)
                    plot( meshMid(:,1), meshMid(:,2), 'b', 'LineWidth', 1)
                    for n = 1: length( mesh)-1
                        plot( mesh( n, [1 3]), mesh( n, [2 4]), 'b', 'LineWidth', 1)
                    end
                    roi = 10;
                    xlim( [pt(1)-roi pt(1)+roi])
                    ylim( [pt(2)-roi pt(2)+roi])
                    
                    % Lpos
                    plot( meshMid(:,1), meshMid(:,2), 'k', 'LineWidth', 10) % total length
                    plot( [ meshMid( 1:bra, 1); D(1)], [ meshMid( 1:bra, 2); D(2)],...
                        'r', 'LineWidth', 3)
                    plot( [IntPt(1) D(1)], [IntPt(2) D(2)], 'm', 'LineWidth', 3)
                    scatter( D(1), D(2), 80, 'c', 'filled') % perpendicular foot
                    scatter( pt(1), pt(2), 200, 'w', 'filled') % signal point
                    
                    % xpos
                    plot( [IntPt(1) D(1)], [IntPt(2) D(2)], 'm', 'LineWidth', 10)
                    plot( [D(1) pt(1)], [D(2) pt(2)], 'w', 'LineWidth', 3)
                    scatter( IntPt(1), IntPt(2), 80, 'c', 'filled') % signal point
                    scatter( D(1), D(2), 80, 'c', 'filled') % perpendicular foot
                    scatter( pt(1), pt(2), 200, 'w', 'filled') % signal point
                    
                    figure( gcf)
                    fprintf( 'track %d, frame %d, LNorm = %.2f, xNorm = %.2f\n', i, frame, lNorm, xNorm)
                    gcf;
                end
            end

%             tfGood(i).spotPosInCell = spotPosInCell; % absolute position inside cell
            tfGood(i).spotPosNorm = spotPosNorm; % normalized position inside cell
        end
%         disp( '        spotNorm Finished ~~~');
    close
    end
    
    function [D, dist] = findPerpFoot( pt, B, C)
    % this function returns the coordinate of the perpendicular foot D so
    % that AD perpendicular to BC (everything in 2D) and the distance of pt
    % to line BC
    % in image, top left [0,0] & bottom right [N,N], the y axis is inverted
    % dist > 0 if pt is on the left side of line BC (pt, B, C: clockwise)
    % dist < 0 if pt is on the right side of line BC (pt, B, C: counterclockwise)

        AB = B - pt; % vector AB
        BC = C - B;  % vector BC
        
        % cross product, right hand rule
        area = AB(:,1).*BC(:,2) - AB(:,2).*BC(:,1); % AB cross BC 
        side = sign( area); % -1: left side, +1: right side

        normVec = [ BC(:,2) -BC(:,1)]; % normal vector of BC in 2D 
        unitNormVec = normVec./ vecnorm( normVec, 2, 2); % unit normal vector of BC
        AD = dot( unitNormVec, AB, 2).* unitNormVec; % AD is perpendicular to BC, dot product
        D = pt + AD; % D point of intersection, Perpendicular Foot
        dist = side.* vecnorm( AD, 2, 2);
    end

    function [IntPt, dist] = findIntersect( pt, k1, B, C)
        
        IntPt = nan( 1, 2);
        b1 = pt(2) - k1*pt(1); % pass through pt
        
        % kx1 + b = y1
        % kx2 + b = y2    
        % k = (y2-y1)/ (x2-x1);  b = y1 - kx1;
        k2 = ( C(2)- B(2))/ (C(1)- B(1));
        b2 = C(2) - k2* C(1);
    
        % k1x + b1 = y
        % k2x + b2 = y
        % x = - (b1-b2)/ (k1-k2)
        % y = k1x + b1 = k2x + b2
        IntPt(1) = -( b1-b2)/ (k1-k2);
        IntPt(2) = k2* IntPt(1) + b2;
        
        if sum( isnan(IntPt)) || sum( isinf(IntPt))
            disp('~~~ Intersection has problem! ~~~')
            IntPt = [nan nan];
        end
        
        % cross product, right hand rule opposite, +1: clockwise, -1: counter        
        AB = B - pt; % vector AB
        BC = C - B;  % vector BC
        side = sign( AB(:,1).*BC(:,2) - AB(:,2).*BC(:,1)); 
        dist = side.* vecnorm( IntPt-pt, 2, 2);
        
    end
    
