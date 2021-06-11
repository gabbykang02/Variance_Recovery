function [idList] = isolateVisibleSomaMod(locs,thresh)

vMid         = 150/2;                                   % Get midpoint of volume (plane of imaging)
nNeur        = 338;                                       % Get number of neural somas in the volume
zLocs        = locs(1:nNeur, end);                                    % Get the z (axial) locations of the neurons
psf = [3.7617712,3.7300000,16.229452];
psf = [3.4506888,3.5908113,14.195392];
thresh       = 5.9000 + psf(end) / 2;   
idList       = find(abs(zLocs - vMid)<thresh);                           % Find all the somas no more than 'thresh' away from the scan plane


end

