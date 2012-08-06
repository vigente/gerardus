function idx = findpeaks3(im, dims, dmin, pkmin)

% check arguments
error(nargchk(1, 4, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

if length(dims) ~= length(dmin)
    error('If minimum distances are provided between peaks, one distance must be provided per dimension in DIMS')
end

% defaults
if (nargin < 2)
    dims = true(1, ndims(im));
else % convert to boolean format
    aux = false(1, ndims(im));
    aux(dims) = true;
    dims = aux;
end
if (nargin < 3) || isempty(dmin)
    dmin = ones(1, ndims(im));
else
    aux = zeros(1, ndims(im));
    aux(dims) = dmin;
    dmin = aux;
end
if (nargin < 4) || isempty(pkmin)
    pkmin = 0;
end

% init peaks locations
idx = true(size(im));

% find peaks along each dimension
for D = 1:ndims(im)
    
    %% find peaks
    
    % if we don't want this dimension to have an effect in the definition
    % of a peak, we skip it
    if (~dims(D))
        im = shiftdim(im, 1);
        idx = shiftdim(idx, 1);
        continue
    end
    
    % we are going to be shifting the dimensions of im, so we can always
    % work along rows
    
    % first we sweep rows in the ascending order
    aux = sign(diff(im, 1));
    
    % assume that the background beyond the image is always 0
    aux = cat(1, aux, -sign(im(end, :, :)));
    
    % second we sweep rows in the descending order
    aux2 = sign(diff(im(end:-1:1, :, :), 1));
    aux2 = cat(1, -sign(im(1, :, :)), aux2(end:-1:1, :, :));
    
    % to mirror Matlab's findpeaks(), we consider peaks only the end of
    % increasing slopes. So if we have a raising slope, a flat, and then a
    % falling slope, we only consider one peak (the end of the raising
    % slope)
    aux = (aux<=0) & (aux2==-1);
    idx = idx & aux;
    
    %% remove peaks that are too small
    idx(im < pkmin) = false;
    
    %% remove peaks that are too close to other peaks, column by column, 
    %% and slice by slice
    
    for C = 1:size(im, 2)
        for S = 1:size(im, 3)
            aux = removeNearbyPeaks(im(idx(:, C, S), C, S), find(idx(:, C, S)), dmin(D));
            idx(:, C, S) = false;
            idx(aux, C, S) = true;
        end
    end
    
    %% prepare for next step
    
    % shift dimensions of im so that we can work always by rows
    im = shiftdim(im, 1);
    idx = shiftdim(idx, 1);
    
end
end

% auxiliary function to remove peaks that are too close to other peaks
function loc = removeNearbyPeaks(val, loc, dmin)

if (isempty(val) || (dmin == 1))
    return
end

% sort peaks according to value in descending order
[~, idx] = sort(val, 'descend');
loc = loc(idx);

% vector to note which peaks need to be deleted
todelete = false(size(loc));

% loop every peak
for I = 1:length(loc)
    if ~todelete(I)
        
        % we have found a peak that doesn't have to be removed. 
        
        % compute the distance of the other peaks to this one
        d = abs(loc(I) - loc);
        
        % peaks that are too close are tagged to be removed, except the
        % current peak itself
        todelete = todelete | (d < dmin);
        todelete(I) = 0;
        
    end
end

% remove peaks
loc(todelete) = [];

end
