function h = image_histogram( im, mask, make_plot )

% h = image_histogram( im, mask )
%
% calculate the histogram of the given image, optionally restricted on the mask
% h(1) is the number of pixels of value 0, etc.
% F. Nedelec, Dec. 2007

%%compatibility with tiffread:
if ( isfield(im,'data') ) 
    im = double( im.data ); 
end

if nargin < 3
    make_plot = 0;
end

%%
if nargin < 2 || isempty(mask)
    
    max_val = ceil(max(reshape( im, numel(im), 1 )));
    
    h = histc( reshape( im, numel(im), 1 ), 0:max_val);
    %h = sum( histc( im, 0:max_val ), 2);
        
else

    if any( size(im) ~= size(mask) )
        error('Image and mask must be have the same size');
    end
    if min(reshape(mask, numel(mask), 1 )) < 0 
        error('Values of the mask should be non-negative');
    end
    if max(reshape(mask, numel(mask), 1 )) > 1
        error('Values of the mask should be lower or equal to 1');
    end

    imm = im .* mask - ( 1-mask );
    val = reshape( imm, numel(imm), 1 );
    max_val = ceil(max(val));
    
    h = histc(val, 0:max_val);
    
end

%% Make a figure to display the histogram
if make_plot
 
    x = (0:size(h)-1)';
    figure('Name',inputname(1), 'Position', [100 150 800 300]);
    axes('Position', [0.05 0.1 0.9 0.8] );
    xlim([0 size(h,1)]);    
    plot(x, h, 'g.' );
    title('Histogram of pixel values');

end



end
