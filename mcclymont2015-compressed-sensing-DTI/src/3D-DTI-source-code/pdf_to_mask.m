function mask = pdf_to_mask(PDF, max_values, fully_sampled_region)

if nargin < 3
    fully_sampled_region = PDF * 0;
end

if nargin < 2% do we allow values of more than 1 on the mask?
    max_values = ceil(PDF);
end


if sum(fully_sampled_region(:) .* max_values(:)) > sum(PDF(:))
    disp('The fully sampled region needs more samples than your PDF allows!')
end

mask = zeros(size(PDF));
mask(fully_sampled_region == 1) = max_values(fully_sampled_region == 1);

pdf_sum = round(sum(PDF(:)));

while sum(mask(:)) < pdf_sum
    % rounding errors
    
    n_remaining = pdf_sum - sum(mask(:));
    
    [~,pos] = histc(rand(1,n_remaining), [0 ; cumsum(PDF(:))] ./ sum(PDF(:))) ;
    
    if sum(pos(:)) == 0
        break
    end
    
    for p = 1:length(pos);
        mask(pos(p)) = mask(pos(p)) + 1;
    end
    
    %mask(pos) = mask(pos) + 1;
    
    mask = min(mask, max_values);
    
    PDF(mask == max_values) = 0;

end

