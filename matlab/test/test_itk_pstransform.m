% test_itk_pstransform.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test 1:
%% 2D grid with 4 landmarks

% source landmarks
x = [...
    .25 .25;...
    .25 .75;...
    .75 .25;...
    .75 .75...
    ];

% target landmarks
y = [...
    .25+.1 .25+.1;...
    .25-.2 .75+.1;...
    .75-.1 .25-.15;...
    .75+.15 .75+.2...
    ];

% grid points to be transformed
[gu, gv] = meshgrid(linspace(0, 1, 11), linspace(0, 1, 11));
xi = [gu(:) gv(:)];

% loop every available transform
for TRANSF = {'elastic', 'elasticr', 'tps', 'tpsr2', 'volume', 'bspline'}
    
    % apply transform to grid
    yi = itk_pstransform(TRANSF{:}, ...
        [x zeros(size(x, 1), 1)], ...
        [y zeros(size(y, 1), 1)], ...
        [xi zeros(size(xi, 1), 1)]);
    
    % plot points
    hold off
    plot(x(:, 1), x(:, 2), 'or')
    hold on
    plot(y(:, 1), y(:, 2), 'xg')
    plot(xi(:, 1), xi(:, 2), '.')
    plot(yi(:, 1), yi(:, 2), '.k')
    title(TRANSF{:})
    
    pause
    
end
