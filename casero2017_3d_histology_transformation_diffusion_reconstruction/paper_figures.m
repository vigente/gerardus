% paper_figures.m
%
% Script to generate the plots/figures used in the paper.
%
% Version: 0.3.2
% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2016 University of Oxford
% 
% University of Oxford means the Chancellor, Masters and Scholars of
% the University of Oxford, having an administrative office at
% Wellington Square, Oxford OX1 2JD, UK. 
%
% This file is part of Gerardus.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. The offer of this
% program under the terms of the License is subject to the License
% being interpreted in accordance with English Law and subject to any
% action against the University of Oxford being under the jurisdiction
% of the English Courts.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% laptop at home
BASEDIR = '/media/rcasero/mvldata/data';
PROJDIR = '/home/rcasero/Software/casero2017_3d_histology_transformation_diffusion_reconstruction';

% workstation at the office
BASEDIR = '/data2';
PROJDIR = '/home/orie1416/Software/casero2017_3d_histology_transformation_diffusion_reconstruction';

% documentation directories
DOCDIR = [PROJDIR filesep 'doc'];
FIGDIR = [DOCDIR filesep 'figures'];
SRCDIR = [PROJDIR filesep 'src'];
RESDIR = [SRCDIR filesep 'results'];

MOUSE = 'Q53';

IMPROCDIR = [BASEDIR '/Mouse/' MOUSE '/ImageProcessing'];
BFDATADIR = [BASEDIR '/Mouse/' MOUSE '/Blockface'];
BFSDATADIR = [IMPROCDIR '/BlockfaceStabilised'];
BFCDATADIR = [IMPROCDIR '/BlockfaceCorrected'];
BFPDATADIR = [IMPROCDIR '/BlockfacePreprocessed'];
% HISTO = 'HistologyTrichrome';
HISTO = 'HistologySiriusRed';
HISTOLO = [HISTO 'LoRes'];
HISTLODIR = [BASEDIR '/Mouse/' MOUSE filesep HISTOLO];
HISTDIR = [BASEDIR '/Mouse/' MOUSE filesep HISTO];
HISTLOREGDIR = [BASEDIR '/Mouse/' MOUSE filesep HISTOLO '_to_BlockfaceCorrectedRegistration'];
%HISTLO2BFREGSIMILARITYDIR = [IMPROCDIR filesep HISTOLO 'GrayToBlockfaceRegistrationSimilarity'];
HISTLO2BFREGRIGID = [HISTOLO '_to_BlockfaceRegistrationRigid'];
HISTLO2BFREGRIGIDDIR = [IMPROCDIR filesep HISTLO2BFREGRIGID];
INTRAHISTLORIGID = ['Intra_' HISTOLO];
INTRAHISTLOBSP = ['Intra_' HISTOLO '_Bsp'];
INTRANOBFHISTLORIGID = ['NoBlockfaceIntra_' HISTOLO];

HISTLOPDIR = [HISTLO2BFREGRIGIDDIR 'Preprocessed'];

HISTREC = [HISTO 'Reconstructed'];
HISTREC_LVFREE = [HISTREC '_LV_FreeWall'];
HISTREF_LVFREE = [HISTREC 'Refined_LV_FreeWall'];
HISTRECDIR = [IMPROCDIR filesep HISTREC];
HISTRECDIR_LVFREE = [HISTRECDIR '_LV_FreeWall'];
HISTREFDIR_LVFREE = [HISTRECDIR 'Refined_LV_FreeWall'];

% initial alignment of histology to blockface, common to all subsequent
% refinements
T_HISTLO2BF = ['t_' HISTOLO '_to_Blockface.mat'];
T_HISTLO2BF_FILE = [IMPROCDIR filesep T_HISTLO2BF];
T_INTRAHISTLO = ['t_Intra_' HISTOLO '.mat'];
T_NOBF_INTRAHISTLO = ['t_NoBlockfaceIntra_' HISTOLO '.mat'];
T_TOTHIST_LVFREE = ['t_Total_' HISTO '_LV_FreeWall.mat'];

% manually traced landmarks
PTS_HIST2BF = ['pts_' HISTOLO '_to_Blockface.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% repeat convolution of solution kernels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of repeat convolutions
for N = [10 40 80];
    
    h2 = cell(0);
    H2 = cell(0);
    
    % values of constant alpha
    for a = [0.50 0.49 0.45];
        
        % create the local neighborhood kernel
        h = [a 1-2*a a];
        
        % repeated convolution
        h2{end+1} = h;
        for I = 1:N; h2{end} = conv(h2{end}, h); end
        
        % in frequency domain
        H2{end+1} = abs(fft(h2{end}));
        
    end

    % convert from cell to matrix format
    h2 = cell2mat(h2');
    H2 = cell2mat(H2');
    
    % plot repeat convolution
    hold off
    N2 = length(h2);
    x2 = (1:N2)-(N2+1)/2;
    plot(x2, h2, 'LineWidth', 3)
    set(gca, 'FontSize', 18)
    xlabel('Relative slice index')
    ylabel('\Psi')
    axis([min(x2) max(x2) 0 max(h2(:))*1.1])
    title(['m = ' num2str(N)])
    if (N == 80)
        legend('\alpha = 0.50', '\alpha = 0.49', '\alpha = 0.45', 'Location', 'NorthEast')
    end
    
    saveas(gca, [FIGDIR filesep 'kernel-spatial-alpha-m-' num2str(N) '.tif'])
    
    % save as CSV for Elsevier's Interactive Plots
    csvfile = [FIGDIR filesep 'kernel-spatial-alpha-m-' num2str(N) '.csv'];
    fid = fopen(csvfile, 'w');
    fprintf(fid, '%s,%s,%s,%s\n', ...
        'Relative slice index', 'Ψ(α = 0.50)', 'Ψ(α = 0.49)', 'Ψ(α = 0.45)');
    fclose(fid);
    dlmwrite(csvfile, ...
        [x2', h2'], '-append', 'delimiter', ',','precision', '%.3E');
    system(['sed -i ''s/E+/E/g'' ' csvfile]);

    % plot frequency domain
    M2 = length(H2);
    theta2 = linspace(0, 2, M2+1);
    theta2(end) = [];
    plot(theta2, H2, 'LineWidth', 3)
    set(gca, 'FontSize', 18)
    xlabel('Angular frequency \Omega/ \pi')
    ylabel('|DFT(\Psi)|')
    axis([0 2 0 1.01])
%     if (N == 80)
%         legend('\alpha = 0.50', '\alpha = 0.49', '\alpha = 0.45', 'Location', 'Best')
%     end
    
    saveas(gca, [FIGDIR filesep 'kernel-freq-alpha-m-' num2str(N) '.tif'])

    % save as CSV for Elsevier's Interactive Plots
    csvfile = [FIGDIR filesep 'kernel-freq-alpha-m-' num2str(N) '.csv'];
    fid = fopen(csvfile, 'w');
    fprintf(fid, '%s,%s,%s,%s\n', ...
        'Angular frequency Ω/π', '|DFT_Ψ(α = 0.50)|', ...
        '|DFT_Ψ(α = 0.49)|', '|DFT_Ψ(α = 0.45)|');
    fclose(fid);
    dlmwrite(csvfile, ...
        [theta2', H2'], '-append', 'delimiter', ',','precision', '%.3E');
    system(['sed -i ''s/E+/E/g'' ' csvfile]);

end

% plot also the legend box so that we can cut it off with Gimp and add it
% to the paper (we cannot add the legend to the figures, because it's too
% big)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualization of transformation diffusion with vertical translation
%% (solution is sine wave)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coefficient
alpha = 0.49;

% number of slices
N = 100;

% place slices along true shape as sine wave
x = (0:(N-1))/N;
t = sin(2*pi*x);

% keep a copy of the true shape
t0 = t;

% random initialization of slice positions away from their perfect
% alignment
rng(10)
noise = rand(1, N);
noise = noise - mean(noise);
t = t + 2.4 * noise;

% plot slice positions
hold off
plot(0:N-1, t, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
hold on
plot(0:N-1, t0, ':k', 'LineWidth', 2)
axis([0 100 -3 2])
xlabel('slice index')
ylabel('slice position')

% save figure
saveas(gca, [FIGDIR filesep 'registration-diffusion-toy-translation-ground-truth.png'])

% save as CSV for Elsevier's Interactive Plots
csvfile = [FIGDIR filesep 'fig4a-registration-diffusion-toy-translation-ground-truth.csv'];
fid = fopen(csvfile, 'w');
fprintf(fid, '%s,%s,%s\n', ...
    'slice index', 'slice position (ground truth)', ...
    'slice position (noisy)');
fclose(fid);   
dlmwrite(csvfile, ...
    [(0:N-1)', t0', t'], '-append', 'delimiter', ',','precision', '%.3E');
system(['sed -i ''s/E+/E/g'' ' csvfile]);

% compute registration between slices
fprev = -diff(t);
fpos = -diff(t(end:-1:1));
fpos = fpos(end:-1:1);

% pad with values so that computing the slice update step is more simple.
% After padding, for the i-th slice we can do alpha*(fprev(i)+fpos(i))
fprev = [NaN fprev];
fpos = [fpos NaN];

% init accumulated transform to apply to each slice
ftot = zeros(1, N);

% NOTE: the optimal reconstruction is ITER = 42

% allocate space for reconstruction error
err = zeros(1, 7000);

%% reduce amount of registration noise, but don't get to optimal reconstruction
for I = 1:5

    % impose Neumann conditions at the extremes, i.e. 2*delta(0,1) and
    % 2*delta(N-1,N-2) respectively
    fprev(1) = fpos(1);
    fpos(end) = fprev(end);
    
    % compute the transform to apply to each slice in this iteration
    f = alpha * (fprev + fpos);
    
    % compose with the previous transforms to obtain a total transform
    ftot = ftot + f;
    
    % update neighbour transformation step (note that we ignore the element
    % for the ghost slices)
    fprev(2:end) = fprev(2:end) - f(2:end) + f(1:end-1);
    fpos(1:end-1) = fpos(1:end-1) - f(1:end-1) + f(2:end);
    
    % compute the error between the reconstruction and the true shape
    err(I) = norm(t + ftot - t0);
    
   
end

% plot slice positions
hold off
plot(0:N-1, t + ftot, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
hold on
plot(0:N-1, t0, ':k', 'LineWidth', 2)
axis([0 100 -3 2])
xlabel('slice index')
ylabel('slice position')

% save figure
saveas(gca, [FIGDIR filesep 'registration-diffusion-toy-translation-iter-5.png'])

% save as CSV for Elsevier's Interactive Plots
csvfile = [FIGDIR filesep 'fig4b-registration-diffusion-toy-translation-iter-5.csv'];
fid = fopen(csvfile, 'w');
fprintf(fid, '%s,%s,%s\n', ...
    'slice index', 'slice position (ground truth)', ...
    'slice position (noisy)');
fclose(fid);   
dlmwrite(csvfile, ...
    [(0:N-1)', t0', (t + ftot)'], '-append', 'delimiter', ',','precision', '%.3E');
system(['sed -i ''s/E+/E/g'' ' csvfile]);

%% run more iterations until achieving optimal reconstruction
for I = 6:42

    % impose Neumann conditions at the extremes, i.e. 2*delta(0,1) and
    % 2*delta(N-1,N-2) respectively
    fprev(1) = fpos(1);
    fpos(end) = fprev(end);
    
    % compute the transform to apply to each slice in this iteration
    f = alpha * (fprev + fpos);
    
    % compose with the previous transforms to obtain a total transform
    ftot = ftot + f;
    
    % update neighbour transformation step (note that we ignore the element
    % for the ghost slices)
    fprev(2:end) = fprev(2:end) - f(2:end) + f(1:end-1);
    fpos(1:end-1) = fpos(1:end-1) - f(1:end-1) + f(2:end);
    
    % compute the error between the reconstruction and the true shape
    err(I) = norm(t + ftot - t0);
    
   
end

% plot slice positions
hold off
plot(0:N-1, t + ftot, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
hold on
plot(0:N-1, t0, ':k', 'LineWidth', 2)
axis([0 100 -3 2])
xlabel('slice index')
ylabel('slice position')

% save figure
saveas(gca, [FIGDIR filesep 'registration-diffusion-toy-translation-iter-42.png'])

% save as CSV for Elsevier's Interactive Plots
csvfile = [FIGDIR filesep 'fig4c-registration-diffusion-toy-translation-iter-42.csv'];
fid = fopen(csvfile, 'w');
fprintf(fid, '%s,%s,%s\n', ...
    'slice index', 'slice position (ground truth)', ...
    'slice position (noisy)');
fclose(fid);   
dlmwrite(csvfile, ...
    [(0:N-1)', t0', (t + ftot)'], '-append', 'delimiter', ',','precision', '%.3E');
system(['sed -i ''s/E+/E/g'' ' csvfile]);

% run more iterations until achieving optimal reconstruction
for I = 43:1000

    % impose Neumann conditions at the extremes, i.e. 2*delta(0,1) and
    % 2*delta(N-1,N-2) respectively
    fprev(1) = fpos(1);
    fpos(end) = fprev(end);
    
    % compute the transform to apply to each slice in this iteration
    f = alpha * (fprev + fpos);
    
    % compose with the previous transforms to obtain a total transform
    ftot = ftot + f;
    
    % update neighbour transformation step (note that we ignore the element
    % for the ghost slices)
    fprev(2:end) = fprev(2:end) - f(2:end) + f(1:end-1);
    fpos(1:end-1) = fpos(1:end-1) - f(1:end-1) + f(2:end);
    
    % compute the error between the reconstruction and the true shape
    err(I) = norm(t + ftot - t0);
    
   
end

% plot slice positions
hold off
plot(0:N-1, t + ftot, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
hold on
plot(0:N-1, t0, ':k', 'LineWidth', 2)
axis([0 100 -3 2])
xlabel('slice index')
ylabel('slice position')

% save figure
saveas(gca, [FIGDIR filesep 'registration-diffusion-toy-translation-iter-1000.png'])

% save as CSV for Elsevier's Interactive Plots
csvfile = [FIGDIR filesep 'fig4d-registration-diffusion-toy-translation-iter-1000.csv'];
fid = fopen(csvfile, 'w');
fprintf(fid, '%s,%s,%s\n', ...
    'slice index', 'slice position (ground truth)', ...
    'slice position (noisy)');
fclose(fid);   
dlmwrite(csvfile, ...
    [(0:N-1)', t0', (t + ftot)'], '-append', 'delimiter', ',','precision', '%.3E');
system(['sed -i ''s/E+/E/g'' ' csvfile]);

%% run more iterations until achieving optimal reconstruction
for I = 1001:7000

    % impose Neumann conditions at the extremes, i.e. 2*delta(0,1) and
    % 2*delta(N-1,N-2) respectively
    fprev(1) = fpos(1);
    fpos(end) = fprev(end);
    
    % compute the transform to apply to each slice in this iteration
    f = alpha * (fprev + fpos);
    
    % compose with the previous transforms to obtain a total transform
    ftot = ftot + f;
    
    % update neighbour transformation step (note that we ignore the element
    % for the ghost slices)
    fprev(2:end) = fprev(2:end) - f(2:end) + f(1:end-1);
    fpos(1:end-1) = fpos(1:end-1) - f(1:end-1) + f(2:end);
    
    % compute the error between the reconstruction and the true shape
    err(I) = norm(t + ftot - t0);
    
   
end

% plot slice positions
hold off
plot(0:N-1, t + ftot, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
hold on
plot(0:N-1, t0, ':k', 'LineWidth', 2)
axis([0 100 -3 2])
xlabel('slice index')
ylabel('slice position')

% save figure
saveas(gca, [FIGDIR filesep 'registration-diffusion-toy-translation-iter-7000.png'])

% save as CSV for Elsevier's Interactive Plots
csvfile = [FIGDIR filesep 'fig4e-registration-diffusion-toy-translation-iter-7000.csv'];
fid = fopen(csvfile, 'w');
fprintf(fid, '%s,%s,%s\n', ...
    'slice index', 'slice position (ground truth)', ...
    'slice position (noisy)');
fclose(fid);   
dlmwrite(csvfile, ...
    [(0:N-1)', t0', (t + ftot)'], '-append', 'delimiter', ',','precision', '%.3E');
system(['sed -i ''s/E+/E/g'' ' csvfile]);

%% plot reconstruction error
hold off
semilogx(err/N, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
xlabel('stack sweeps')
ylabel('mean square error')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% synthetic problem, comparison with Gaffling 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coefficient
alpha = 0.49;

% number of slices
N = 100;

% place slices along true shape as sine wave
x = (0:(N-1))/N;
t = sin(2*pi*x);

% keep a copy of the true shape
t0 = t;

% random initialization of slice positions away from their perfect
% alignment
rng(10)
noise = rand(1, N);
noise = noise - mean(noise);
t = t + 2.4 * noise;

% plot slice positions
hold off
plot(0:N-1, t, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
hold on
plot(0:N-1, t0, ':k', 'LineWidth', 2)
axis([0 100 -3 2])
xlabel('slice index')
ylabel('slice position')

% compute registration between slices
fprev = -diff(t);
fpos = -diff(t(end:-1:1));
fpos = fpos(end:-1:1);

% pad with values so that computing the slice update step is more simple.
% After padding, for the i-th slice we can do alpha*(fprev(i)+fpos(i))
fprev = [NaN fprev];
fpos = [fpos NaN];

% init accumulated transform to apply to each slice
ftot_casero = zeros(1, N);

% NOTE: the optimal reconstruction for our method is ITER = 42

% allocate space for reconstruction error
err_casero = zeros(1, 7000);

%% run a lot of iterations for our method
for I = 1:7000

    % impose Neumann conditions at the extremes, i.e. 2*delta(0,1) and
    % 2*delta(N-1,N-2) respectively
    fprev(1) = fpos(1);
    fpos(end) = fprev(end);
    
    % compute the transform to apply to each slice in this iteration
    f = alpha * (fprev + fpos);
    
    % compose with the previous transforms to obtain a total transform
    ftot_casero = ftot_casero + f;
    
    % update neighbour transformation step (note that we ignore the element
    % for the ghost slices)
    fprev(2:end) = fprev(2:end) - f(2:end) + f(1:end-1);
    fpos(1:end-1) = fpos(1:end-1) - f(1:end-1) + f(2:end);
    
    % compute the error between the reconstruction and the true shape
    err_casero(I) = norm(t + ftot_casero - t0);
    
   
end

% init accumulated transform to apply to each slice
ftot_gaffling = zeros(1, N);

% allocate space for reconstruction error
err_gaffling = zeros(1, 7000);

% expand stack with dummy slices for Neumman boundary
fgaff = [NaN t NaN];

%% run a lot of iterations for Gaffling's Gauss-Seidel method
for I = 1:7000

    % impose Neumann conditions at the extremes
    fgaff(1) = fgaff(2);
    fgaff(end) = fgaff(end-1);
    
    % iteratively apply Gauss Seidel iterations (one sweep forth, one sweep
    % back)
    if (mod(I,2)) % odd sweeps
        for J = 2:N
            
            fgaff(J) = fgaff(J-1) + 0.5 * (fgaff(J+1) - fgaff(J-1));
            
        end
    else % even sweeps
        for J = N:-1:2
            
            fgaff(J) = fgaff(J+1) + 0.5 * (fgaff(J-1) - fgaff(J+1));
            
        end
    end
    
    % compute the error between the reconstruction and the true shape
    err_gaffling(I) = norm(fgaff(2:end-1) - t0);
   
end

%% plot reconstruction error
hold off
semilogx(err_casero/N, 'k', 'LineWidth', 2)
hold on
semilogx(err_gaffling/N, 'k--', 'LineWidth', 2)
set(gca, 'FontSize', 16)
xlabel('stack sweeps')
ylabel('mean square error')
legend('TDR', 'Gafflin2015', 'Location', 'NorthWest')

% save figure
saveas(gca, [FIGDIR filesep 'fig4f-registration-diffusion-toy-translation-err.png'])

% save as CSV for Elsevier's Interactive Plots
% too many rows for Elsevier's visualizer
csvfile = [FIGDIR filesep 'fig4f-registration-diffusion-toy-translation-err.csv'];
fid = fopen(csvfile, 'w');
fprintf(fid, '%s,%s,%s\n', ...
    'stack sweeps', 'MSE TDR', 'MSE Gafflin2015');
fclose(fid);   
dlmwrite(csvfile, ...
    [(1:7000)', err_casero'/N, err_gaffling'/N], '-append', 'delimiter', ',','precision', '%.3E');
system(['sed -i ''s/E+/E/g'' ' csvfile]);


% find how many iterations until solution worsens by 10%
[minerr_casero, Iopt] = min(err_casero);
idx = find(minerr_casero ./ err_casero > .9);
idx = idx(end);
disp(['Casero interval: ' num2str(idx) ' - ' num2str(Iopt) ' = ' num2str(idx-Iopt)])

[minerr_gaffling, Iopt] = min(err_gaffling);
idx = find(minerr_gaffling ./ err_gaffling > .9);
idx = idx(end);
disp(['Gaffling interval: ' num2str(idx) ' - ' num2str(Iopt) ' = ' num2str(idx-Iopt)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% synthetic problem adding drift, comparison with Gaffling 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coefficient
alpha = 0.49;

% number of slices
N = 100;

% place slices along true shape as sine wave
x = (0:(N-1))/N;
t = sin(2*pi*x);

% keep a copy of the true shape
t0 = t;

% % add a small drift component with zero mean
t = t + linspace(0, 1, N) - 0.5;

% random initialization of slice positions away from their perfect
% alignment
rng(10)
noise = rand(1, N);
noise = noise - mean(noise);
t = t + 2.4 * noise;

% compute registration between slices
fprev = -diff(t);
fpos = -diff(t(end:-1:1));
fpos = fpos(end:-1:1);

% pad with values so that computing the slice update step is more simple.
% After padding, for the i-th slice we can do alpha*(fprev(i)+fpos(i))
fprev = [NaN fprev];
fpos = [fpos NaN];

% init accumulated transform to apply to each slice
ftot_casero = zeros(1, N);

% NOTE: the optimal reconstruction for our method is ITER = 42

% allocate space for reconstruction error
err_casero = zeros(1, 7000);

%% run a lot of iterations for our method
best_err = Inf;
for I = 1:7000

    % impose Neumann conditions at the extremes, i.e. 2*delta(0,1) and
    % 2*delta(N-1,N-2) respectively
    fprev(1) = fpos(1);
    fpos(end) = fprev(end);
    
    % compute the transform to apply to each slice in this iteration
    f = alpha * (fprev + fpos);
    
    % compose with the previous transforms to obtain a total transform
    ftot_casero = ftot_casero + f;
    
    % update neighbour transformation step (note that we ignore the element
    % for the ghost slices)
    fprev(2:end) = fprev(2:end) - f(2:end) + f(1:end-1);
    fpos(1:end-1) = fpos(1:end-1) - f(1:end-1) + f(2:end);
    
    % compute the error between the reconstruction and the true shape
    err_casero(I) = norm(t + ftot_casero - t0);

    % keep a copy of the best reconstruction
    if (err_casero(I) < best_err)
        best_err = err_casero(I);
        best_casero = t + ftot_casero;
    end
   
end

% init accumulated transform to apply to each slice
ftot_gaffling = zeros(1, N);

% allocate space for reconstruction error
err_gaffling = zeros(1, 7000);

% expand stack with dummy slices for Neumman boundary
fgaff = [NaN t NaN];

%% run a lot of iterations for Gaffling's Gauss-Seidel method
best_err = Inf;
for I = 1:7000

    % impose Neumann conditions at the extremes
    fgaff(1) = fgaff(2);
    fgaff(end) = fgaff(end-1);
    
    % iteratively apply Gauss Seidel iterations (one sweep forth, one sweep
    % back)
    if (mod(I,2)) % odd sweeps
        for J = 2:N
            
            fgaff(J) = fgaff(J-1) + 0.5 * (fgaff(J+1) - fgaff(J-1));
            
        end
    else % even sweeps
        for J = N:-1:2
            
            fgaff(J) = fgaff(J+1) + 0.5 * (fgaff(J-1) - fgaff(J+1));
            
        end
    end
    
    % compute the error between the reconstruction and the true shape
    err_gaffling(I) = norm(fgaff(2:end-1) - t0);
   
    % keep a copy of the best reconstruction
    if (err_gaffling(I) < best_err)
        best_err = err_gaffling(I);
        best_gaffling = fgaff(2:end-1);
    end
   
end

% plot initial slice positions
hold off
plot(0:N-1, t, 'k', 'LineWidth', 2)
set(gca, 'FontSize', 16)
hold on
plot(0:N-1, t0, ':k', 'LineWidth', 2)
axis([0 100 -3 2])
xlabel('slice index')
ylabel('slice position')

% save figure
saveas(gca, [FIGDIR filesep 'fig5a-registration-diffusion-toy-translation-ground-truth-drift.png'])

% save as CSV for Elsevier's Interactive Plots
csvfile = [FIGDIR filesep 'fig5a-registration-diffusion-toy-translation-ground-truth-drift.csv'];
fid = fopen(csvfile, 'w');
fprintf(fid, '%s,%s,%s\n', ...
    'slice index', 'slice position (ground truth)', ...
    'slice position (noisy)');
fclose(fid);   
dlmwrite(csvfile, ...
    [(0:N-1)', t0', t'], '-append', 'delimiter', ',','precision', '%.3E');
system(['sed -i ''s/E+/E/g'' ' csvfile]);



%% plot reconstruction error
hold off
semilogx(err_casero/N, 'k', 'LineWidth', 2)
hold on
semilogx(err_gaffling/N, 'k--', 'LineWidth', 2)
set(gca, 'FontSize', 16)
xlabel('stack sweeps')
ylabel('mean square error')
legend('TDR', 'Gafflin2015', 'Location', 'NorthWest')

% save figure
saveas(gca, [FIGDIR filesep 'fig5b-registration-diffusion-toy-translation-err-drift.png'])

% save as CSV for Elsevier's Interactive Plots
% too many rows for Elsevier's visualizer
csvfile = [FIGDIR filesep 'fig5b-registration-diffusion-toy-translation-err-drift.csv'];
fid = fopen(csvfile, 'w');
fprintf(fid, '%s,%s,%s\n', ...
    'stack sweeps', 'MSE TDR', 'MSE Gafflin2015');
fclose(fid);   
dlmwrite(csvfile, ...
    [(1:7000)', err_casero'/N, err_gaffling'/N], '-append', 'delimiter', ',','precision', '%.3E');
system(['sed -i ''s/E+/E/g'' ' csvfile]);


% find how many iterations until solution worsens by 10%
[minerr_casero, Iopt] = min(err_casero);
idx = find(minerr_casero ./ err_casero > .9);
idx = idx(end);
disp(['Casero interval: ' num2str(idx) ' - ' num2str(Iopt) ' = ' num2str(idx-Iopt)])

[minerr_gaffling, Iopt] = min(err_gaffling);
idx = find(minerr_gaffling ./ err_gaffling > .9);
idx = idx(end);
disp(['Gaffling interval: ' num2str(idx) ' - ' num2str(Iopt) ' = ' num2str(idx-Iopt)])
