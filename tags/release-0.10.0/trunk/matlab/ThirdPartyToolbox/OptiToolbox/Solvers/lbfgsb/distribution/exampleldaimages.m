% Script demonstrating the latent Dirichlet allocation model on synthetic
% data.

% Script parameters.
K      = 3;   % The number of topics.
D      = 25;  % The number of documents (images).
L      = 16;  % The length of each document.
height = 4;   % The height of the images (in pixels).
width  = 4;   % The width of the images (in pixels).
nr     = 4;   % Number of rows of documents to display.
nc     = 6;   % Number of columns of documents to display.
ns     = 3;   % Number of samples to display.

% The size of the vocabulary (W) is equal to the number of pixels.
W = height*width;
L = repmat(L,1,D);

% Set a uniform prior distribution over the topics.
nu = ones(1,K)/K;

% Set a uniform prior distribution over the words.
eta = ones(1,W)/W;

% Randomly generate the topics and the documents.
[topics w] = createsyntheticdata(nu,eta,K,W,L);
wc = zeros(D,W);
for d = 1:D
  wc(d,:) = hist(w{d},1:W)';
end

% Display the topics.
figure(1);
set(gcf,'MenuBar','none');
set(gcf,'NumberTitle','off');
set(gcf,'Name','True topics');
set(gcf,'Color','white');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) 228 65]);
colormap(gray);
maxval = max(topics(:));
for k = 1:K
  subplot(1,K,k);
  imagesc(reshape(topics(:,k),height,width),[0 maxval]);
  axis off
end

% Show some of the documents that were generated.
figure(2);
set(gcf,'MenuBar','none');
set(gcf,'NumberTitle','off');
set(gcf,'Name','Data');
set(gcf,'Color','white');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) 222 355]);
colormap(gray);
maxval = max(wc(:));
d = 0;
for c = 1:nc
  for r = 1:nr
    d = d + 1;
    subplot(nc,nr,d);
    imagesc(reshape(wc(d,:),height,width),[0 maxval]);
    axis off
  end
end
drawnow

% Compute the mean field estimate.
tic;
fprintf('Finding solution to the mean field lower bound.\n');
[xi0 gamma0 Phi0]  = mfldainit(W,K,w);
[LnZ xi gamma Phi] = mflda(nu,eta,w,xi0,gamma0,Phi0);
fprintf('Optimization took %i seconds.\n', round(toc));

% Show the estimated topics.
figure(3);
set(gcf,'MenuBar','none');
set(gcf,'NumberTitle','off');
set(gcf,'Name','Topics (mean field)');
set(gcf,'Color','white');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) 197 226]);
colormap(gray);

% Repeat for every sample to display.
for i = 1:ns
  beta = zeros(W,K);
  for k = 1:K
    beta(:,k) = dirichletrnd(xi(:,k));
  end
  maxval = max(beta(:));

  for k = 1:K
    subplot(ns,K,K*(i-1)+k);
    imagesc(reshape(beta(:,k),height,width),[0 maxval]);
    axis off
    title(' ');
  end
  subplot(ns,K,K*(i-1)+1);
  title(sprintf('Sample #%i',i));
end
