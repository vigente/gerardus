function th2bf = aux_histology2blockface(filenameh, filenamebf, outdir, polymask, ellipmask)

% use auto estimation of background colour in transformix
optTransformix.AutoDefaultPixelValue = true;

% load histology and blockface
imh0 = scimat_load(filenameh);
imbf = scimat_load(filenamebf);

% plot images
subplot(2, 1, 1)
hold off
imagesc(imbf.data)
subplot(2, 1, 2)
hold off
imagesc(squeeze(imh0.data))

% preprocessing of blockface and histology to prepare them for
% registration
imh.axis = imh0.axis;
imh.rotmat = imh0.rotmat;
[imh.data, imbf.data] = histology_blockface_preprocessing(...
    squeeze(imh0.data), imbf.data, polymask, ellipmask);

% plot images
subplot(2, 1, 1)
hold off
imagesc(imbf.data)
subplot(2, 1, 2)
hold off
imagesc(imh.data)

% convolution registration of histology to blockface
th2bf = regmatchedfilt(imbf, imh, (-45:45)/180*pi);

% apply rigid transform to original histology
imh2bf = transformix(th2bf, imh0, optTransformix);

% plot histology and blockface overlap
subplot(2, 1, 1)
hold off
imagesc(imfuse(imbf.data, squeeze(imh2bf.data)))
drawnow

% save registered image
[~, name] = fileparts(filenameh);
scimat_save([outdir filesep name '.mha'], imh2bf);
