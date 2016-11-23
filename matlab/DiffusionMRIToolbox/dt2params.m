function [FA, ADC, VectorField, EigVals] = dt2params(DT)

sz = size(DT);
M = reshape(DT, [prod(sz(1:end-1)), sz(end)]);
mask = M(:,1) ~= 0;   

M = M(mask(:), :);

% initialise variables
FA_mask = zeros(size(M,1), 1);
ADC_mask = zeros(size(M,1),1);
VectorF_mask = zeros(size(M,1), 3);
VectorF2_mask = zeros(size(M,1), 3);
VectorF3_mask = zeros(size(M,1), 3);
EigVals_mask = zeros(size(M,1), 3);

parfor i = 1:size(M,1)

    Mi = M(i,:);

    % The DiffusionTensor (Remember it is a symetric matrix,
    % thus for instance Dxy == Dyx)
    DiffusionTensor=[Mi(1) Mi(2) Mi(3); Mi(2) Mi(4) Mi(5); Mi(3) Mi(5) Mi(6)];

    % Calculate the eigenvalues and vectors, and sort the 
    % eigenvalues from small to large
    [EigenVectors,D]=eig(DiffusionTensor); 
    EigenValues=diag(D);
    [~,index]=sort(EigenValues); 
    EigenValues=EigenValues(index); 
    EigenVectors=EigenVectors(:,index);
    EigenValues_old=EigenValues;

    % Regulating of the eigen values (negative eigenvalues are
    % due to noise and other non-idealities of MRI)
    EigenValues=abs(EigenValues);

    % Apparent Diffuse Coefficient
    ADCv = mean(EigenValues);%(EigenValues(1)+EigenValues(2)+EigenValues(3))/3;

    % Fractional Anistropy (2 different definitions exist)
    % First FA definition:
    %FAv=(1/sqrt(2))*( sqrt((EigenValues(1)-EigenValues(2)).^2+(EigenValues(2)-EigenValues(3)).^2+(EigenValues(1)-EigenValues(3)).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
    % Second FA definition:
    FA_mask(i)=sqrt(1.5)*( sqrt((EigenValues(1)-ADCv).^2+(EigenValues(2)-ADCv).^2+(EigenValues(3)-ADCv).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
    ADC_mask(i)=ADCv;
    VectorF_mask(i,:)=EigenVectors(:,3)*EigenValues_old(3);
    VectorF2_mask(i,:)=EigenVectors(:,2)*EigenValues_old(2);
    VectorF3_mask(i,:)=EigenVectors(:,1)*EigenValues_old(1);

    EigVals_mask(i,:) = EigenValues_old([3 2 1]);


end



% reshape to match input dimensions, with zero filling outside the mask
VectorF = zeros([prod(sz(1:end-1)), 3]);
VectorF(mask(:), :) = VectorF_mask;
VectorF = reshape(VectorF, [sz(1:end-1), 3]);
VectorF2 = zeros([prod(sz(1:end-1)), 3]);
VectorF2(mask(:), :) = VectorF2_mask;
VectorF2 = reshape(VectorF2, [sz(1:end-1), 3]);
VectorF3 = zeros([prod(sz(1:end-1)), 3]);
VectorF3(mask(:), :) = VectorF3_mask;
VectorF3 = reshape(VectorF3, [sz(1:end-1), 3]);
VectorField = cat(ndims(VectorF)+1, VectorF, VectorF2, VectorF3);

EigVals = zeros([prod(sz(1:end-1)), 3]);
EigVals(mask(:), :) = EigVals_mask;
EigVals = reshape(EigVals, [sz(1:end-1), 3]);

% quick hack because reshape sz needs to be at least length 2
if length(sz) == 2
    sz = [sz(1), 1, sz(2)];
end

FA = zeros([prod(sz(1:end-1)), 1]);
FA(mask(:)) = FA_mask;
FA = reshape(FA, sz(1:end-1));
ADC = zeros([prod(sz(1:end-1)), 1]);
ADC(mask(:)) = ADC_mask;
ADC = reshape(ADC, sz(1:end-1));
