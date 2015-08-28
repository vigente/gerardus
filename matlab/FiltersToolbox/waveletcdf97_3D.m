function X = waveletcdf97_3D(X, Level)
%WAVELETCDF97  Cohen-Daubechies-Feauveau 9/7 wavelet transform in 3D.
%   Y = WAVELETCDF97(X, L) decomposes X with L stages of the
%   Cohen-Daubechies-Feauveau (CDF) 9/7 wavelet.  For the
%   inverse transform, WAVELETCDF97(X, -L) inverts L stages.
%   Filter boundary handling is half-sample symmetric.
%
%   X may be of any size; it need not have size divisible by 2^L.
%   For example, if X has length 9, one stage of decomposition
%   produces a lowpass subband of length 5 and a highpass subband
%   of length 4.  Transforms of any length have perfect
%   reconstruction (exact inversion).
%
%
%   Example:
%   Y = waveletcdf97(X, 5);    % Transform image X using 5 stages
%   R = waveletcdf97(Y, -5);   % Reconstruct from Y

% This file is a 3D extension of waveletcdf97 by Pascal Getreuer 2004-2006

% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2015 University of Oxford
% Version: 0.1.0
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



if nargin < 2, error('Not enough input arguments.'); end
if ndims(X) ~= 3, error('Input must be a 3D array.'); end
if any(size(Level) ~= 1), error('Invalid transform level.'); end

N1 = size(X,1);
N2 = size(X,2);
N3 = size(X,3);

% Lifting scheme filter coefficients for CDF 9/7
LiftFilter = [-1.5861343420693648,-0.0529801185718856,0.8829110755411875,0.4435068520511142];
ScaleFactor = 1.1496043988602418;

S1 = LiftFilter(1);
S2 = LiftFilter(2);
S3 = LiftFilter(3);
ExtrapolateOdd = -2*[S1*S2*S3,S2*S3,S1+S3+3*S1*S2*S3]/(1+2*S2*S3);

LiftFilter = LiftFilter([1,1],:);

if Level >= 0   % Forward transform
   for k = 1:Level
      M1 = ceil(N1/2);
      M2 = ceil(N2/2);
      M3 = ceil(N3/2);
      
      %%% Transform along columns %%%
      if N1 > 1         
         RightShift = [2:M1,M1];
         X0 = X(1:2:N1,1:N2,1:N3);

         % Apply lifting stages
         if rem(N1,2)
            X1 = [X(2:2:N1,1:N2,1:N3);X0(M1-1,:,:)*ExtrapolateOdd(1)...
                  + X(N1-1,1:N2,1:N3)*ExtrapolateOdd(2)...
                  + X0(M1,:,:)*ExtrapolateOdd(3)]...
               + filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
               X0(1,:,:)*LiftFilter(1,1),1);
         else
            X1 = X(2:2:N1,1:N2,1:N3) ...
               + filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
               X0(1,:,:)*LiftFilter(1,1),1);
         end

         X0 = X0 + filter(LiftFilter(:,2),1,...
            X1,X1(1,:,:)*LiftFilter(1,2),1);
         X1 = X1 + filter(LiftFilter(:,3),1,...
            X0(RightShift,:,:),X0(1,:,:)*LiftFilter(1,3),1);
         X0 = X0 + filter(LiftFilter(:,4),1,...
            X1,X1(1,:,:)*LiftFilter(1,4),1);

         if rem(N1,2)
            X1(M1,:,:) = [];
         end

         X(1:N1,1:N2,1:N3) = [X0*ScaleFactor;X1/ScaleFactor];
      end

      %%% Transform along rows %%%
      if N2 > 1
         RightShift = [2:M2,M2];
         X0 = permute(X(1:N1,1:2:N2,1:N3),[2,1,3]);

         % Apply lifting stages
         if rem(N2,2)
            X1 = permute([X(1:N1,2:2:N2,1:N3),X(1:N1,N2-2,1:N3)*ExtrapolateOdd(1)...
                  + X(1:N1,N2-1,1:N3)*ExtrapolateOdd(2) ...
                  + X(1:N1,N2,1:N3)*ExtrapolateOdd(3)],[2,1,3])...
               + filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
               X0(1,:,:)*LiftFilter(1,1),1);
         else
            X1 = permute(X(1:N1,2:2:N2,1:N3),[2,1,3]) ...
               + filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
               X0(1,:,:)*LiftFilter(1,1),1);
         end

         X0 = X0 + filter(LiftFilter(:,2),1,...
            X1,X1(1,:,:)*LiftFilter(1,2),1);
         X1 = X1 + filter(LiftFilter(:,3),1,...
            X0(RightShift,:,:),X0(1,:,:)*LiftFilter(1,3),1);
         X0 = X0 + filter(LiftFilter(:,4),1,...
            X1,X1(1,:,:)*LiftFilter(1,4),1);

         if rem(N2,2)
            X1(M2,:,:) = [];
         end

         X(1:N1,1:N2,1:N3) = permute([X0*ScaleFactor;X1/ScaleFactor],[2,1,3]);
      end
      
      %%% Transform along slices %%%
      if N3 > 1
         RightShift = [2:M3,M3];
         X0 = permute(X(1:N1,1:N2,1:2:N3),[3,2,1]);

         % Apply lifting stages
         if rem(N3,2)
            X1 = permute(cat(3, X(1:N1,1:N2,2:2:N3),X(1:N1,1:N2,N3-2)*ExtrapolateOdd(1)...
                  + X(1:N1,1:N2,N3-1)*ExtrapolateOdd(2) ...
                  + X(1:N1,1:N2,N3)*ExtrapolateOdd(3)),[3,2,1])...
               + filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
               X0(1,:,:)*LiftFilter(1,1),1);
         else
            X1 = permute(X(1:N1,1:N2,2:2:N3),[3,2,1]) ...
               + filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
               X0(1,:,:)*LiftFilter(1,1),1);
         end

         X0 = X0 + filter(LiftFilter(:,2),1,...
            X1,X1(1,:,:)*LiftFilter(1,2),1);
         X1 = X1 + filter(LiftFilter(:,3),1,...
            X0(RightShift,:,:),X0(1,:,:)*LiftFilter(1,3),1);
         X0 = X0 + filter(LiftFilter(:,4),1,...
            X1,X1(1,:,:)*LiftFilter(1,4),1);

         if rem(N3,2)
            X1(M3,:,:) = [];
         end

         X(1:N1,1:N2,1:N3) = permute([X0*ScaleFactor;X1/ScaleFactor],[3,2,1]);
      end
      
      N1 = M1;
      N2 = M2;
      N3 = M3;
   end
else           % Inverse transform
   for k = 1+Level:0
      M1 = ceil(N1*pow2(k));
      M2 = ceil(N2*pow2(k));
      M3 = ceil(N3*pow2(k));
      

      %%% Inverse transform along rows %%%
      if M2 > 1
         Q = ceil(M2/2);
         RightShift = [2:Q,Q];
         X1 = permute(X(1:M1,Q+1:M2,1:M3)*ScaleFactor,[2,1,3]);

         if rem(M2,2)
            X1(Q,1,1) = 0;
         end

         % Undo lifting stages
         X0 = permute(X(1:M1,1:Q,1:M3)/ScaleFactor,[2,1,3]) ...
            - filter(LiftFilter(:,4),1,X1,X1(1,:,:)*LiftFilter(1,4),1);
         X1 = X1 - filter(LiftFilter(:,3),1,X0(RightShift,:,:),...
            X0(1,:,:)*LiftFilter(1,3),1);
         X0 = X0 - filter(LiftFilter(:,2),1,X1,...
            X1(1,:,:)*LiftFilter(1,2),1);
         X1 = X1 - filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
            X0(1,:,:)*LiftFilter(1,1),1);

         if rem(M2,2)
            X1(Q,:,:) = [];
         end

         X(1:M1,[1:2:M2,2:2:M2],1:M3) = permute([X0;X1],[2,1,3]);
      end

      %%% Inverse transform along columns %%%
      if M1 > 1
         Q = ceil(M1/2);
         RightShift = [2:Q,Q];
         X1 = X(Q+1:M1,1:M2,1:M3)*ScaleFactor;

         if rem(M1,2)
            X1(Q,1,1) = 0;
         end

         % Undo lifting stages
         X0 = X(1:Q,1:M2,1:M3)/ScaleFactor ...
            - filter(LiftFilter(:,4),1,X1,X1(1,:,:)*LiftFilter(1,4),1);
         X1 = X1 - filter(LiftFilter(:,3),1,X0(RightShift,:,:),...
            X0(1,:,:)*LiftFilter(1,3),1);
         X0 = X0 - filter(LiftFilter(:,2),1,X1,...
            X1(1,:,:)*LiftFilter(1,2),1);
         X1 = X1 - filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
            X0(1,:,:)*LiftFilter(1,1),1);

         if rem(M1,2)
            X1(Q,:,:) = [];
         end

         X([1:2:M1,2:2:M1],1:M2,1:M3) = [X0;X1];
      end
      
      %%% Inverse transform along slices %%%
      if M3 > 1
         Q = ceil(M3/2);
         RightShift = [2:Q,Q];
         X1 = permute(X(1:M1,1:M2,Q+1:M3)*ScaleFactor,[3,2,1]);

         if rem(M3,2)
            X1(Q,1,1) = 0;
         end

         % Undo lifting stages
         X0 = permute(X(1:M1,1:M2,1:Q)/ScaleFactor,[3,2,1]) ...
            - filter(LiftFilter(:,4),1,X1,X1(1,:,:)*LiftFilter(1,4),1);
         X1 = X1 - filter(LiftFilter(:,3),1,X0(RightShift,:,:),...
            X0(1,:,:)*LiftFilter(1,3),1);
         X0 = X0 - filter(LiftFilter(:,2),1,X1,...
            X1(1,:,:)*LiftFilter(1,2),1);
         X1 = X1 - filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
            X0(1,:,:)*LiftFilter(1,1),1);

         if rem(M3,2)
            X1(Q,:,:) = [];
         end

         X(1:M1,1:M2,[1:2:M3,2:2:M3]) = permute([X0;X1],[3,2,1]);
      end
      
   end
end
