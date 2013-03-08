function thres_img=SB_intrp(data_SBintrp, inter_slice_no)
%Function that performs the shape-base interpolation following Herman_2005
%paper.
%INPUT
%data_SBintrp = the 3D image in form of a 3D matlab matrix
%inter_slice_no= the number of slices to generate by interpolation between
%each couple of original slices
%OUTPUT
%thres_img = 3D matlab matrix with the 

% First distance map

%Rules
% 1) Replace all 0's by -99
% 2) Replace all 1's by +99
%Exceptions
% 3) If a 0 shares an edge with a 1 replace it with -5
% 4) If a 1 shares an edge with a 0 replace it with +5

%Pad the image in the xy plane with zeros [x y]
pad_data_SBintrp=padarray(data_SBintrp,[1 1], 0);

new_data_SBintrp=pad_data_SBintrp;

[size_x,size_y,size_z]=size(new_data_SBintrp);

display('SB initialisation')
tic
for slice=1:size_z
    pad_data_slice=pad_data_SBintrp(:,:,slice);
    for k=2:(size_x-1) %excluding the padding border
        
        for m=2:(size_y-1) %excluding the padding border
            
            tmp_pixel=pad_data_slice(k,m);
            neighbours_tmp=pixel_neighbour_vector(pad_data_slice,k,m,4);
            
            %zero pixel
            if tmp_pixel==0
                neighbour_ones=find(neighbours_tmp==1);
                %if there is at least one nonzero entry in the neighbourhood of
                %the zero pixel
                if isempty(neighbour_ones)==0 % this means it's not empty
                    
                    new_data_SBintrp(k,m,slice)=-5;
                else
                    % all the neighbours of the zero pixel all zeros
                    new_data_SBintrp(k,m,slice)=-99;
                end
                
                %one-pixel
            elseif tmp_pixel==1
                neighbour_zeros=find(neighbours_tmp==0);
                % if there is at least one zero entry in the neighbourhood of
                % the 1 pixel
                if isempty(neighbour_zeros)==0 %this means is not empty
                    
                    new_data_SBintrp(k,m,slice)=5;
                else
                    % all the neighbours of the 1 pixel all ones
                    new_data_SBintrp(k,m,slice)=99;
                    
                end
            end
        end
    end
end

toc

%borders are still zero, convert them to -99
new_data_SBintrp(new_data_SBintrp==0)=-99;


% CHAMFER PROCESS STEP 1

chamfer_matrix1=new_data_SBintrp;
%in the first step of the process we use the 3x3 upper near optimal city
%block
block_type=1;

display('First chamfer step')
tic
for slice=1:size_z
    pad_data_slice=chamfer_matrix1(:,:,slice);
    
    for k=1:size_x
        for m=1:size_y
            
            mynewpixel = city_block(pad_data_slice,k,m,block_type);
            %updating the slice
            pad_data_slice(k,m)=mynewpixel;
            
            
        end
    end
    % put the updated 3D slice back in the original chamfer stack
    chamfer_matrix1(:,:,slice)=pad_data_slice;
end
toc


%CHAMFER STEP 2

chamfer_matrix2=chamfer_matrix1;
%in the second step of the process we use the 3x3 lower near optimal city
%block
block_type=-1;

display('Second chamfer step')
tic
for slice=1:size_z
    pad_data_slice=chamfer_matrix2(:,:,slice);
    
    for k=1:size_x
        for m=1:size_y
            
            %in order to reverse the span of the image
            kk=size_x-k+1;
            mm=size_y-m+1;
            mynewpixel = city_block(pad_data_slice,kk,mm,block_type);
            %updating the slice
            pad_data_slice(kk,mm)=mynewpixel;
            
            
        end
    end
    % put the updated 3D slice back in the original chamfer stack
    chamfer_matrix2(:,:,slice)=pad_data_slice;
end
toc


% CUBIC INTERPOLATION STEP




inter_slice_sp=zeros(size_x,size_y,1);

%linear inbetween interpolation
x0=(1:size_x)';
y0=(1:size_y)';
z0= linspace(1,2,inter_slice_no)';
[xf,yf,zf]=meshgrid(x0,y0,z0);


display('Two-slices-per-time cubic interpolation step')
for slice=1:size_z-1
display(slice)
tic
initial_final_slice_sp=chamfer_matrix2(:,:,(slice:slice+1));
tmp_inter_slice_sp=interp3(initial_final_slice_sp,xf,yf,zf,'spline');
inter_slice_sp=cat(3,inter_slice_sp,tmp_inter_slice_sp(:,:,2:end));
toc
end

% GET THE BINARY

thres_img=inter_slice_sp;
thres_img(thres_img<=5)=0;
thres_img(thres_img>5)=1;

end


% Complementary functions

function neighbours=pixel_neighbour_vector(matrix,k,m, neighbour_type)


[size_x,size_y]=size(matrix);

if k>1 & k < size_x & m>1 & m < size_y
    %inner pixels
    if neighbour_type==4
        neighbours=[matrix(k-1,m), matrix(k, m-1), matrix(k, m+1),...
            matrix(k+1,m) ];
        
    elseif neighbour_type==8
        neighbours=[matrix(k-1,m-1:m+1), matrix(k, m-1), matrix(k, m+1),...
            matrix(k+1,m-1:m+1) ];
    end
    
else
    
    %upper-left vertex (1,1)
    if k==1 & m==1
        
        if neighbour_type==4
            neighbours=[ matrix(k,m+1), matrix(k+1,m) ]; %only 2 neighbours
            
        elseif neighbour_type==8
            neighbours=[0,0,0,0, matrix(k,m+1),0, matrix(k+1,m:m+1) ]; %only 3 neighbours
        end
   
        
       
        
    %upper-right vertex (1,end)
    elseif k==1 & m==size_y
        
        if neighbour_type==4
            neighbours=[ matrix(k, m-1), matrix(k+1,m) ]; %only 2 neighbours
            
        elseif neighbour_type==8
            neighbours=[ 0,0,0,matrix(k, m-1),0,...
                matrix(k+1,m-1:m), 0]; %only 3 neighbours
        end
     
        
        %upper border (1,1:end-1)
    elseif k==1 & m>1 & m< size_y
        
         if neighbour_type==4
            neighbours=[ matrix(k, m-1), matrix(k, m+1), matrix(k+1,m) ];
            
         elseif neighbour_type==8   
            neighbours=[0,0,0 matrix(k, m-1), matrix(k, m+1),...
                    matrix(k+1,m-1:m+1) ]; %5 neighbours
         end
        
    
    
    
%bottom-left vertex (end,1)   
    elseif k==size_x & m==1
    
        if neighbour_type==4
            neighbours=[matrix(k-1,m), matrix(k, m+1) ];
            
        elseif neighbour_type==8
            neighbours=[0,matrix(k-1,m:m+1),0, matrix(k, m+1), 0, 0,0];
        end
    
        
     
        
        %bottom-right vertex (end,end)
    elseif k==size_x & m==size_y
 
        if neighbour_type==4
            neighbours=[matrix(k-1,m), matrix(k, m-1)];
            
        elseif neighbour_type==8
            neighbours=[matrix(k-1,m-1:m),0, matrix(k, m-1), 0, 0,0,0];
        end
        %bottom border (end,1:end-1)
    elseif  k==size_x & m>1 ^ m<size_y
        
     if neighbour_type==4
        neighbours=[matrix(k-1,m), matrix(k, m-1), matrix(k, m+1)];
        
    elseif neighbour_type==8
        neighbours=[matrix(k-1,m-1:m+1), matrix(k, m-1), matrix(k, m+1),...
            0, 0,0 ];
    end
        
    
    
    %all the inner values of k=2:end-1
    elseif k>1 && k < size_x && m==1
    
        if neighbour_type==4
        neighbours=[matrix(k-1,m), matrix(k, m+1),  matrix(k+1,m) ];
        
    elseif neighbour_type==8
        neighbours=[0,matrix(k-1,m:m+1), 0, matrix(k, m+1),...
           0, matrix(k+1,m:m+1) ];
        end  
       
        %final pixel of the row
    elseif k>1 && k < size_x && m==size_y
        
        if neighbour_type==4
            neighbours=[matrix(k-1,m), matrix(k, m-1),  matrix(k+1,m) ];
            
        elseif neighbour_type==8
            neighbours=[matrix(k-1,m-1:m), 0,matrix(k, m-1), 0,...
                matrix(k+1,m-1:m), 0 ];
        end
        
        
  
    end
end  
end



function mynewpixel = city_block(matrix,k,m,block_type)

%function that applies the city-block chamfer distance calculation to the
%pixel (k,m) of the 2D image matrix. The block_type = 1 if the 3x3
%upper near optimal city block matrix is applied, while block_type =-1 if
%the lower near optimal one is used. 
operation=0; 

%define the current pixel
mypixel=matrix(k,m);

[size_x,size_y]=size(matrix);

%decide if at mypixel the sum or difference has to be applied
if mypixel >0 
    
    if mypixel ==5
        mynewpixel=5;
    else
    operation=1;
    end
            
elseif mypixel <=0 
    
    if mypixel ==-5
        mynewpixel=-5;
    else
    operation=-1;
    end
end

%get the neighbours by applying the function submatrix. The neighbours
%vector follows the order shown below to go around the given pixel (k,m):
%(k-1,m-1),(k-1,m),(k-1,m+1),(k, m-1),(k, m+1),(k+1,m-1),(k+1:m),(k+1:m+1)

neighbours=pixel_neighbour_vector(matrix,k,m,8); %mypixel=matrix(k,m)

if operation == 1 % sum
    if block_type == 1 %upper near optimal template
        
        mynewpixel=cityblock_sum_up(mypixel,neighbours,k,m,size_x,size_y);
        
    elseif block_type == - 1 %lower near optimal template
        
        mynewpixel=cityblock_sum_down(mypixel,neighbours,k,m,size_x,size_y);
        
    end

elseif operation == -1 %difference
    if block_type == 1 %upper near optimal template
       
       mynewpixel=cityblock_diff_up(mypixel,neighbours,k,m,size_x,size_y);
       
    elseif block_type == - 1 %lower near optimal template
        
      mynewpixel=cityblock_diff_down(mypixel,neighbours,k,m,size_x,size_y); 
        
    end
end
end

function mynewpixel=cityblock_diff_down(mypixel,neighbours,k,m,size_x,size_y)

% city block applied from top to bottom from left to right
%it expects a neighbours vector of 8 elements
if length(neighbours)~=8
    display('not 8')
end


if k==1
    
    %upper-left vertex (1,1)
    if m==1
        
         city_b= [neighbours(5)-10;neighbours(7)-10;...
             neighbours(8)-14; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
        %upper-right vertex (1,end)
    elseif m==size_y
        
        city_b= [neighbours(6)-14; neighbours(7)-10;mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
        
        %upper border (1,1:end-1)
    else
        
        city_b= [neighbours(5)-10;neighbours(6)-14; neighbours(7)-10;...
            neighbours(8)-14; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
    end
    
    
    
elseif k==size_x
    
    %bottom-left vertex (end,1)
    if m==1
        
        city_b= [neighbours(5)-10; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
        %bottom-right vertex (end,end)
    elseif m==size_y
       
        mynewpixel=mypixel;
        
        %bottom border (end,1:end-1)
    else
        
        city_b= [neighbours(5)-10; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
    end
    
    %all the inner values of k=2:end-1
else
    
    %beginning pixel of the row
    if m==1
        
        city_b= [neighbours(5)-10;neighbours(7)-10;...
                 neighbours(8)-14; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
        %final pixel of the row
    elseif m==size_y
        
        city_b= [neighbours(6)-14;neighbours(7)-10; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
        
        %all the inner pixels of the image
    else
        
        city_b= [neighbours(5)-10;neighbours(6)-14;neighbours(7)-10;...
            neighbours(8)-14; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
    end
    
end
end

function mynewpixel=cityblock_diff_up(mypixel,neighbours,k,m,size_x,size_y)

% city block applied from top to bottom from left to right
if length(neighbours)~=8
    display('not 8')
end

if k==1
    
    %upper-left vertex (1,1)
    if m==1
        
        mynewpixel=mypixel;
        
        %upper-right vertex (1,end)
    elseif m==size_y
        
        city_b= [neighbours(4)-10; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
        
        %upper border (1,1:end-1)
    else
        
        city_b= [neighbours(4)-10; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
    end
    
    
    
elseif k==size_x
    
    %bottom-left vertex (end,1)
    if m==1
        
        city_b= [neighbours(2)-10;neighbours(3)-14; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
        %bottom-right vertex (end,end)
    elseif m==size_y
        
        city_b= [neighbours(1)-14;neighbours(2)-10;...
            neighbours(4)-14; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
        %bottom border (end,1:end-1)
    else
        
        city_b= [neighbours(1)-14;neighbours(2)-10;neighbours(3)-14;...
            neighbours(4)-10; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
    end
    
    %all the inner values of k=2:end-1
else
    
    %beginning pixel of the row
    if m==1
        
        city_b= [neighbours(2)-10;neighbours(3)-14; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
        %final pixel of the row
    elseif m==size_y
        
        city_b= [neighbours(1)-14;neighbours(2)-10;neighbours(4)-10; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
        
        
        %all the inner pixels of the image
    else
        
        city_b= [neighbours(1)-14;neighbours(2)-10;neighbours(3)-14;...
            neighbours(4)-10; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=max(city_b(:));
    end
    
end
end

function mynewpixel=cityblock_sum_down(mypixel,neighbours,k,m,size_x,size_y)

% city block applied from top to bottom from left to right
%it expects a neighbours vector of 8 elements
if length(neighbours)~=8
    display('not 8')
end


if k==1
    
    %upper-left vertex (1,1)
    if m==1
        
         city_b= [neighbours(5)+10;neighbours(7)+10;...
             neighbours(8)+14; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
        %upper-right vertex (1,end)
    elseif m==size_y
        
        city_b= [neighbours(6)+14; neighbours(7)+10;mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
        
        %upper border (1,1:end-1)
    else
        
        city_b= [neighbours(5)+10;neighbours(6)+14; neighbours(7)+10;...
            neighbours(8)+14; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
    end
    
    
    
elseif k==size_x
    
    %bottom-left vertex (end,1)
    if m==1
        
        city_b= [neighbours(5)+10; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
        %bottom-right vertex (end,end)
    elseif m==size_y
       
        mynewpixel=mypixel;
        
        %bottom border (end,1:end-1)
    else
        
        city_b= [neighbours(5)+10; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
    end
    
    %all the inner values of k=2:end-1
else
    
    %beginning pixel of the row
    if m==1
        
        city_b= [neighbours(5)+10;neighbours(7)+10;...
                 neighbours(8)+14; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
        %final pixel of the row
    elseif m==size_y
        
        city_b= [neighbours(6)+14;neighbours(7)+10; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
        
        %all the inner pixels of the image
    else
        
        city_b= [neighbours(5)+10;neighbours(6)+14;neighbours(7)+10;...
            neighbours(8)+14; mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
    end
    
end
end

function mynewpixel=cityblock_sum_up(mypixel,neighbours,k,m,size_x,size_y)

% city block applied from top to bottom from left to right
%it expects a neighbours vector of 8 elements
if length(neighbours)~=8
    display('not 8')
end


if k==1
    
    %upper-left vertex (1,1)
    if m==1
        
        mynewpixel=mypixel;
        
        %upper-right vertex (1,end)
    elseif m==size_y
        
        city_b= [10+neighbours(4); mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
        
        %upper border (1,1:end-1)
    else
        
        city_b= [10+neighbours(4); mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
    end
    
    
    
elseif k==size_x
    
    %bottom-left vertex (end,1)
    if m==1
        
        city_b= [10+neighbours(2);14+neighbours(3); mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
        %bottom-right vertex (end,end)
    elseif m==size_y
        city_b= [14+ neighbours(1);10+neighbours(2);...
            10+neighbours(4); mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
        %bottom border (end,1:end-1)
    else
        
        city_b= [14+ neighbours(1);10+neighbours(2);14+neighbours(3);...
            10+neighbours(4); mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
    end
    
    %all the inner values of k=2:end-1
else
    
    %beginning pixel of the row
    if m==1
        
        city_b= [10+neighbours(2);14+neighbours(3); mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
        %final pixel of the row
    elseif m==size_y
        
        city_b= [10+neighbours(2);14+neighbours(3); mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
        
        
        %all the inner pixels of the image
    else
        
        city_b= [14+ neighbours(1);10+neighbours(2);14+neighbours(3);...
            10+neighbours(4); mypixel];
        %we update the image-pixel with the smallest of the sums
        mynewpixel=min(city_b(:));
    end
    
end
end