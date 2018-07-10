function normalised_dw_image = normalise_to_b0(dw_image,b)
%normalise diffusion weighted images so that the mean of the b=0 images is
%equal to 1.
%
%inputs 
%dw_image - the diffusion weighted image to normalise
%b - corresponding b-values
%
%output
%normalised_dw_image - the normalised diffusion weighted image
%
% Author
% Paddy Slator (p.slator@ucl.ac.uk)

%check that the number of b values matches the number of images 
if size(dw_image,4)~=length(b)
   disp('can''t normalise dw image: number of b-values doesn''t match number of volumes')
   normalised_dw_image=[];
   return
end

%add smallest number to the image to prevent divide by zero errors
dw_image = dw_image + eps;
%make sure it's a double
dw_image = double(dw_image);
%find b0 image (or images)
b0_index=find(b==0);

%if more than one b0 image, take the mean of them
if length(b0_index)>1
    b0_image=mean(dw_image(:,:,:,b0_index),4);
elseif length(b0_index)==1
    b0_image=dw_image(:,:,:,b0_index);
else
    disp('can"t normalise dw image: no b0 volumes')
    normalised_dw_image=[];
    return
end

normalised_dw_image=zeros(size(dw_image));

for i=1:length(b)
    normalised_dw_image(:,:,:,i)=dw_image(:,:,:,i)./b0_image;
end

   


end