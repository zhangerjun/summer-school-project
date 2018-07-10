function [name,ext] = remove_ext_from_nifti(filename)
%remove the file extension (.nii or .nii.gz) from a nifti file 
% input
% filename - filename of the input nifti file
%
% output
% name -  nifti filename with extension removed
% ext -  file extension that has been removed
%
% Author
% Paddy Slator (p.slator@ucl.ac.uk)


niigz_pos=strfind(filename,'.nii.gz');
nii_pos=strfind(filename,'.nii');

if ~isempty(niigz_pos)
    name=filename(1:(niigz_pos-1));
    ext=filename(niigz_pos:end);
else
    name=filename(1:(nii_pos-1));
    ext=filename(nii_pos:end);
end



end