%Microstructure model selection in fetal diffusion MRI
%This is the main file for you to run. There are also a number of sub-functions which you will need to write.
%First of all, make sure the 'summer-school-project' directory and sub-directories are on the path.

%% PART 1

clear 

%TASK - load the dmri image "dmri_2909.nii.gz" from the data folder using the "load_untouch_nii" function
%assign it to a variable called "dmri" - this will be a structure containing
%the the image in the "img" field.
dmri = 0;
%TASK - use imagesc to view the first diffusion-weighted volume of the 7th z-slice of the image 

%TASK - convert the image to double - this helps with lots of things!

%TASK - zero measurements are bad news! so add the smallest floating point number (eps) to the image


%find all mask files
mask_filenames = dir('data/*mask*');
n_masks = length(mask_filenames);

%load the masks in a structure called masks - with fields given by the names in "mask_filenames"
for i=1:n_masks %loop over the mask filenames
	%TASK - load the current mask to a variable called "this_mask" 
	this_mask = 0;
	%we use a utilities function to make a nice name for the mask 
	mask_names{i} = remove_ext_from_nifti(mask_filenames(i).name);
	%TASK - store this_mask in a field called "mask_names{i}" in a structure called masks
	masks.(mask_names{i}) = 0;
end

%TASK - load the dmri gradient table "dmri_2909_grads.txt" using importdata
grads = 0;
%TASK - extract the b-value and b-vector (i.e. gradient directions) from grads 
bvecs = 0;
bvals = 0;

%now we will normalise the image
%use dmri image as a template
normdmri = dmri;
%TASK - normalise the image to the b=0 volumes using the utilities file "normalise_to_b0"



%TASK - calculate the mean normalised signal (use calculate_mean_signal) in each of the masks
for i=1:n_masks
	mean_norm_sig.(mask_names{i}) = 0;
end


%TASK - plot the mean signal (y-axis) against b-values (x-axis) within each mask
figure;hold on;
for i=1:n_masks
	%plot here	
end
legend(make_nice_figure_string(mask_names))
xlabel('b-values s/mm^2')
ylabel('signal')





%Now we will fit some models to the mean signal


%TASK - estimate the signal to noise ratio for each mask by taking the 1 over the standard deviation of the b=0 volumes
for i=1:length(mask_names)
	SNR.(mask_names{i}) = 0;
end

%choose the initial, minimum, and maximum values for each model
%ADC params - S0, ADC
initvals.ADC = [1  0.005];
minvals.ADC = [0 0];
maxvals.ADC = [10 1];

%IVIM params - S0, f, dp, d
initvals.IVIM = [1 0.5 0.05 0.002];
minvals.IVIM = [0 0 0 0];
maxvals.IVIM = [10 1 1 .005];

%triexp params - S0, f1, f2, d1, d2, d3
initvals.triexp = [1 0.3 0.3 0.05 0.01 0.002];
minvals.triexp = [0 0 0 0 0 0];
maxvals.triexp = [10 1 1 1 .05 .005];

%quadexp params - S0, f1, f2, f3, d1, d2, d3, d4
initvals.quadexp = [1 0.3 0.3 0.3 0.1 0.05 0.01 0.002];
minvals.quadexp = [0 0 0 0 0 0 0 0];
maxvals.quadexp = [10 1 1 1 1 1 .05 .005];



%list of models to fit 
models = {'ADC','IVIM','triexp','quadexp'};
%TASK - fill in the missing "synth_" files. 
%These take the model params and gradient table as input and return sythetic dmri signal. 
%I have given you synth_triexp as an example.


%now we do the fitting
for i=1:length(mask_names)%loop over ROI masks
	%TASK - get the dmri measurement to fit to - this ROI's normalised mean signal 
	meas = 0;
	%get the sigma (standard deviation of Rician noise distribution) for this mask
	sigma = 1/SNR.(mask_names{i});
	for j=1:length(models)%loop over models	

		%get this model's synthetic measurement function for passing to sumres
		synthfun = str2func(['synth_' models{j}]);			

		%define inline function for calculating sum of residuals for this model - we use an inline function for readability and simplicity 
		%otherwise we would need to pass extra parameters (i.e. meas,params,grads) to fmincon,
		%(see - http://uk.mathworks.com/help/optim/ug/passing-extra-parameters.html if you really want to do this)
		%TASK - fill in RicianLogLik
		sumres = @(params) -RicianLogLik(meas, synthfun(params,grads), sigma);
	 	
		%TASK - minimize the sumres function for this model, subject to initvals, minvals and maxvals, using the fmincon function
		%[fitted_params.(mask_names{i}).(models{j}),fval.(mask_names{i}).(models{j})] = fmincon();

		%TASK - calculate the BIC for this model 
		loglik = 0;
		nparams = 0;
		ndata = 0;
		BIC.(mask_names{i}).(models{j}) = 0;

	end
end


%now we plot the measured against synthetic signal for each mask - to see if the fits worked ok (make sure they look reasonable!)
%note that we add the Rician noise offset to the measurements - using the calculate_noise_floor function
for i=1:length(mask_names)
	figure;hold on;
	%plot the data
	plot(bvals, mean_norm_sig.(mask_names{i}),'o')
	%uncomment to plot the ADC model fit
	%plot(bvals,synth_ADC(fitted_params.(mask_names{i}).ADC,grads) + fitted_params.(mask_names{i}).ADC(1)*calculate_noise_floor(SNR.(mask_names{i})),'x')
	%TASK - add the plots for the other models  



	legend({'signal','ADC'})
	title(make_nice_figure_string(mask_names{i}))
end


%TASK - have a look at the BIC values for each model in each ROI (lowest BIC indicates the model that explains the data best)
%Do the results make sense?

%put the BIC values in a matrix so it's easy to do a barplot
BIC_matrix = zeros(length(mask_names),length(models));
for i=1:length(mask_names)
	for j=1:length(models)
		BIC_matrix(i,j) = BIC.(mask_names{i}).(models{j});
	end
end
%TASK - make a barplot of the BIC values for all models - grouped by ROI masks






 
%% PART 2 - Now we'll fit some models voxelwise to the data, starting with IVIM
% (I've given you skeleton code for a voxelwise IVIM fit - if you have time you 
% can fit more models)
% Warning - this will take a few minutes to run!


%for voxelwise fits we'll just define the signal to noise ratio
SNR = 20;
sigma = 1/SNR;


%TASK - get the dimensions of the image
Nx = 0;
Ny = 0;
Nz = 0;


%we'll store the parameters and BIC values in arrays so we can plot them later 
%define the number of parameters in the IVIM model
nparam_IVIM = 4;
%TASK - initialise these arrays as zero arrays -the BIC arrays should be Nx x Ny x Nz, 
%the parameter arrays Nx x Ny x Nz x number of model parameters
param_maps.IVIM = 0;
BIC_map.IVIM = 0;



%you can add more models here later if you have time!
models = {'IVIM'};

%do the voxelwise fitting for each model 
for j=1:length(models)
	%set up a counter
	l=1;

	%takes a long time so we'll only fit a single slice 
	zslice = 7;
	%TASK - work out the number of voxels we need to fit using nnz
	nvox = 0;

	for x = 1:Nx
		for y = 1:Ny
			for z = zslice
				if masks.placenta_and_uterine_wall_mask(x,y,z) %fit for the placenta and uterine wall mask
					%TASK - get the dmri measurement for this voxel
					meas = 0;

					%get this model's synthetic measurement function for passing to sumres
					synthfun = str2func(['synth_' models{j}]);	

					%define inline function for calculating sum of residuals for this model  
					sumres = @(params) -RicianLogLik(meas, synthfun(params,grads), sigma);

					%TASK - minimize the sumres function for this model in this voxel, subject to initvals, minvals and maxvals, using the fmincon function
					%[fitted_params_vox,fval_vox,exitflag] = fmincon();

					%store the parameter values
					param_maps.(models{j})(x,y,z,:) = fitted_params_vox;

					%TASK - calculate the BIC for this model in this voxel
					loglik = 0;
					nparams = 0;
					ndata = 0;
					BIC_vox = 0;

					%store the BIC value	
					BIC_map.(models{j})(x,y,z) = BIC_vox;

					disp(['fit voxel' num2str(l) ' of ' num2str(nvox)])
					l=l+1;
				end
						 
			end
		end
	end
end


%TASK - use "plot_overlayed_images" to make a nice picture!





%% harder stuff! Fit some anisotropic models. These are trickier to fit.
%We'll start with Stick-ball, which is like IVIM, but models the fast-attenuating signal
%with a "stick" compartment.
% Again I'll give you skeleton code for just one model - if you have time you 
% can fit more.

%define the number of parameters in the stick-ball model
nparam_StickBall = 6;
%TASK - initialise the parameter and BIC map arrays
param_maps.StickBall = 0;
BIC_map.StickBall = 0;


%choose the initial, minimum, and maximum values for stick-ball models
%stick ball params - S0, f, dp, d, theta, pi
initvals.StickBall = [1 0.5 0.05 0.002 0 0];
minvals.StickBall = [0 0 0 0 -100 -100];
maxvals.StickBall = [10 1 1 .01 100 100];


models = {'StickBall'};

for j=1:length(models)
	%set up a counter
	l=1;

	%takes a long time so we'll only fit a single slice 
	zslice = 7;
	%TASK - work out the number of voxels to fit using nnz
	nvox = 0;

	for x = 1:Nx
		for y = 1:Ny
			for z = zslice
				if masks.placenta_and_uterine_wall_mask(x,y,z)
					%TASK - get the measurement for this voxel
					meas = 0;

					%get this model's synthetic measurement function for passing to sumres
					synthfun = str2func(['synth_' models{j}]);	

					%define inline function for calculating sum of residuals for this model  
					sumres = @(params) -RicianLogLik(meas, synthfun(params,grads), sigma);
										
					%TASK - minimize the sumres function for this model in this voxel, subject to initvals, minvals and maxvals, using the fmincon function
					%[fitted_params_vox,fval_vox,exitflag] = fmincon();

					%store the parameter values	
					param_maps.(models{j})(x,y,z,:) = fitted_params_vox;

					%TASK - calculate the BIC for this model in this voxel
					loglik = 0;
					nparams = 0;
					ndata = 0;
					BIC_vox = 0;

					%store the BIC value	
					BIC_map.(models{j})(x,y,z) = BIC_vox;

					disp(['fit voxel' num2str(l) ' of ' num2str(nvox)])
					l=l+1;
				end						
			end
		end
	end
end






