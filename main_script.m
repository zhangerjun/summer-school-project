%Microstructure model selection in fetal diffusion MRI

%This is the main file for you to run. There are also a number of sub-functions which you will need to write.

%First of all, make sure the 'summer-school-project' directory and sub-directories are on the path.

%PART 1

clear 

%load the dmri image "dmri_2909.nii.gz" from the data folder using the "load_untouch_nii" function
dmri = load_untouch_nii('dmri_2909.nii.gz');
%use imagesc to view the first diffusion-weighted volume of the 7th z-slice of the image 
figure;imagesc(dmri.img(:,:,7,1))
%convert the image to double - this helps with lots of things!
dmri.img = double(dmri.img);
%zero measurements are bad news! so add the smallest floating point number (eps) to the image
dmri.img  = dmri.img + eps;

%find all mask files
mask_filenames = dir('data/*mask*');
n_masks = length(mask_filenames);

%load the masks in a structure called masks - with fields given by the names in "mask_filenames"
for i=1:n_masks
	this_mask = load_untouch_nii(mask_filenames(i).name);

	mask_names{i} = remove_ext_from_nifti(mask_filenames(i).name);

	masks.(mask_names{i}) = this_mask.img;
end

%load the dmri gradient table "dmri_2909_grads.txt" using importdata
grads = importdata('dmri_2909_grads.txt');
%extract the b-value and b-vector
bvecs = grads(:,1:3);
bvals = grads(:,4);


%normalise the image to the b=0 volumes
normdmri = dmri;
normdmri.img = normalise_to_b0(dmri.img,bvals);


%calculate the mean normalised signal (use calculate_mean_signal) in each of the masks
for i=1:n_masks
	mean_norm_sig.(mask_names{i}) = calculate_mean_signal(normdmri.img,masks.(mask_names{i}));
end


%plot the mean signal against b-values within each mask
figure;hold on;
for i=1:n_masks
	plot(bvals, mean_norm_sig.(mask_names{i}),'o')
end
legend(make_nice_figure_string(mask_names))
xlabel('b-values s/mm^2')
ylabel('signal')





%Now we will fit some models to the mean signal

%first we choose the matlab nonlinear fitting function options
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');

%next estimate the signal to noise ratio for each mask by taking the standard deviation of the b=0 volumes
for i=1:length(mask_names)
	this_mask_mean_sig = mean_norm_sig.(mask_names{i});
	SNR.(mask_names{i}) = 1/std(this_mask_mean_sig(bvals==0));
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


%function which calculates the sum of the residuals - we use an inline function for readability and simplicity 
%otherwise would need to pass extra parameters (i.e. meas,params,grads) to fmincon,
%(see - http://uk.mathworks.com/help/optim/ug/passing-extra-parameters.html if you really want to do this)

%list of models to fit 
models = {'ADC','IVIM','triexp','quadexp'};
%TASK - fill in the missing "synth_" files. 
%These take the model params and gradient table as input and return sythetic dmri signal. 
%I have given you synth_triexp as an example.


%now we do the fitting
for i=1:length(mask_names)%loop over ROI masks
	%get the dmri measurement to fit to - this ROI's normalised mean signal 
	meas = mean_norm_sig.(mask_names{i})';
	%get the sigma (standard deviation of Rician noise distribution) for this mask
	sigma = 1/SNR.(mask_names{i});
	for j=1:length(models)%loop over models	

		%get this model's synthetic measurement function for passing to sumres
		synthfun = str2func(['synth_' models{j}]);			

		%define inline function for calculating sum of residuals for this model  
		sumres = @(params) -RicianLogLik(meas, synthfun(params,grads), sigma);
	 	
		%minimize the sumres function for this model, subject to initvals, minvals and maxvals
		[fitted_params.(mask_names{i}).(models{j}),fval.(mask_names{i}).(models{j})] = fmincon(sumres,initvals.(models{j}),[],[],[],[],minvals.(models{j}),maxvals.(models{j}),[],[]);

		%TASK calculate the BIC for this model 
		loglik = -fval.(mask_names{i}).(models{j});
		nparams = length(fitted_params.(mask_names{i}).(models{j}));
		ndata = length(meas);
		BIC.(mask_names{i}).(models{j}) = nparams*log(ndata) - 2*loglik;
	end
end


%plot the measured against synthetic signal for each mask - to see if the fits worked ok (make sure they look reasonable!)
for i=1:length(mask_names)
	figure;hold on;
	plot(bvals, mean_norm_sig.(mask_names{i}),'o')
	plot(bvals,synth_ADC(fitted_params.(mask_names{i}).ADC,grads) + fitted_params.(mask_names{i}).ADC(1)*calculate_noise_floor(SNR.(mask_names{i})),'x')
	plot(bvals,synth_IVIM(fitted_params.(mask_names{i}).IVIM,grads) + fitted_params.(mask_names{i}).IVIM(1)*calculate_noise_floor(SNR.(mask_names{i})),'^')
	plot(bvals,synth_triexp(fitted_params.(mask_names{i}).triexp,grads) + fitted_params.(mask_names{i}).triexp(1)*calculate_noise_floor(SNR.(mask_names{i})),'sq')
	legend({'signal','ADC','IVIM','triexp'})
	title(make_nice_figure_string(mask_names{i}))
end


%now have a look at the BIC values for each model in each ROI (lowest BIC indicates the model that explains the data best)
%Do the results make sense?

%put the BIC values in a matrix so it's easy to do a barplot
BIC_matrix = zeros(length(mask_names),length(models));
for i=1:length(mask_names)
	for j=1:length(models)
		BIC_matrix(i,j) = BIC.(mask_names{i}).(models{j});
	end
end
%TASK - make a barplot of the BIC values for all models - grouped by ROI masks






%PART 2 
%% Now we'll fit some models voxelwise to the data, starting with IVIM
% (I've given you skeleton code for a voxelwise IVIM fit - if you have time you 
% can fit more models)
% Warning - this will take a few minutes to run!


%for voxelwise fits we'll just define the signal to noise ratio
SNR = 20;
sigma = 1/SNR;


%get the dimensions of the image
Nx = size(dmri.img,1);
Ny = size(dmri.img,2);
Nz = size(dmri.img,3);


%we'll store the parameters and BIC values in arrays so we can plot them later 
nparam_IVIM = 4;
%TASK - initialise these arrays as zeros -the BIC arrays should be Nx x Ny x Nz, 
%the parameter arrays Nx x Ny x Nz x number of model parameters
param_maps.IVIM = zeros([size(masks.all_masks) nparam_IVIM]);
BIC_map.IVIM = zeros(size(masks.all_masks));



%(add more models here if you have time!)
models = {'IVIM'};


for j=1:length(models)

	%set up a counter
	l=1;

	%takes a long time so we'll only fit a single slice 
	zslice = 7;
	%work out the number of voxels to fit
	nvox = nnz(masks.all_masks(:,:,zslice))

	for x = 1:Nx
		for y = 1:Ny
			for z = zslice
				if masks.placenta_and_uterine_wall_mask(x,y,z)
					%get the measurement for this voxel
					meas = double(squeeze(normdmri.img(x,y,z,:)));

					%get this model's synthetic measurement function for passing to sumres
					synthfun = str2func(['synth_' models{j}]);	

					%define inline function for calculating sum of residuals for this model  
					sumres = @(params) -RicianLogLik(meas, synthfun(params,grads), sigma);


					[fitted_params_vox,fval_vox,exitflag] = fmincon(sumres,initvals.(models{j}),[],[],[],[],minvals.(models{j}),maxvals.(models{j}),[]);

					%store the parameter values	
					param_maps.(models{j})(x,y,z,:) = fitted_params_vox;

					%calculate the BIC for this model in this voxel
					loglik = -fval_vox;
					nparams = length(fitted_params_vox);
					ndata = length(meas);
					BIC_vox = nparams*log(ndata) - 2*loglik;

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



%harder stuff! Fit some anisotropic models. These are trickier to fit.
%We'll start with Stick-ball, which is like IVIM, but models the fast-attenuating signal
%with a "stick" compartment.
% Again I'll give you skeleton code for just one model - if you have time you 
% can fit more 

%initialise the parameter and BIC map arrays
nparam_StickBall = 6;
param_maps.StickBall = zeros([size(masks.all_masks) nparam_StickBall]);
BIC_map.StickBall = zeros(size(masks.all_masks));


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
	%work out the number of voxels to fit
	nvox = nnz(masks.all_masks(:,:,zslice))

	for x = 1:Nx
		for y = 1:Ny
			for z = zslice
				if masks.placenta_and_uterine_wall_mask(x,y,z)
					%get the measurement for this voxel
					meas = double(squeeze(normdmri.img(x,y,z,:)));

					%get this model's synthetic measurement function for passing to sumres
					synthfun = str2func(['synth_' models{j}]);	

					%define inline function for calculating sum of residuals for this model  
					sumres = @(params) -RicianLogLik(meas, synthfun(params,grads), sigma);
										
					%do the fit
					[fitted_params_vox,fval_vox,exitflag] = fmincon(sumres,initvals.(models{j}),[],[],[],[],minvals.(models{j}),maxvals.(models{j}),[]);

					%store the parameter values	
					param_maps.(models{j})(x,y,z,:) = fitted_params_vox;

					%calculate the BIC for this model in this voxel
					loglik = -fval_vox;
					nparams = length(fitted_params_vox);
					ndata = length(meas);
					BIC_vox = nparams*log(ndata) - 2*loglik;

					%store the BIC value	
					BIC_map.(models{j})(x,y,z) = BIC_vox;

					disp(['fit voxel' num2str(l) ' of ' num2str(nvox)])
					l=l+1;
				end						
			end
		end
	end
end






%stuff that I don't think I'll use

% models={'DT'};
% j=1;

% %set up a counter
% l=1;

% %takes a long time so we'll only fit a single slice 
% zslice = 7;
% %work out the number of voxels to fit
% nvox = nnz(masks.placenta_and_uterine_wall_mask(:,:,zslice))

% nparam_DT = 7;
% param_maps.DT = zeros([size(masks.placenta_and_uterine_wall_mask) nparam_DT]);
% initTheta = zeros(size(masks.placenta_and_uterine_wall_mask));
% initPhi = zeros(size(masks.placenta_and_uterine_wall_mask));


% %convert the gradient table to protocol object
% protocol = grad2protocol(grads);

% for x = 1:Nx
% 	for y = 1:Ny
% 		for z = zslice
% 			if masks.placenta_and_uterine_wall_mask(x,y,z)
% 				%get the measurement for this voxel
% 				meas = double(squeeze(normdmri.img(x,y,z,:)));

% 				%Fit the linear DT model 
% 				fitted_params_vox = FitLinearDT(meas,protocol,1);
				
% 				%store the DT fit parameter values	
% 				param_maps.(models{j})(x,y,z,:) = fitted_params_vox;

% 				%put into 3x3 tensor form
% 				dt = D_vector_to_D(fitted_params_vox(2:end));
% 				%get the eigenvalues
% 				[eigvec, eigval] = eig(dt);
%     			%find the principle (highest) eigenvalue
%     			[~, ind] = sort(abs([eigval(1, 1) eigval(2, 2) eigval(3, 3)]));
%     			%fibre direction is the eigenvector corresponding to the highest eigenvalue    
%     			fibredir = eigvec(:,ind(3));
%     			%convert these to polar and azimuthal angles
% 				Theta = acos(fibredir(3));
%     			Phi = atan2(fibredir(2), fibredir(1));
%     			%store these angles - will use as initial values for the fitting
%     			initTheta(x,y,z) = Theta;
% 				initPhi(x,y,z) = Phi;
				

% 				%calculate the BIC for this model in this voxel
% 				loglik = -fval_vox;
% 				nparams = length(fitted_params_vox);
% 				ndata = length(meas);
% 				BIC_vox = nparams*log(ndata) - 2*loglik;

% 				%store the BIC value	
% 				BIC_map.(models{j})(x,y,z) = BIC_vox;

% 				disp(['fit voxel' num2str(l) ' of ' num2str(nvox)])
% 				l=l+1;
% 			end
					 
% 		end
% 	end
% end


