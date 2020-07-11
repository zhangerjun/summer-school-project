%Microstructure model selection in fetal diffusion MRI

%This is the main file for you to run. There are also a number of sub-functions which you will need to write.

%First of all, make sure the 'summer-school-project' directory and sub-directories are on the path.

%PART 1

clear 

%load the dmri image "dmri_2909.nii.gz" from the data folder using the "load_untouch_nii" function
dmri = load_untouch_nii('data/dmri_2909.nii.gz');
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
%make a barplot of the BIC values for all models - grouped by ROI masks
figure;
bar(categorical(make_nice_figure_string(mask_names)),BIC_matrix)
legend(models)






