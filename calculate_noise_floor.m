function noise_floor = calculate_noise_floor(SNR)
%calculate the mean of the magnitude reconstructed MR signal in the absense
%of any true signal.

%e.g. Jones, D. K., & Basser, P. J. (2004). ?Squashing peanuts and smashing pumpkins?: How noise distorts diffusion-weighted MR data. Magnetic Resonance in Medicine, 52(5), 979?993. http://doi.org/10.1002/mrm.20283

%Author:
%Paddy Slator (p.slator@ucl.ac.uk)


sigma=1/SNR;
noise_floor=sigma*sqrt(pi/2);





end