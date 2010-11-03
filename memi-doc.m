% Simulation of the MEMI system for partially incoherent light.
% 2010-10-31 Martin Kielhorn

%		            ky				                   -+-----
%    +--------------+--------------+	                |     \---
%    |           ---+---           |	                |         \-
%    |        --/   |   \--        |	            	|           \-
%    |      -/      |      \-      |	                |             \
%    |     /        |        \     |	                |              \
%    |     |        |        |     |	                |              |
%    |    /         |         \    |                    |      	       	\
%    +----+---------+---------+----+  kx                ----------------+--- kz
%    |    \         |         /    |                    |\-   alpha    	/
%    |     |        |        |     |	                |  \-          |
%    |     \        |        /     |	                |    \-        /
%    |      -\      |      /-      |	                |      \      /
%    |        --\   |   /--        |	                |       \-  /-
%    |           ---+---           |	        kxmax  -+---------X--
%    +--------------+--------------+	                |     /---  \-
%	     	                   	                    k  -+-----	      \
%    |------ s pixels -------------|		           	kx
%
%	       |- r pixels --------|			                       |-----|
%								                                     dkz
%


%  The objective is defined by its NA and its refractive index n, with
%  NA=n*sin(alpha). The sample distance in real space must be at least
%  dx=1/kxmax and dz=1/dkz if electric fields should be represented
%  without aliasing. We are interested in the intensity in real
%  space. Its calculation is a convolution in Fourierspace, making a
%  donut out of the sphere cap that is twice as wide and twice as
%  high. The sampling in real space should therefore be dx=1/(2*kxmax)
%  and dz=1/(2*dkz).

%  Looking at the right image we see: sin(alpha)=kxmax/k. Therefore we
%  have: kxmax=NA*k/n. We further know k=n/lambda, so we have:
%  kxmax=NA/lambda.

%  We can also express kz with alpha: cos(alpha)=kz/k. Using
%  cos(asin(x))=sqrt(1-x^2) we find:
%  dkz=k-k*cos(alpha)=(n-sqrt(n^2-NA^2))/lambda.

%  With NA=1.38, n=1.515 and lambda=.5 um we obtain: dx=181nm and
%  dz=281nm. Lets further set the size of the calculation domain to
%  s=256. Then the (rectangular) field is xmax=dx*s=46um wide in real
%  space. The diameter r of the back focal plane for this calculation
%  parameters would be 128 pixels.

%  To define all the parameters of the simulation, we start in the
%  BFP. We choose s, the width (and height) of the calculation
%  domain. On this array we determine kx and ky such, that kx=0 in the
%  center and kx=2*kxmax on the right border. For each of these points
%  we can calculate kz=sqrt(k-(kx^2+ky^2)). These choices inherently
%  define the width of the field in sample space. If we want to
%  simulate a bigger field of view we have to increase s. Note: one
%  could reduce kx on the border to kxmax (instead of 2*kxmax) but we
%  need some padding there. Otherwise artefacts due to aliasing will
%  deteriorate results.


%  Lets look at the simplest model of the MEMI system:

%	       	      	relais 	          tubelens
%	       	      	       	       	       	       |   	       sample
%      	  ^    	       	  ^    		      ^    	   v	  ^
% +-------+-   	  |   	  |    	  |    	  |    	         /+------+
% |	      | \--       	  |    	   	      |           /-- |
% |	      |    \- |   	  |	      |	      |         /-    |      |
% +-------+------\+-------+-------+--------------/-+------+----------> z
% |	      |    	  | \--   |	      |       |   /--      	  |	     |
% |	      |    	       \- |	              | /-         	  |
% |    	  |    	  |   	 \+-------|-------+-   	       	  |	     |
%      	  v    	       	  v          	  v    	   ^   	  v
% angle	     	 MMA  		 LCOS	      	       |    objective
%
%	       	      				  BFP
% 	       	      	          ^
% Lenses are drawn like this: |
%  		      	              v

% Front and back focal planes of neighboring lenses coincide.  Note:
% The focal lengths are different for each lens. In the simulation we
% account for that by changing the scale factor for each of the planes
% (i.e. what is the size of one pixel in the real world in this
% plane).

% Example: The focal length of the tubelens can be adjusted between
% 224mm to 448mm. To what size is the field of 46um magnified for a
% choice of 224mm and a 63x objective? The focal length of the
% objective with standard Zeiss tubelength fobj=164.5/63=2.6. The size
% of the 256x256 array in the LCOS plane is therefore:
% 0.046*224/2.6=3.9mm. Note: increasing s (from 256) will increase the
% size of the simulated LCOS.

% In this model the MMA and LCOS are assumed to be transmissive
% devices.

% To simulate the MEMI system for incoherent light we sum the
% intensity images of many coherent beams in sample space, i.e. we put
% a laser point source (think of mono mode fiber) plane
% 'angle'. Ideally we would do one simulation for each of the 256x256
% pixels in the calculation domain 'angle'. Probably we can limit
% ourselves to pixels that are image points of bright regions on the
% LCOS. Assuming the the MMA contains only one small 'hole' (bright
% region) - which is the standard operation condition of our device -
% we can reduce the number of pixels in 'angle' even more. For a small
% transmissive region of the MMA the image of the laser spot from
% angle will be blurred and illuminate several pixels on the LCOS. So
% we will not need to run simulations for neighboring pixels in
% 'angle'.

% Example: A gaussian with spot width = 8 pixels in the MMA plane,
% leads to a spot width of 10.2 pixels in the LCOS plane:

% s=[256 256];
% r=8;
% mma=exp(-rr(s).^2/r^2);
% lcos=exp(-rr(s).^2/(256/(pi*r))^2);
% lcos-real(ft(mma))/max(real(ft(mma))) % is practially zero

% For this example it would be sufficient to run 25x25 coherent
% propagations, if the full LCOS is transparent. In a real example the
% LCOS image will be sparse as well, reducing the necessary number of
% coherent propagations even more.

% Refocusing: When the following circular symmetric phase function is
% multiplied to the complex field distribution in the pupil plane,
% images at planes other than the focus can be calculated.
% exp(kz_of_kxky(:,:)*(1i*k/dx*Dz)); The defocus Dz is given in the
% same coordinates as dx. We would like to calculate the defocus up to
% 40 um deep, so that we can predict the light intensity into the
% deepest layer of the embryo. 


% The disk-gauss.jpg contains z slices at positions
% EDU>> double(z)
% ans =

%  Columns 1 through 8
%
%        0    0.2809    0.5619    0.8428    1.1238    1.4047    1.6857    1.9666
%
%  Columns 9 through 13
%
%    2.2476    2.5285    2.8095    3.0904    3.3714

% The maximum difference is 5.18% but one can see the diffraction
% pattern of the disk rasterizing the LCOS. Probably the sampling has
% to be increased. Assuming the pattern doesn't affect out-of-focus
% slices very much (good approximation as the pattern is out of focus)
% then we can take a look at the intensity distribution in those
% slices: The disk shows a bit more energy in higher angles.

% I conclude: gray values for the mma are not necessary.


% This is the code:

%%
% all lengths are in micrometer
NA=1.38; % numerical aperture
ri=1.515;% refractive index of immersion oil
lambda=.5; % emission vacuum wavelength
kxmax=NA/lambda; % radius of bfp for oil immersion objective
s=[256 256]; % image size
cosTheta=newim(s)+1;                  
kz_of_kxky=cosTheta;

klen=ri/lambda; % radius of ewald sphere

kx=xx(s,'freq')*4*kxmax;
ky=yy(s,'freq')*4*kxmax;              

tmp2=klen^2-(kx.^2+ky.^2);  
aperture=tmp2>0; % only the cap of the ewald sphere contributes
cosTheta(aperture)=sqrt(tmp2(aperture));        
                    % Theta, being the angle to the optic axis
kz_of_kxky=cosTheta;
%% prepare phase to propagate to different slices of the stack
nz=13; % number of slices
dz=lambda/(2*(ri-sqrt(ri^2-NA^2))); % step size in um
z=xx(nz,'corner')*dz; % defocus in um
sxy=1/(2*kxmax); % um per pixel
vol=[s nz];
fpropmat=newim(vol,'dcomplex');
pupil=newim(s,'dcomplex');
pupil_mask=berosion(aperture); % make aperture one pixel smaller, for division
for i=0:length(z)-1
    pupil=exp(kz_of_kxky(:,:)*(1i*klen/sxy*z(i)));
    pupil(pupil_mask)=pupil(pupil_mask)./sqrt(cosTheta(pupil_mask));
    fpropmat(:,:,i)=pupil;
end
clear pupil pupil_mask
%% spot width for lcos and mma
rl=.1; 
rm=.1;
%% define either smooth or b/w lcos image
%lcos=exp(-rr(s,'freq').^2/(rl/2)^2)
lcos=1.*(rr(s,'freq')<rl/2)
%% define a stack of mma images that get gradually smoother
%mma=exp(-((xx(s,'freq')-.2/2).^2 + (yy(s,'freq')-.3/2).^2)/(rm/2)^2);
ns=6; % number of smoothing steps
mma=newim([s ns],'dfloat');
% start with disk
mma(:,:,0)=((xx(s,'freq')-.2/2).^2 + (yy(s,'freq')-.3/2).^2)<(rm/4)^2;
for i=1:ns-1
    mma(:,:,i)=gaussf(squeeze(mma(:,:,0)),i);
end
% overlay mma image and bfp diameter for inspection
mma_overlay=mma+repmat(~aperture,[1 1 ns]);

%% spot size of mma in k-space, the smaller the mma, the less illumination
%% angles are needed for the incoherent image
rmk=2/(s(1)*pi*rm);
%abs(ft(mma))/max(abs(ft(mma)))-exp(-rr(s,'freq').^2/rmk^2);
rmk_pixels=rmk*256;
%% reduce sampling of lcos image according to spot size generated by mma
%% aperture, illuminate pixels that contribute with 3% to the result
small_lcos=resample(lcos,[1/rmk_pixels 1/rmk_pixels],[0 0]);
ill_mask=small_lcos>.03;
ill_mask=bdilation(ill_mask,2); % to be sure everything will be illuminated, grow mask
% kx must have kxmax on the border
ill_kx=xx(ill_mask,'freq')*rmk_pixels;
ill_ky=yy(ill_mask,'freq')*rmk_pixels;

%% compare mask and scaled lcos image
overlay(255*small_lcos,~ill_mask)
%% parallel coherent illumination
coh=newim([s ns nz],'dfloat');
% ft of each of the mma images with different smoothness
kmma=dip_fouriertransform(mma.*repmat(aperture,[1 1 ns]),'forward',[1 1 0]);
% multiply with lcos transmission, transform and filter in bfp
bfp=dip_fouriertransform(kmma.*repmat(lcos,[1 1 ns]),'forward',[1 1 0]).*repmat(aperture,[1 1 ns]);
for i=0:nz-1
    % defocused image in sample space for differently smoothed mma transmissions
    coh(:,:,:,i)=abs(dip_fouriertransform(bfp.*...
            repmat(squeeze(fpropmat(:,:,i)),[1 1 ns]),...
            'forward',[1 1 0])).^2;
end
clear kmma bfp

%% incoherent illumination
incoh=newim([s ns nz],'dfloat'); % sample space
lcos_incoh=newim([s ns],'dfloat'); % accumulates light in lcos plane
kmma=newim([s ns],'dcomplex'); % field on lcos
bfp=newim([s ns],'dcomplex'); % field in bfp
[row col]=size(ill_mask);
count=0;
all=sum(ill_mask);
for j=0:row-1
    for i=0:col-1
        if ill_mask(i,j)
            shifter=exp(1i*(ill_kx(i,j)*xx(s)+ill_ky(i,j)*yy(s)));
            kmma(:,:,:)=dip_fouriertransform(mma.*...
                repmat(aperture.*shifter,[1 1 ns]),'forward',[1 1 0]);
            lcos_incoh(:,:,:)=lcos_incoh(:,:,:)+abs(kmma(:,:,:)).^2;
            bfp(:,:,:)=dip_fouriertransform(kmma.*repmat(lcos,[1 1 ns]),...
                    'forward',[1 1 0]).*repmat(aperture,[1 1 ns]);
            for k=0:nz-1
                incoh(:,:,:,k)=squeeze(incoh(:,:,:,k))+...
                        abs(dip_fouriertransform(bfp.*...
                        repmat(squeeze(fpropmat(:,:,k)),[1 1 ns]),...
                        'forward',[1 1 0])).^2;
            end
            count=count+1;
            [count all]
        end
    end
end
clear bfp kmma shifter;

%%
permute(reshape(permute(incoh),[s(1) s(2) nz*ns]),[2 1 3])


%% For the next to run you need to do several simulations with different MMA
%% and LCOS images.
%% compare gauss and sharp edged mma hole
comp=reshape(gincoh/max(gincoh(:,:,0))-dincoh/max(dincoh(:,:,0)),[256 256*nz]);
out=255*(comp-min(comp))/(max(comp)-min(comp));
writeim(out,'/home/martin/1031/memi/disk-gauss-big.jpg','JPEG');
%%
mosaic=[reshape(gincoh,[256 256*nz]) reshape(dincoh,[256 256*nz])];
out=255*(mosaic-min(mosaic))/(max(mosaic)-min(mosaic));
writeim(out,...
    '/home/martin/1031/memi/gauss_disk.jpg'...
    ,'JPEG');
%% gaussian has only 60% of the intensity compared to disk
% with more sampling even only 56%
sum(gincoh(:,:,0))/sum(dincoh(:,:,0))