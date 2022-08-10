%% PIV main
%
% Author: Riccardo Giordani
%
% Main for PIV analyzer
% Read two PIV images taken at different times and calculate the flow field

close all
clear 
clc; 

tic

%% LOAD PIV IMAGES

image1 = imread('NACA23012_1.bmp');
image2 = imread('NACA23012_2.bmp');

% Save a copy of the imported original images
image1_orig = image1;
image2_orig = image2;


%% APPLY MASK TO REMOVE AIRFOIL AREA 
image1 = rm_airfoil(image1);
image2 = rm_airfoil(image2);


%% PRE-PROCESS

% Filter execution order
order = [1,4,2,3];

% 1) Contrast Limit Adaptive Histogram Equalization 
% CLAHE parameters
clahe = 1; % 1=activated, 0=deactivated
clahe_window = 32;

% 2) High Pass Filter
% HPF parameters
hpf = 1;  % 1=activated, 0=deactivated
hpf_size = 15;

% 3) Min-max Filter
minmax = 0;  % 1=activated, 0=deactivated
minmax_window = 15; %15

% 4) Intensity campting
int_capt = 1;  % 1=activated, 0=deactivated
capt_scaling = 1.5;  % Scaling factor: 0.5 < n < 2

% Image enhancement
image1 = preproc(image1,order,clahe,clahe_window,hpf,hpf_size,minmax,minmax_window,int_capt,capt_scaling);
image2 = preproc(image2,order,clahe,clahe_window,hpf,hpf_size,minmax,minmax_window,int_capt,capt_scaling);


%% IMAGE ANALYSIS

% Interrogation window size [w_size x w_size]
w_size = 32;  % [pixel]

% Distance between interrogation window center points
% (overlapping (w_size-wc_step))
wc_step = 16; % [pixel]

% Overlapping percentage
overlap = (w_size - wc_step)/w_size * 100;
fprintf('Percentuale di sovrapposizione delle finestre: %d %% \n', overlap);

% Subpixel interpolation 3pt Gauss
subpx = 1;  %  % 1=activated, 0=deactivated

% Cross-correlation analysis
% [x, y, dx, dy] = CC (image1,image2,w_size, wc_step, subpx);
[x, y, dx, dy] = XC(image1,image2,w_size, wc_step,subpx);

% Graphical representation - Displacement field
figure(1)
imagesc(double(image1_orig)+double(image2_orig));
colormap('gray');
hold on
quiver(x,y,dx,dy,'r','AutoScaleFactor', 1.5);
hold off;
axis image;
title('Displacement Field [pixel]')
drawnow;


%% VELOCITY FIELD COMPUTATION

% Length scale
L = 103e-3; % [m]
H = 82e-3;  % [m]

yscale = L/size(image1,1);  % [m/px]
xscale = H/size(image1,2);  % [m/px]

% Temporal scale
Dt = 10e-6; % [sec]
dt = 1e-6;  % [sec]

%tscale = Dt + dt;                                                   
tscale = Dt;

% Conversion displacement -> velocity
u = dx.*(xscale/tscale);
v = dy.*(yscale/tscale);

% Velocity magnitude
V = sqrt(u.^2+v.^2);

% Scatter plot displacements
% dx_res = reshape(dx,[],1);
% dy_res = reshape(dy,[],1);
% 
% scatter(dy_res,dx_res,'+')

% Scatter plot velocity components
u_res = reshape(u,[],1);
v_res = reshape(v,[],1);

scatter(v_res,u_res,'+');
xlabel('v [m/s]');
ylabel('u [m/s]');

% Graphical representation - Velocity field
figure(2)
imagesc(double(image1_orig)+double(image2_orig));
colormap('gray');
hold on
quiver(x,y,u,v,'r','AutoScaleFactor',3);
hold off;
axis image;
title('Velocity Field [m/s]')
drawnow;

% parameter.u_lim(1) = input(' Limite inferiore velocità orizzontale: ');
% parameter.u_lim(2) = input(' Limite superiore velocità orizzontale: ');
% parameter.v_lim(1) = input(' Limite inferiore velocità verticale: ');
% parameter.v_lim(2) = input(' Limite superiore velocità verticale: ');

%% POST-PROCESS

% Velocity limit - 1st limit method 
method.choice = 0;   % (0=manual limits, 1=Vector difference test, 2=std dev test, 3=Normalized median test)
method.on_uv = 1;  % (0=test V, 1=test uv)

% Velocity limit - 2nd limit method
method2.active = 1;
method2.choice = 1; % (0=manual limits, 1=Vector difference test, 2=std dev test, 3=Normalized median test)
method2.on_uv = 1;  % (0=test V, 1=test uv)

% Velocity limit - 3rd limit method
method3.active = 1;
method3.choice = 3; % (0=manual limits, 1=Vector difference test, 2=std dev test, 3=Normalized median test)
method3.on_uv = 1;  % (0=test V, 1=test uv)

% Velocity test - parameters
parameter.tile_size = 3; % Working cell dimension

% Manual velocity limits - parameters
parameter.u_lim = [-33,17];
parameter.v_lim = [-52,18];
parameter.V_lim = [-50 50];

% Velocity difference limit - parameters
parameter.u_diff = 2;  
parameter.v_diff = 2;  
parameter.V_diff = 2;  

% Velocity normalized median limit - parameters
parameter.epsilon = 0.1; % offset (account for residual fluctuations obtained from correlation analysis)
parameter.res_toll = 2;  % Residual limit
parameter.builtin = 1;

% Missing velocity filled by interpolation
interp.uv = 1;  % (0=deactivated, 1=activated)
interp.V = 1;

% Outliers removal
[u,v,V] = postproc(u,v,V,method,interp,parameter);
if method2.active==1
    [u,v,V] = postproc(u,v,V,method2,interp,parameter);
end
if method3.active==1
    [u,v,V] = postproc(u,v,V,method3,interp,parameter);
end

% Graphical representation - Velocity field
figure(3)
imagesc(double(image1_orig)+double(image2_orig));
colormap('gray');
hold on
quiver(x,y,u,v,'r','AutoScaleFactor',2);
hold off;
axis image;
title('Velocity Field [m/s]')
drawnow;

figure(4)
hold on
contourf(x,y,V,':')
quiver(x,y,u,v,'k','AutoScaleFactor',1.5);
hold off
axis('ij')

figure(5)
hold on
contourf(x,y,V,':')
streamline = streamslice(x,y,u,v,7,'cubic');
set(streamline, 'Color', 'w' )
hold off
axis('ij')

% Print Velocity statistics
% fprintf('Average velocity = %f m/s \n', nanmean(V, 'all'));
% fprintf('Max velocity = %f m/s \n', nanmax(V, [], 'all'));

toc

