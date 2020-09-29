clear all
close all
clc; 


%% CARICAMENTO DELLA COPPIA DI IMMAGINI

image1 = imread('NACA23012_1.bmp');
image2 = imread('NACA23012_2.bmp');

image1_old = image1;

figure(1)
imshow(image1)

%% PRE-PROCESS
% Ordine esecuzione filtri
order = [1,2,4,3];

% 1) Contrast Limit Adaptive Histogram Equalization 
% Parametri CLAHE
clahe = 1; % 1=attivato, 0=disattivato
clahe_window = 32;

% 2) High Pass Filter
% Parametri HPF
hpf = 0;  % 1=attivato, 0=disattivato
hpf_size = 5;

% 3) Min-max Filter
minmax = 1;  % 1=attivato, 0=disattivato
minmax_window = 15; %15

% 4) Intensity campting
int_capt = 1;  % 1=attivato, 0=disattivato
capt_scaling = 1.5;  % Scaling factor: 0.5 < n < 2

% 5) Low Pass filter

% Miglioramento delle due immagini
image1 = preproc(image1,order,clahe,clahe_window,hpf,hpf_size,minmax,minmax_window,int_capt,capt_scaling);
image2 = preproc(image2,order,clahe,clahe_window,hpf,hpf_size,minmax,minmax_window,int_capt,capt_scaling);

figure(2)
imshow(image1)

figure(3)
imhist(image1_old)

figure(4)
imhist(image1)