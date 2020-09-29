%% PIV main

close all
clear 
clc; 

tic

%% CARICAMENTO DELLA COPPIA DI IMMAGINI

image1 = imread('NACA23012_1.bmp');
image2 = imread('NACA23012_2.bmp');

% Salva le immagini originali
image1_orig = image1;
image2_orig = image2;


%% RIMOZIONE PROFILO NACA
image1 = rm_airfoil(image1);
image2 = rm_airfoil(image2);


%% PRE-PROCESS

% Ordine esecuzione filtri
order = [1,4,2,3];

% 1) Contrast Limit Adaptive Histogram Equalization 
% Parametri CLAHE
clahe = 1; % 1=attivato, 0=disattivato
clahe_window = 32;

% 2) High Pass Filter
% Parametri HPF
hpf = 1;  % 1=attivato, 0=disattivato
hpf_size = 15;

% 3) Min-max Filter
minmax = 0;  % 1=attivato, 0=disattivato
minmax_window = 15; %15

% 4) Intensity campting
int_capt = 1;  % 1=attivato, 0=disattivato
capt_scaling = 1.5;  % Scaling factor: 0.5 < n < 2

% Miglioramento delle due immagini
image1 = preproc(image1,order,clahe,clahe_window,hpf,hpf_size,minmax,minmax_window,int_capt,capt_scaling);
image2 = preproc(image2,order,clahe,clahe_window,hpf,hpf_size,minmax,minmax_window,int_capt,capt_scaling);


%% ANALISI DELLE IMMAGINI

% Dimenzione finestra di interrogazione [w_size x w_size]
w_size = 32;  % [pixel]

% Distanza centri finestre di interrogazione 
% (sovrapposizione di (w_size-wc_step))
wc_step = 16; % [pixel]

% Percentuale di sovrapposizione
overlap = (w_size - wc_step)/w_size * 100;
fprintf('Percentuale di sovrapposizione delle finestre: %d %% \n', overlap);

% Subpixel interpolation 3pt Gauss
subpx = 1;  %  (1=attivato, 0=disattivato)

% Analisi cross-correlazioni
% [x, y, dx, dy] = CC (image1,image2,w_size, wc_step, subpx);
[x, y, dx, dy] = XC(image1,image2,w_size, wc_step,subpx);

% Rappresentazione grafica del campo degli spostamenti vettoriali
figure(1)
imagesc(double(image1_orig)+double(image2_orig));
colormap('gray');
hold on
quiver(x,y,dx,dy,'r','AutoScaleFactor', 1.5);
hold off;
axis image;
title('Campo di velocità [pixel]')
drawnow;


%% DETERMINAZIONE DEL CAMPO DI VELOCITA'

% Scala di lunghezza
L = 103e-3; % [m]
H = 82e-3;  % [m]

yscale = L/size(image1,1);  % [m/px]
xscale = H/size(image1,2);  % [m/px]

% Scala temporale
Dt = 10e-6; % [sec]
dt = 1e-6;  % [sec]

%tscale = Dt + dt;                                                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tscale = Dt;

% Conversione spostamenti -> velocità
u = dx.*(xscale/tscale);
v = dy.*(yscale/tscale);

% Calcolo modulo della velocità
V = sqrt(u.^2+v.^2);

% Scatter plot spostamenti
% dx_res = reshape(dx,[],1);
% dy_res = reshape(dy,[],1);
% 
% scatter(dy_res,dx_res,'+')

% Scatter plot componenti velocità
u_res = reshape(u,[],1);
v_res = reshape(v,[],1);

scatter(v_res,u_res,'+');
xlabel('v [m/s]');
ylabel('u [m/s]');

% Rappresentazione grafica del campo di velocità
figure(2)
imagesc(double(image1_orig)+double(image2_orig));
colormap('gray');
hold on
quiver(x,y,u,v,'r','AutoScaleFactor',3);
hold off;
axis image;
title('Campo di velocità [m/s]')
drawnow;

% parameter.u_lim(1) = input(' Limite inferiore velocità orizzontale: ');
% parameter.u_lim(2) = input(' Limite superiore velocità orizzontale: ');
% parameter.v_lim(1) = input(' Limite inferiore velocità verticale: ');
% parameter.v_lim(2) = input(' Limite superiore velocità verticale: ');

%% POST-PROCESS

% Metodo scelta dei limiti sulle velocità
method.choice = 0;   % (0=limiti manuali, 1=Vector difference test, 2=std dev test, 3=Normalized median test)
method.on_uv = 1;  % (0=test su V, 1=test su uv)

% aggiunta sencondo metodo di selezione
method2.active = 1;
method2.choice = 1; % (0=limiti manuali, 1=Vector difference test, 2=std dev test, 3=Normalized median test)
method2.on_uv = 1;  % (0=test su V, 1=test su uv)

% aggiunta terzo metodo di selezione
method3.active = 1;
method3.choice = 3; % (0=limiti manuali, 1=Vector difference test, 2=std dev test, 3=Normalized median test)
method3.on_uv = 1;  % (0=test su V, 1=test su uv)

% Impostazione parametri per i test sulla velocità
parameter.tile_size = 3; % Dimensione della cella su operare

% Limiti di velocità manuali
parameter.u_lim = [-33,17];
parameter.v_lim = [-52,18];
parameter.V_lim = [-50 50];

% Limite per test sulla differenza di velocità
parameter.u_diff = 2;  
parameter.v_diff = 2;  
parameter.V_diff = 2;  

% Parametri test normalizzato della mediana
parameter.epsilon = 0.1; % offset che tiene conto delle fluttuazioni rimanenti ottenute dall'analisi delle correlazioni
parameter.res_toll = 2;  % Limite per il residuo
parameter.builtin = 1;

% Interpolazione velocità mancanti
interp.uv = 1;  % (0=disattivato, 1=attivato)
interp.V = 1;

% Rimozione degli outliers e interpolazione
[u,v,V] = postproc(u,v,V,method,interp,parameter);
if method2.active==1
    [u,v,V] = postproc(u,v,V,method2,interp,parameter);
end
if method3.active==1
    [u,v,V] = postproc(u,v,V,method3,interp,parameter);
end

% Rappresentazione grafica del campo di velocità
figure(3)
imagesc(double(image1_orig)+double(image2_orig));
colormap('gray');
hold on
quiver(x,y,u,v,'r','AutoScaleFactor',2);
hold off;
axis image;
title('Campo di velocità [m/s]')
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

% Stampa valori delle velocità
% fprintf('Velocità media = %f m/s \n', nanmean(V, 'all'));
% fprintf('Velocità massima = %f m/s \n', nanmax(V, [], 'all'));

toc

