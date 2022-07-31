
% read the receiver function dataset(or other seismic dataset)
% preprocess: transform the data from t-x domain to tau-p domain by using
% the high resolution Radon transform based on L1 norm
% To improve efficiency, the matlab parallel pool option has been turned on
% Written by Yifan Lu, more details can be found in https://doi.org/10.1093/gji/ggac260

clear;clc;close all;
set(0,'defaultfigurecolor','w');

% read the dataset including RFs baz deg ...
in_file = 'data_org.mat';
load(in_file);

% some parameters
shot = RFs(1:10:2400,:);
nx = size(shot,2);                                                              % Number of stations
nt = size(shot,1);                                                              % Sampling points
dx = 0.5;                                                                       % Spatial sampling rate£¬distance:km
dt = 0.005*10;                                                                  % time sampling rate£¬time:s
x = (0:nx-1)*dx;
t = ((0:nt-1)*dt)'-2;
xoff = x - 50;

% pad the data
pad_length = 15;
pad_npoint = round(pad_length / dt);
shot = [zeros(pad_npoint,nx);shot;zeros(pad_npoint,nx)];                         % RFs after pad
t = [(min(t) - (pad_npoint:-1:1)*dt)';t;(max(t) + (1:pad_npoint)*dt)'];          % time series after pad

%% Radon transform

% theoretical range£º Pmax = 1/(2*¡÷x*fmax)
% theoretical interval£º deltap = 1/((Xmax - Xmin)*fmax)

pmin = -1.5;
pmax =  1.5;
dp = 0.003;
p = pmin:dp:pmax;

maxiter = 20;                                                                     % number of iterations
rthresh = 0.001;                                                                  % threshold of iterations
method = 'cgg';
[ Rfft,f,iF,rfft,allniter,allfiter ] = Radon_forward(p,t,shot',xoff,maxiter,rthresh,method);
R = ifft(Rfft,size(Rfft,2),2,'symmetric');
disp('pre_process finished...');

save('pre_processing.mat');

