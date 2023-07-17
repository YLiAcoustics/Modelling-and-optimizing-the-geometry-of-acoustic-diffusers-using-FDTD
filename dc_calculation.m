%
% Program info
% Program author: Yuqing Li, August 2019
% Program details: This matlab script calls the function FDTD_func.m
% to measure the polar response for a diffuser with a known well depth sequence. 
% It calculate the diffusion coefficient and generate a polar plot for a known well depth sequence. 
% The user can test the performance of any type of diffuser over the bandwidth of interest.

close all
clearvars

%%%% EDIT %%%%
%% global options                      
sys_opts.freq_range = '1/3 octave';                     % bandwidth of diffusion coefficient measurement: choose between '1/3 octave' and 'wideband'
sys_opts.polarplot = 1;                                % plot polar response: choose between 1(yes) and 0(no)

%% FDTD options (to be sent to FDTD function)
% simulation options                   
sim_opts.use_gpuarray = 0;                              % use GPU calculation
sim_opts.plot_on = 0;                                   % plot FDTD simulation
sim_opts.calculate_efficiency = 1;                      % calculate simulation efficiency
sim_opts.spatial_operater = 'centered';                 % spatial difference operator for boundary update: choose between 'centered' and 'noncentered'
sim_opts.source_type = 'sin';                           % sound source type: choose between 'point' and 'line' for impulsive excitation, choose 'sin' for a signal excitation
sim_opts.diffuser_type = 'test';                        % type of diffuser geometry: choose between 'random', 'QRD', 'plane' and 'test'
sim_opts.incidence = 'normal';                          % incidence angle: choose between 'normal', 'oblique left' and 'oblique right'

% simulation parameters
sim_par.Lx = 20;                                        % room width (m)
sim_par.Ly = 20;                                        % room length (m)
sim_par.dur = 0.05;                                     % simulation duration (s) 
sim_par.xi = 0.5;                                       % coordinate of excitation (normalised, 0-1)
sim_par.yi = 1;
sim_par.c = 340;                                        % speed of sound (m/s)
sim_par.fmax = 4000;                                    % maximum frequency of accurate sampling
sim_par.Np = 8;                                         % number of sample points per spatial wavelength
sim_par.lambda = 1/sqrt(2);                             % Courant number (<=1/sqrt(2))
sim_par.tarfreq = 630;                                  % simulation target frequency (Hz) (center frequency of the 1/3 octave band of interest)

% diffuser Parameters
diff_par.desfreq = 500;                                 % diffuser design frequency (Hz)
diff_par.w0 = 0.05;                                     % well width (m)
diff_par.rR = 5;                                        % distance of receivers from the diffuser center (m)
diff_par.sR = 10;                                       % distance of sources from the diffuser center (m)
diff_par.dx = 0.5;                                      % normalised coordinate of diffuser center (dx - x coordinate)
diff_par.dy = 0.5;     
diff_par.reso = 5;                                      % angular resolution of receivers (degrees): >0,<=5
diff_par.nwell = 7;                                     % number of wells for a QRD diffuser (a prime number)
diff_par.period = 1;                                    % number of periods for a QRD diffuser
%%%% END EFIT %%%%

%%%% DO NOT EDIT %%%%
% diffuser option: simulate with/without the diffuser
diff_opts(1).diffuser_on = 1;                              
diff_opts(2) = diff_opts(1);
diff_opts(2).diffuser_on = 0;

% derived parameters
sim_par.SR = round(sim_par.Np*sim_par.fmax/sim_par.lambda);                     % sample rate (Hz)
sim_par.T = 1/sim_par.SR;                                                       % time step
sim_par.hmin = sim_par.c*sim_par.T/sim_par.lambda;                              % minimal space grid (for stability)
diff_par.maxwl = sim_par.c/diff_par.desfreq;                                    % maximum wavelength (corresponding to design frequency)
diff_par.maxwd = diff_par.maxwl/2/diff_par.nwell*max(qrs(diff_par.nwell,diff_par.period));    % maximum QRD well depth
diff_par.theta = deg2rad([-90:diff_par.reso:90]');                              % angular distribution of receivers
sim_par.ltheta = length(diff_par.theta);                                        % number of receivers

% ensure the diffuser size dose not exceed the measurement polar diameter
assert(diff_par.w0*diff_par.nwell*diff_par.period<2*diff_par.rR,'Diffuser size exceeds maximum') 
%%%% END DO NOT EDIT %%%%

%% define the well depth sequence according to diffuser type (all types have the same length and maximum well depth)
if strcmp(sim_opts.diffuser_type,'random')   % DO NOT EDIT % generate random well depth sequence
    welldepth = diff_par.maxwd*rand(diff_par.nwell*diff_par.period,1);
elseif strcmp(sim_opts.diffuser_type,'QRD')  % DO NOT EDIT % generate QRD well depth sequence
    welldepth = diff_par.maxwl/2*qrs(diff_par.nwell,diff_par.period)/diff_par.nwell;
elseif strcmp(sim_opts.diffuser_type,'plane')   % DO NOT EDIT % generate plane surface
    welldepth = zeros(diff_par.nwell*diff_par.period,1);
elseif strcmp(sim_opts.diffuser_type,'test')  % EDIT % user-defined well depth sequence
      welldepth = [0.175995 0.127410 0.006693 0.113718 0.048786 0.110522 0.160456];
%     welldepth = [0.1760 0.1274 0.0067 0.1137 0.0488 0.1105 0.1605];                                                                                          % N = 7 wideband-optimized diffuser
%     welldepth = [0.1173 0.2291 0.0331 0.0677 0.1477 0.0233 0.1200 0.0668 0.1600 0.1740 0.2437];                                                              % N = 11 wideband-optimized diffuser
%     welldepth = [0.1323 0.2584 0.0373 0.0764 0.1666 0.0263 0.1354 0.0753 0.1805 0.1963 0.2749 0.307 0.2915];                                                 % N = 13 wideband-optimized diffuser
%     welldepth = [0.2923 0.2717 0.1221 0.2404 0.1515 0.0038 0.3188 0.2782 0.1566 0.0138 0.1395 0.1507 0.0839 0.3082 0.2308 0.2578 0.2616];                    % N = 17 wideband-optimized diffuser
%     welldepth = [0.0842 0.0422 0.0947 0.1763 0.2498 0.1674 0.2437 0.1017 0.2960 0.2868 0.1877 0.0989 0.1407 0.0363 0.1113 0.1291 0.0340 0.2220 0.2595];      % N = 19 wideband-optimized diffuser
% or type your own depth sequence
end

%%%% DO NOT EDIT %%%%
%% polar response measurement
h1 = FDTD_func(sim_opts,sim_par,diff_par,diff_opts(1),welldepth);             % test 1: with diffuser
h2 = FDTD_func(sim_opts,sim_par,diff_par,diff_opts(2),welldepth);             % test 2: without diffuser
IR = h1-h2;                                                                          % polar response matrix of the diffuser(1 column for each receiver)

% calculate sound pressure levels at reach receiver
if strcmp(sys_opts.freq_range,'1/3 octave')
    octFilt = octaveFilter('CenterFrequency',sim_par.tarfreq,'Bandwidth','1/3 octave','SampleRate',sim_par.SR); 
    FilIR = octFilt(IR);                                                              % filter the impulse responses to the 1/3 octave band of interest
    I = sum(FilIR.^2);                                                                % sound pressure level                                                       
elseif strcmp(sys_opts.freq_range,'wideband') 
    I = sum(IR.^2);                                                                              
end 

%% calculate directional diffusion coefficient
dc = (sum(I)^2-sum(I.^2))/((sim_par.ltheta-1)*sum(I.^2))                                      % calculate the diffusion coefficient
assert(0<=dc<=1);

% produce polar plot 
xx = deg2rad([-90:90]);
yy = interp1(diff_par.theta,20*log10(I),xx,'spline');                                 % use spline interpolation to connect the discrete values
figure;
ax = polaraxes;
polarplot(diff_par.theta,20*log10(I),'o',xx,yy);                                      % plot the sound pressure level in dB
ax.ThetaDir = 'clockwise'; ax.ThetaZeroLocation = 'top'; ax.ThetaLim = [-90 90],ax.RLim = [min(20*log10(I))-10,max(20*log10(I))+10],ax.RAxis.Label.String = 'Sound pressure level(dB)',thetaticks(ax,[-90:5:90]);

fprintf('target frequency = %d\n',sim_par.tarfreq);
fprintf('diffuser type = %s\n',sim_opts.diffuser_type);
fprintf('diffusion coefficient = %f\n',dc);
%%%% END DO NOT EDIT %%%%
