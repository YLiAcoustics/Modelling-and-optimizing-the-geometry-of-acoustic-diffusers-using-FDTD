%
% Program info
% Program author: Yuqing Li, August 2019
% Program details: This matlab function takes the input well depth sequence from opti.m and pass it 
% to function opti_FDTD_func.m to calculate the diffusion coefficient dc. 
% It returns 1-dc to opti.m for optimization.
function oneminusdc = opti_dc_calculation_func(welldepth)

%%%% EDIT %%%%
%% global options
sys_opts.freq_range = '1/3 octave';                     % frequency range of diffusion coefficient /optimization: choose between '1/3 octave' and 'wideband'

%% simulation options
sim_opts.use_gpuarray = 0;                              % use GPU calculation
sim_opts.plot_on = 0;                                   % plot FDTD simulation
sim_opts.calculate_efficiency = 1;                      % calculate simulation efficiency
sim_opts.spatial_operater = 'centered';                 % form of spatial difference operator for boundary update: choose between 'centered' and 'noncentered'
sim_opts.source_type = 'line';                          % type of the sound source: choose between 'point', 'sin' and 'line'
sim_opts.incidence = 'normal';                          % incidence angle: choose between 'normal', 'oblique left' and 'oblique right'

%% simulation parameters
sim_par.Lx = 20;                                        % room width (m)
sim_par.Ly = 20;                                        % room length (m)
sim_par.dur = 0.05;                                     % simulation duration (s) 
sim_par.c = 340;                                        % speed of sound (m/s)
sim_par.fmax = 4000;                                    % maximum frequency of accurate sampling
sim_par.Np = 8;                                         % number of sample points per spatial wavelength
sim_par.lambda = 1/sqrt(2);                             % Courant number (<=1/sqrt(2))
sim_par.tarfreq = 630;                                 % simulation target frequency (Hz) (center frequency of the 1/3 octave band for optimization)

%% diffuser parameters
diff_par.desfreq = 500;                                  % diffuser design frequency (Hz)
diff_par.w0 = 0.05;                                      % well width (m)
diff_par.rR = 5;                                         % distance of receivers from the diffuser center (m)
diff_par.sR = 10;                                        % distance of sources from the diffuser center (m)
diff_par.dx = 0.5;                                       % normalised coordinate of diffuser center (dx - x coordinate)
diff_par.dy = 0.5;    
diff_par.reso = 5;                                       % angular resolution of receivers (degrees): >0,<=5
%%%% END EDIT %%%%

%%%% DO NOT EDIT %%%%
% derived simulation parameters
sim_par.SR = round(sim_par.Np*sim_par.fmax/sim_par.lambda);           % sample rate (Hz)
sim_par.T = 1/sim_par.SR;                                             % time step
sim_par.hmin = sim_par.c*sim_par.T/sim_par.lambda;                    % minimal space grid (stability)
diff_par.theta = deg2rad([-90:diff_par.reso:90]');                    % angular distribution of receivers
sim_par.ltheta = length(diff_par.theta);                              % number of receivers

% diffuser option: simulate with/without the diffuser
diff_opts(1).diffuser_on = 1;   
diff_opts(2) = diff_opts(1);
diff_opts(2).diffuser_on = 0;

%% polar response measurement
h1 = opti_FDTD_func(sim_opts,sim_par,diff_par,diff_opts(1),welldepth);
h2 = opti_FDTD_func(sim_opts,sim_par,diff_par,diff_opts(2),welldepth);
IR = h1-h2;                                              % the receiver output matrix (Nf rows, nrec columns


if strcmp(sys_opts.freq_range,'1/3 octave')
    octFilt = octaveFilter('CenterFrequency',sim_par.tarfreq,'Bandwidth','1/3 octave','SampleRate',sim_par.SR); 
    FilIR = octFilt(IR);                                                              % filter the impulse responses to the 1/3 octave band of interest
    I = sum(FilIR.^2);                                                                % sound pressure level                                                       
elseif strcmp(sys_opts.freq_range,'wideband') 
    I = sum(IR.^2);                                                                              
end 

% calculate directional diffusion coefficient
dc = double((sum(I)^2-sum(I.^2))/((sim_par.ltheta-1)*sum(I.^2)));          % convert dc to double type because fmincon only takes double type argument
assert(0<=dc<=1);
oneminusdc = 1-dc;
%%%% END DO NOT EDIT %%%%
end