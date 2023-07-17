%
% Program info
% Program author: Yuqing Li, August 2019
% Program details: This matlab script performs diffuser geometry optimization by searching for minimum
% (1-dc) value returned by function Yuqing_opti_dc_calculation_func.m for the test
% frequency, outputs maximum dc and the corresponding well depth sequence.

%%%% NOTICE:
% The optimization algorithm in this script only allows one input variable,
% 'welldepth', to be sent to the objective function Yuqing_opti_dc_calculation.m. 
% Therefore, the type of frequency range for optimization 'sys_opts.freq_range'
% and the target 1/3 octave center frequency 'sim_par.tarfreq' need to be set in Yuqing_opti_dc_calculation_func.m.

% If 'sys_opts.freq_range' is set to 'wideband', there is no need to edit 'sim_par.tarfreq'.

% Also, if the user would like to change the diffuser design frequency, remember to change both 'desfreq' in
% this script and 'diff_par.desfreq'  in Yuqing_opti_dc_calculation_func.m.

close all
clear all

%% define parameters
%%%% DO NOT EDIT %%%%%
% fixed parameter
c = 340;                                        % speed of sound (m/s)
%%%% END DO NOT EDIT %%%%

%%%% EDIT %%%%
% diffuser parameters
nwell = 7;                                      % number of wells
desfreq = 500;                                  % diffuser design frequency (Hz)

% diffuser type
dtype = 'QRD';                                  % type of initial diffuser geometry: choose between 'random' and 'QRD'
% optimization parameter
k = 20;                                         % number of start points
%%%% END EDIT %%%%

%%%% DO NOT EDIT %%%%
%% call functions and search for minimum value of (1-dc)
% initialize well depths according to initial diffuser type
if strcmp(dtype,'random')
    welldepth = c/desfreq/2*rand(nwell*period,1);
elseif strcmp(dtype,'QRD')
    welldepth = c/desfreq/2*Yuqing_qrs(nwell,1)/nwell;
end 
maxwd = max(welldepth);
dc0 = 1-Yuqing_opti_dc_calculation_func(welldepth);             % initial diffusion coefficient


% use multistart to find global minimum
rng default              %for reproducibility
opts = optimoptions(@fmincon,'Algorithm','sqp');            % use fincon to search for local minima under bounded conditions (maximum well depth of the optimized diffuser
                                                            % is constrained by that of the corresponding QRD)
tic
problem = createOptimProblem('fmincon','objective',...
    @Yuqing_opti_dc_calculation_func,'x0',welldepth,'lb',zeros(length(welldepth),1),'ub',maxwd*ones(length(welldepth),1),'options',opts);
ms = MultiStart('Display','iter','FunctionTolerance',1e-4,'StartPointsToRun','bounds-ineqs','UseParallel',true);      

% optimization result
[owd,minoneminusdc,eflag,output,manymins] = run(ms,problem,k) % owd: optimal well depth sequence found; minoneminusdc: minimal (1-dc)

%% output
totalruntime = toc;
maxdc = 1-minoneminusdc;

fprintf('number of wells = %d\n', nwell);
fprintf('total time for finding global minimum = %f\n', totalruntime);

%%%% END DO NOT EDIT %%%%