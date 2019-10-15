%
% Program info
% Program author: Yuqing Li, August 2019
% Program details: This matlab function simulates sound propagation and
%                  diffusion using FDTD room acoustics modelling in 2D. It
%                  returns a matrix of recorded signals from all receivers
%                  (1 column for each reciever).

function rec = Yuqing_opti_FDTD_func(sim_opts,sim_par,diff_par,diff_opts,welldepth)    

%% Define Parameters
% 1. copy over parameters
% simulation parameters
SR = sim_par.SR;                                       % sample rate
c = sim_par.c;                                         % speed of sound
lambda = sim_par.lambda;                               % Courant number
hmin = sim_par.hmin;                                   % minimal space grid (for stability)

dur = sim_par.dur;                                     % simulation duration    
Lx = sim_par.Lx;                                       % room dimensions (m)        
Ly = sim_par.Ly;                  

% diffuser parameters
dx = diff_par.dx;                                       % normalised coordinate of diffuser center (dx - x coordinate)
dy = diff_par.dy;                                                
w0 = diff_par.w0;                                       % well width (m)
desfreq = diff_par.desfreq;                             % diffuser design frequency

% measurement parameters
rR = diff_par.rR;                                       % receiver radius (m)
theta = diff_par.theta;                                 % angular resolution of receivers
sR = diff_par.sR;                                       % source radius(m)
nrec = sim_par.ltheta ;                                 % number of receivers

% 2. derived parameters
% simulation parameters
Nf = round(SR*dur);                                    % simulation duration in samples
N = floor(Lx/hmin)-1;                                  % number of grid points to update in x direction
M = floor(Ly/hmin)-1;                                  % number of grid points to update in y direction

% diffuser geometry parameters
nw = length(welldepth);                                 % number of wells
w = round(w0/hmin);                                     % well width in number of grid points
maxwd = c/desfreq/2/length(welldepth)*max(Yuqing_qrs(length(welldepth),1));           % 
MAXWD = round(maxwd/hmin);                              % maximum well depth in number of grid points
H = MAXWD+3;                                            % maximum bar height/diffuser area height (at least 3 grid points to make sure the inner/outer point determiner works)

%% 2. Boundary classification 
% Boundary conditions                                    
gamma = 1;                                              % wall admittance (fully absorptive at normal incidence))
gammad = 0;                                             % diffuser surface admittance (fully reflective)

% diffuser position
dcenter = [dx*(N+2),dy*(M+2)];                                % diffuser center coordinate (index 1: column number). 
dlength = nw*w;                                               % diffuser length in numbers of points

dstartx = round(dcenter(1)-dlength/2);                         % diffuser start point (dstartx: x coordinate (column number)
dendx = dstartx+nw*w-1;                                        % diffuser end point
dbottom = ceil(dcenter(2)-H/2);                                % corrected dbottom position
dtop = dbottom+H-1;

% node classification mask  
alpha = 4*ones(M+2,N+2,'uint8');                               % 4 = air/outer corners, 3 = walls, 2 = inner corners
d = ones(M+2,N+2);                                             % inner/outer point determiner. A layer of outer points are added to the outside of the domain.
D = ones(H,nw*w);                                              % inner/outer point determiner for the diffuser
    
if diff_opts.diffuser_on
    WD = round(welldepth/hmin);                                % well depths in number of grid points
    BH = H-WD;                                                 % bar heights in number of grids
    dg = round(repelem(BH,w));                                 % diffuser geometry taking into account well widths
    for a = 2:nw*w-1
        D(2:dg(a)-1,a) = 0;                                    % set inner points to 0
     end
    d(dbottom:dbottom+H-1,dstartx:dendx) = D;
end

% boundary inner points
d(1,:) = 0;
d(M+2,:) = 0;
d(:,1) = 0;
d(:,N+2) = 0;

% sum up inner/outer points to get the node classification mask  
alpha(2:M+1,2:N+1) = d([2:M+1]+1,2:N+1)+d([2:M+1]-1,2:N+1)+d([2:M+1],[2:N+1]+1)+d([2:M+1],[2:N+1]-1);

% set alpha at inner points to 0
zero_d = find(~d);                                             % find inner points (d = 0)
alpha(zero_d) = 0;                                        
zero_alpha = find(alpha==0);


%% 3. FDTD Simulation
% 3.1 state variables
u0 = zeros(M+2,N+2,'single');                                  % state at time index n+1
u1 = zeros(M+2,N+2,'single');                                  % state at time index n
u2 = zeros(M+2,N+2,'single');                                  % state at time index n-1

% 3.2 impulsive excitation
% source position determined by angle of incidence 
if strcmp(sim_opts.incidence,'normal')
   ai = 0;
elseif strcmp(sim_opts.incidence,'oblique left')
   ai = -55;    
else
   ai = 55;
end
ai = deg2rad(ai);
sp = floor([dcenter(1)+sR/hmin*sin(ai) dcenter(2)+sR/hmin*cos(ai)]);    % source position. index 1: column number
li = sp(1);                                                       % source position coordinate
mi = sp(2);
assert((1<li<N+2) && (1<mi<M+2), 'source position beyond simulation domain');

% impulse source type
if strcmp(sim_opts.source_type,'point')                           % a point impulse
    u1(mi,li) = 1e5;
elseif strcmp(sim_opts.source_type,'line')                        % an impulse array parallel to the diffuser bottom  
    assert(strcmp(sim_opts.incidence,'normal'),'linear impulse not applicable for oblique incidence') 
    u1(mi,3:N) = 100;                                             
end

% 3.3 output
rp = round([dcenter(1)+rR/hmin*sin(theta) dcenter(2)+rR/hmin*cos(theta)]);   % receiver position matrix  (index 1: column number) 
lo = rp(:,1);                                                                % receiver ouput coordinate (lo: column number)
mo = rp(:,2);

assert((min(lo)>2) && (max(lo)<N+1) && (min(mo)>2) && (max(mo)<M+1), 'at least one receive position beyond simulation domain')

RP = (M+2)*(lo-1)+mo;                                                        % indices of receiver positions in the domain

rec = zeros(Nf,nrec,'single');                                               % initialize the output matrix

% plot simulation domain
if sim_opts.plot_on
    view_alpha = alpha;
    view_alpha(view_alpha==4) = 0;
    surface(view_alpha);
    xlim([1 N+2]),ylim([1 M+2]);
    hold on;
    han = surface(u0);                                             % create sound pressure surface plot
    view(2);
    colormap(flipud(bone));
    hold on; colorbar; shading interp;
    plot(dcenter(1),dcenter(2),'xk');                              % plot diffuser center position
    
    % plot receiver positions    
    plot(lo,mo,'*r'); hold off;        
    axis square
end

% use GPU calculation
if sim_opts.use_gpuarray
    gpudev = gpuDevice;
    u0 = gpuArray(u0);
    u1 = gpuArray(u1);
    u2 = gpuArray(u2);
    alpha = gpuArray(alpha);
end


%3.4 main loop
tic

for n = 1:Nf

    % update equation (frequency-independent impedance)
    % centered time, non-centered space grid
    u0(2:M+1,2:N+1) = (lambda^2*(u1([2:M+1]+1,2:N+1)+u1([2:M+1]-1,2:N+1)+u1(2:M+1,[2:N+1]+1)+u1(2:M+1,[2:N+1]-1))+ (2-lambda^2*single(alpha(2:M+1,2:N+1))).*u1(2:M+1,2:N+1)+(single(4-alpha(2:M+1,2:N+1))*gamma*lambda/2-1).*u2(2:M+1,2:N+1))./(1+single(4-alpha(2:M+1,2:N+1))*lambda*gamma/2); 

    if diff_opts.diffuser_on
    % diffuser area update (different gamma)
    u0(dbottom:dtop,dstartx:dendx) = (lambda^2*(u1(dbottom:dtop,[dstartx:dendx]+1)+u1(dbottom:dtop,[dstartx:dendx]-1)+u1([dbottom:dtop]+1,dstartx:dendx)+u1([dbottom:dtop]-1,dstartx:dendx))+(2-lambda^2*single(alpha(dbottom:dtop,dstartx:dendx))).*u1(dbottom:dtop,dstartx:dendx)+(single(4-alpha(dbottom:dtop,dstartx:dendx))*gammad*lambda/2-1).*u2(dbottom:dtop,dstartx:dendx))./(1+single(4-alpha(dbottom:dtop,dstartx:dendx))*lambda*gammad/2); 
    end
   

    if strcmp(sim_opts.spatial_operater,'centered')
    % centered time & space operater
    % walls
    % right
        u0(3:M,N+1) = (lambda^2*(2*u1(3:M,N)+u1([3:M]+1,N)+u1([3:M]-1,N))+2*(1-2*lambda^2)*u1(3:M,N+1)+(lambda*gamma-1)*u2(3:M,N+1))/(lambda*gamma+1); 
    % left
        u0(3:M,2) = (lambda^2*(2*u1(3:M,3)+u1([3:M]+1,2)+u1([3:M]-1,2))+2*(1-2*lambda^2)*u1(3:M,2)+(lambda*gamma-1)*u2(3:M,2))/(lambda*gamma+1);
    % top
        u0(M+1,3:N) = (lambda^2*(2*u1(M,3:N)+u1(M,[3:N]+1)+u1(M+1,[3:N]-1))+2*(1-2*lambda^2)*u1(M+1,3:N)+(lambda*gamma-1)*u2(M+1,3:N))/(lambda*gamma+1);
    % bottom
        u0(2,3:N) = (lambda^2*(2*u1(3,3:N)+u1(2,[3:N]+1)+u1(2,[3:N]-1))+2*(1-2*lambda^2)*u1(2,3:N)+(lambda*gamma-1)*u2(2,3:N))/(lambda*gamma+1);     

    % corners
    % top right
        u0(M+1,N+1) = (2*lambda^2*(u1(M+1,N)+u1(M,N+1))+2*(1-2*lambda^2)*u1(M+1,N+1)+(2*lambda*gamma-1)*u2(M+1,N+1))/(2*lambda*gamma+1);
    % top left
        u0(M+1,2) = (2*lambda^2*(u1(M+1,3)+u1(M,2))+2*(1-2*lambda^2)*u1(M+1,2)+(2*lambda*gamma-1)*u2(M+1,2))/(2*lambda*gamma+1);
    % bottom right 
        u0(2,N+1) = (2*lambda^2*(u1(2,N)+u1(3,N+1))+2*(1-2*lambda^2)*u1(2,N+1)+(2*lambda*gamma-1)*u2(2,N+1))/(2*lambda*gamma+1);
    % bottom left
        u0(2,2) = (2*lambda^2*(u1(3,2)+u1(2,3))+2*(1-2*lambda^2)*u1(2,2)+(2*lambda*gamma-1)*u2(2,2))/(2*lambda*gamma+1);
    end    

    % zero pressure at inner points
     u0(zero_alpha) = 0;

    % frequency-dependent impedance

    % read output 
    rec(n,:) = gather(u0(RP));   

    % shift states to step forward in time
    u2 = u1; 
    u1 = u0;
end

if sim_opts.use_gpuarray
    wait(gpudev);
end

runtime = toc; 
fprintf('runtime = %f\n',runtime);
 
%% 4. calculate computational speed (points/sec)
if sim_opts.calculate_efficiency
    fprintf('points per second = %d\n', N*M*Nf/runtime);
end


end