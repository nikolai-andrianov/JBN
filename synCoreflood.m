function synCoreflood
% Set up and run a numerical simulation using MRST (https://www.sintef.no/projectweb/mrst)
% of water injecting in a rock sample, which is composed of two core plugs
% with similar porosity and permeability. The sample is initially saturated 
% with oil and water. Water is injected with a constant rate at the left
% face of the sample, and constant pressure is maintained at the right face
% of the sample. 
% 
% Numerically, the sample is represented with a 1D geometry with constant 
% porosity and permeability in all grid blocks. Both water and oil are
% considered to be incompressible fluids with constant viscosities. 
% The relative permebility functions are taken in the Brooks-Corey form;
% zero capillary pressure is used.
%
% The simulated values of differential pressure and average water saturation 
% in the sample are interpreted with the Johnson–Bossler–Naumann (JBN) method, 
% and the obtained relative permeabilities are compared with the ground truth 
% relative permeabilities, used to run the numerical simulation.

%{
Copyright 2024 Nikolai Andrianov, nia@geus.dk

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%}

% Path to MRST
addpath('C:\GEUS\Matlab\mrst-2023b'); startup

mrstModule add ad-core ad-props ad-blackoil

%% Fluid properties

% Densities and viscosities for water 
rhow = 1.0514*gram/(centi*meter)^3; % Brine density @ 200 bara, 60 degC (Rowe&Chou)
muw = 0.565*centi*poise;

% Oil parameters
rhoo = 755.04*kilo*gram/(meter)^3; 
muo = 1.249*centi*poise;


%% Sample data
Lx = (7.489 + 7.457) *centi*meter; 
D = 3.74*centi*meter;
A = pi*D^2/4;

phi = mean([35.4, 35.4])/100;                   % Average of porosities for core plugs
Kg = harmmean([1233.17, 1212.29])*milli*darcy;  % Average gas permeability
pv = phi * A * Lx;

%% Define the simulation geometry

N = 50; 
M = 1;      % Number of gridblocks in the transverse direction (y & z)

% Represent the core as a rectangular parallelepiped 
Ly = sqrt(A);
Lz = Ly;

% Refine the cells towards the edges of the core
dxi = Lx / (N - 5);
dx = [0.1, 0.2, 0.4, 0.8, ones(1, N - 8), 0.8, 0.4, 0.2, 0.1]' * dxi;
x = [0; cumsum(dx)];
y = [0: Ly/M: Ly]';
z = [0: Lz/M: Lz]';

% Add two ghost cells at the inflow and outlfow
dx = [0.1*dxi; dx; 0.1*dxi];
x = [0; cumsum(dx)] - 0.1*dxi;

G = computeGeometry( tensorGrid(x, y, z) );

%% Rock properties

% Create the rock structure
rock = makeRock(G, Kg, phi);

% Set rock compressibility
pref = 1 * barsa;
cr = 0;
rock.cr = cr;
rock.pref = pref;

% plotCellData(G, rock.perm/milli/darcy), view(3), colorbar

%% Fluid model

fluid = initSimpleADIFluid('phases','WO',    ... % Fluid phases: water and oil
                           'mu',  [muw, muo],      ... % Viscosities
                           'rho', [rhow,rhoo],     ... % Surface densities [kg/m^3]
                           'c',   [0, 0],          ... % Fluid compressibility[Cm, Cm],              
                           'pRef', 0,              ... % Reference pressure 
                           'cR',  0                ... % Rock compressibility
                           );


% Setup the Brooks-Corey relative permeabilities 
nw = 1;
no = 3;
Swr = 0;
Sor = 0.2;
krw = 0.9;
kro = 0.7;

f.krW  = myCoreyPhaseRelpermAD(nw, Swr, krw, Swr + Sor);
f.krO = myCoreyPhaseRelpermAD(no, Sor, kro, Swr + Sor);  
f.pcOW = @(p) 0;

figure(1)
sw = [Swr:0.01:1-Sor];
plot(sw, f.krW(sw))
hold on
plot(sw, f.krO(1-sw))
title('Relative permeabilities')
legend('True k_{rw}', 'True k_{ro}')
xlabel('S_w')

% % Plot the inverse of the total mobility
% lw = f.krW(sw) / muw;
% lo = f.krO(1-sw) / muo;
% lt = lw + lo;
% 
% figure
% plot(sw, 1 ./ lt)
% hold on
% xlabel('S_w')
% legend('\lambda_{t}')
% plot(Swr, muo/kro, 'ok')
% plot(1-Sor, muw/krw, 'ok')


%% Introduce a second saturation region for the boundary cells

% The 1st region corresponds to the interior cells, the second one - to the
% boundary cells with the same relative permeabilities as for the interior,
% but with pc = 0
fluid.krW = {f.krW, f.krW};
fluid.krO = {f.krO, f.krO};
fluid.pcOW = {f.pcOW, @(sw)0};

% Get the indices of the inflow and outflow boundaries
ib_inj = boundaryFaceIndices(G, 'West');
ib_prod = boundaryFaceIndices(G, 'East');

% Indices of cells adjacent to the inflow and outflow boundaries
ic_inj = max(G.faces.neighbors(ib_inj, :));  
ic_prod = max(G.faces.neighbors(ib_prod, :)); 

% Mark cells with the corresponding saturation region
reg = ones(G.cells.num, 1);
reg(ic_inj) = 2;
reg(ic_prod) = 2;

% Set up the different regions 
rock.regions = struct('saturation', reg);


%% Set up the two-phase oil-water model

model = TwoPhaseOilWaterModel(G, rock, fluid);

%% Initial conditions

% Set initial pressure to 200 bar
p0 = 200*barsa;
state.pressure = ones(G.cells.num,1) * p0;  

% Initially saturated with oil and residual brine 
state.s = repmat([Swr 1-Swr], G.cells.num,1);

state.s = repmat([0 1], G.cells.num,1);

%% Set up the injection schedule

dt_fp = 10 * hour;
qinj = 500; % ml/h
wat_inj = 1;

% The volumetric rates change at specified time instants
t0 = 0;
ti = t0 + [0; cumsum(dt_fp)];
qw = qinj .* wat_inj       * (milli*litre)/hour;
qn = qinj .* (1 - wat_inj) * (milli*litre)/hour;

qt = qw + qn;

% Volumetric composition of the injected fluid  
sat = [qw; qn] ./ [qt; qt];

% Time steps, corresponding to different flow periods
dt = diff(ti);

ndt = numel(dt);

% Dirichlet BC at the right, no-flow conditions otherwise
% The initial pressure is interpreted as the right BC pressure 
bc = pside([], G, 'East', p0, 'sat', [1 0]);  
%bc = pside([], G, 'East', p0, 'sat', []);  

% A template for the Neumann BC at the left
bc = fluxside(bc, G, 'West', 1* (centi*meter)^3/hour, 'sat', [1 0]);

% Create the schedule with the template BC
schedule = simpleSchedule(dt, 'bc', bc);

% Replicate the control structure to accomodate for changing BC
schedule.control = repmat(schedule.control, ndt, 1);

% Clear the steps in the schedule as they will be overwritten with the
% small time steps within each flow period 
schedule = rmfield(schedule, 'step');

% Max time step within a FP
maxdts = 10*minute;
%maxdts = 60*minute;

% Number and length of min time steps at a start of a FP
nmin = 30;
dtmin = 1; 

% Number of linearly increasing time steps in a ramp-up
Nr = 10;

for n = 1:ndt
    
    % Adjust the injected fluid rate and composition to the current flow period
    schedule.control(n).bc.value = [p0, qt(n)]';
    schedule.control(n).bc.sat = [[1 0]; sat(:, n)'];
    
    % Introduce smaller timesteps for each flow period
    %dts = rampupTimesteps(dt(n), maxdts);
    
    % Partitition the first time interval within a FP into nmin timesteps 
    % of length dtmin, followed by Nr linearly increasing time steps 
    dtmax = min(dt(n), maxdts);
    dtn = 2*(dtmax - dtmin*nmin) / Nr - dtmin;
    dd = (dtn - dtmin) / (Nr - 1);
    dtinc = ones(Nr, 1) * dd;
    dtt = (cumsum(dtinc) - dd) + dtmin;
    dts1 = [ones(nmin, 1)*dtmin; dtt];
    dts1(end) = dtmax - sum(dts1(1:end-1));
    
    % Even steps
    dt_rem = repmat(maxdts, floor((dt(n) - maxdts)/maxdts), 1);
    % Final ministep if present
    dt_final = dt(n) - dtmax - sum(dt_rem);
    if dt_final <= 0
        dt_final = [];
    end
    % Combined timesteps
    dts = [dts1; dt_rem; dt_final];    

    % Save the small time steps in the schedule for each flow period and
    % provide an index to the corresponding entry in schedule.control
    if ~isfield(schedule, 'step')
        schedule.step.val = dts;
        schedule.step.control = ones(numel(dts), 1) * n;
    else
        schedule.step.val = [schedule.step.val; dts];
        schedule.step.control = [schedule.step.control; ones(numel(dts), 1) * n];
    end
        
end

%% Run the simulation

% Set a max number of iterations greater than default (=6)
maxTimestepCuts = 10;

solver = NonLinearSolver('maxTimestepCuts', maxTimestepCuts);

[~, states] = simulateScheduleAD(state, model, schedule, ...
                                 'NonLinearSolver', solver);


%% Get the simulation results
ns = numel(states);
t = t0 + cumsum(schedule.step.val);

% Get the time-dependent parameters
[pw_inj, pw_prod, pn_inj, pn_prod, ...
 qw_inj, qn_inj, qw_prod, qn_prod, sw_aver, pvi] = deal(zeros(ns, 1));

% Initial water volume in core
sw0 = mean(state.s(:, 1));

t_prev = 0;
for n = 1:ns
    
    % Influxes and outfluxes
    qw_inj(n) = sum(states{n}.flux(ib_inj, 1));
    qn_inj(n) = sum(states{n}.flux(ib_inj, 2));
    qw_prod(n) = sum(states{n}.flux(ib_prod, 1));
    qn_prod(n) = sum(states{n}.flux(ib_prod, 2));
    
    % Pore volumes injected
    pvi(n) = (qw_inj(n) + qn_inj(n)) * (states{n}.time - t_prev) / pv;
    t_prev = states{n}.time;  
    
    % Average water saturation
    sw_aver(n) = mean(states{n}.s(:, 1));
        
    % Pressures
    pn_inj(n) = mean(states{n}.pressure(ic_inj));
    sw = mean(states{n}.s(ic_inj, 1));
    if isfield(fluid, 'pcOW')
        pc = fluid.pcOW{2}(sw);
    else
        pc = 0;
    end
    pw_inj(n) = pn_inj(n) - pc;

    pn_prod(n) = mean(states{n}.pressure(ic_prod));
    sw = mean(states{n}.s(ic_prod, 1));
    if isfield(fluid, 'pcOW')
        pc = fluid.pcOW{2}(sw);
    else
        pc = 0;
    end
    pw_prod(n) = pn_prod(n) - pc;
    
end

pvi = cumsum(pvi);

%% Plot the simulation results

figure(5)
plot(t/hour, (pw_inj - pw_prod)/barsa)
hold on
xlabel('Time from start (hours)')
ylabel('Differential pressure (bar)')
legend('Simulation')

figure(6)
plot(t/hour, sw_aver)
hold on
xlabel('Time from start (hours)')
ylabel('Average S_w')

figure(7)
plot(t/hour, qw_prod / ((milli*litre)/hour))
hold on
xlabel('Time from start (hours)')
ylabel('q_{w,prod} (ml/h)')


%% Save the simulated data

T = array2table([t/hour, pvi, sw_aver, (pw_inj - pw_prod)/barsa, ...
    qw_prod / ((milli*litre)/hour), qn_prod / ((milli*litre)/hour)], ...
    'VariableNames', {'Time (hour)', 'PV_inj', 'Sw_aver', 'dP (bar)', ...
    'qw_prod (ml/h)', 'qn_prod (ml/h)'});

% Leave only the rows for FP #1
ind_fp2 = find(schedule.step.control == 2, 1);
T([ind_fp2:end], :) = [];

out_fname = 'sim_results.csv';
writetable(T, out_fname, 'Delimiter', ',');
disp(['Simulation results are saved in ' out_fname])

%% Compare with JBN rel perms

% Total injection velocity
vt = qw / A;

% Single-phase permeability to oil with Swc
Ko = kro * Kg;  

% The pressure drop, which corresponds to the single-phase permeability to
% oil with Swr
dPo = vt * (muo / Ko) * Lx / barsa;

% Assign short names for the simulation data
T.Properties.VariableNames = {'Time', 'PV_inj', 'Sw_aver', 'dP', 'qw_prod', 'qn_prod'}; 

% Get the JBN rel perms
[Swo, krwJBN, kroJBN] = JBN(T.Sw_aver, T.PV_inj, T.dP, T.qw_prod, dPo, ...
    muw, muo, krw, kro);

% Plot the JBN rel perm with the same colors as the true rel perms
figure(1)
set(gca,'ColorOrderIndex', 1)
plot(Swo, krwJBN, 'o', 'DisplayName', 'JBN k_{rw}')
plot(Swo, kroJBN, '*', 'DisplayName', 'JBN k_{ro}')

end





