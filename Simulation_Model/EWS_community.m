%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main file for Solving SES System. 
% For Evaluating the SES system.
% This is abstracted from original reference
% Material: .
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: R. Kett (2017)
% Original Authors: A. Richter, V. Dakos (2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clears current memory.
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Global variables must be declared for 
% each scope they occupy 
% See: SES_dynamics_vanilla as an example.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global  a b n p q r w1 ema mvp K alpha beta gamma noiseRes noisePrice RnoiseProfc RnoiseProfd RnoisePrice sd_resource eopt ...
        nSON   pSON   w1SON   gammaSON   betaSON   emaSON   noisePriceSON ...
        nOMNRF pOMNRF w1OMNRF gammaOMNRF betaOMNRF emaOMNRF noisePriceOMNRF ...
        harvestOMNRFTEST harvestOMNRFTESTCoop harvestOMNRFTESTDef T v eOMNRF ecOMNRF ...
        initYear

tic
rng('shuffle')

% Model switch sets which model is being simulated: 0 for vanilla, 1 for redux.
modelSwitch = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vanilla Model Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Their correpsonding meaning is commented 
% beside each variable.
% At q = 0.0093, cannot use store anymore!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n     = 100;   % number of agents
p     = 5000;  % resource price
q     = 0.01;  % catchability coefficient
w1    = 1;     % wage cost
r     = 0.2;   % intrinsic growth rate
gamma = 0.1;   % strength of moral pressure
beta  = 0.2;   % strength of temptation
K     = 10;    % carrying capacity TODO: Change to carryingCapacity
mvp   = 0.1;   % minimum viable poulation
ema   = 0.6;   % effort endowment
a     = -0.3;  % parameter in harvest control rule
alpha = 1;     % stock effect
msy   = - r*(sqrt(K^2 - K*mvp + mvp^2) + K + mvp)*(sqrt(K^2 - K*mvp + mvp^2) + K - 2*mvp)*(sqrt(K^2 - K*mvp + mvp^2) - 2*K + mvp)/(27*K); % maximum sustainable yield

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ross Model Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Their correpsonding meaning is commented 
% beside each variable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSON       = 100;  % Number of SON agents.
pSON       = 5000; % SON Resource Price.
w1SON      = 1;    % SON Wage Cost.
gammaSON   = 0.5;  % SON Moral Pressure Strength. not van atm.
betaSON    = 0.1;  % SON Temptation Strength.
emaSON     = 0.6;  % SON Effort Endowment.
nOMNRF     = 100;  % Number of OMNRF agents.
pOMNRF     = 5000; % OMNRF Resource Price.
w1OMNRF    = 1;    % OMNRF Wage Cost.
gammaOMNRF = 0.4;  % OMNRF Moral Pressure Strength. %0.5 and 0.2 for beta give stable with 0.3 and 0.05 for OMNRF
betaOMNRF  = 0.1;  % OMNRF Temptation Strength.
emaOMNRF   = 0.6;  % OMNRF Effort Endowment.
delta      = 0.5;  % OMNRF harvest range modifier.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stochastic Model Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Their correpsonding meaning is commented 
% beside each variable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sd_resource = 0.075; % standard deviation resource
sd_price    = 20;    % standard deviation price

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Their correpsonding meaning is commented 
% beside each variable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MC      = 10;                    %number of Monte Carlo simulations TODO: Change to monteCarlo
options = odeset('RelTol',1e-6); % input for ODE solvee
Tspan   = 1:100;                 % length of simulation TODO: Change to timeSpan.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Storage Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Their correpsonding meaning is commented 
% beside each variable.
% Used for figure creation and data analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
permStore = [];     % Initialize empty array for permStore.
tempStore = [];     % Initialize empty array for tempStore.
storeOc   = [];     % Initialize empty array for StoreOc. 
figure1   = figure; % Initialize figure 1.
figure2   = figure; % Initialize figure 2.
figure3   = figure; % Initialize figure 3.

% Creating File Directories.
dirCount = 0;
dirExist = 0;
while dirExist ~= 1 
    if exist(strcat('/Users/woof-redux/Documents/SESY_V2/data/',datestr(date),'_','VANILLA_run_', num2str(dirCount))) == 0 && modelSwitch == 0 
        dirExist = 1;
        mkdir(strcat('/Users/woof-redux/Documents/SESY_V2/data/',datestr(date),'_','VANILLA_run_', num2str(dirCount)))
    elseif exist(strcat('/Users/woof-redux/Documents/SESY_V2/data/',datestr(date),'_','VANILLA_run_', num2str(dirCount)))  == 7 && modelSwitch == 0
        dirExist = 0;
        dirCount = dirCount+1;
    elseif exist(strcat('/Users/woof-redux/Documents/SESY_V2/data/',datestr(date),'_','ROSS_run_', num2str(dirCount)))  == 0 && modelSwitch == 1
        dirExist = 1;
        mkdir(strcat('/Users/woof-redux/Documents/SESY_V2/data/',datestr(date),'_','ROSS_run_', num2str(dirCount)))
    elseif exist(strcat('/Users/woof-redux/Documents/SESY_V2/data/',datestr(date),'_','ROSS_run_', num2str(dirCount)))  == 7 && modelSwitch == 1
        dirExist = 0;
        dirCount = dirCount+1;
    end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Ross SDE System.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Initialize arrays for current 
%    monte carlo simulation.
% 2. Generate noise.
% 3. Solve ODE by yearly integrand bounds.
% 4. Store yearly biomass, coop and step data.
% 5. Plot all MC runs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : MC
    endOfYear      = 0; % Resets endOfYear to zero for next MC runs.
    storeCoop      = [];
    storeCoopSON   = [];
    storeCoopOMNRF = [];
    storeBiomass   = [];
    storeTime      = [];
    tempStore      = [];

    % Choice of initial solutions.
    if modelSwitch == 0 
        solutionVanilla = [K, n/2];
    elseif modelSwitch == 1
        solutionRoss = [0.01, nSON/2, nOMNRF/2];
    end

    % All noise is generated for the span of the current simulation. 
    % IE. if Tspan = 100, we have 100 noise points.
    RnoiseRes = randn(length(Tspan),1);

    % Each Managements Placeholder Price Noise.
    RnoisePrice      = sd_price*randn(length(Tspan),1); % Vanilla Noise Input Placeholder Price.
    RnoisePriceSON   = sd_price*randn(length(Tspan),1); % SON Noise Input Placeholder Price.
    RnoisePriceOMNRF = sd_price*randn(length(Tspan),1); % OMNRF Noise Input Placeholder Price 
    xopt = (sqrt(K^2 - K*mvp + mvp^2) + K + mvp)/3;     % maximum sustainable biomass, note that msy=r*(k - mvp)^2/(4*K).
    eopt = r*xopt*(xopt-mvp)*(1-xopt/K)/(q*xopt*n);     % optimal effort.
    b    = (eopt-a)/xopt;                               % routine to determine socially optimally extraction level.

    % Initial eOMNRF setting.
    eOMNRF = -1;

    % This loop must contain the total harvest and coop harvest for OMNRF.
    for j = 1 : length(Tspan)
        v = 0;        % Tricky global scheme for tracking yearly harvest.
        initYear = 1; % Track when we enter a new range -> [j, j+1]
        harvestOMNRFTESTDef = [];
        harvestOMNRFTESTCoop = [];
        harvestOMNRFTEST = [];

        % Initial solution -> becomes last entry of previous j step.
        if modelSwitch == 0       
            initstate = solutionVanilla(end,:);
        elseif modelSwitch == 1
            initstate = solutionRoss(end,:);       
        end

        noiseRes = RnoiseRes(j); % Current years resource noise.

        % Each Managements Price Noise.
        noisePrice      = RnoisePrice(j);
        noisePriceSON   = RnoisePriceSON(j); % Current Years SON Price Noise.
        noisePriceOMNRF = RnoisePriceOMNRF(j); % Current Years OMNRF Price Noise.

        if modelSwitch == 0 
            [T, solutionVanilla] = ode45('SES_dynamics_vanilla', [j,j+1], initstate, options); % Run ODE routine
            storeBiomass         = [storeBiomass; solutionVanilla(:,1)];
            storeCoop            = [storeCoop; solutionVanilla(:,2)];
            fprintf('j = in MS0\n', j)
        elseif modelSwitch == 1
            [T, solutionRoss] = ode45('SES_dynamics_ross', [j,j+1], initstate, options); % Run ODE routine
            storeBiomass      = [storeBiomass; solutionRoss(:,1)];
            storeCoopSON      = [storeCoopSON; solutionRoss(:,2)];
            storeCoopOMNRF    = [storeCoopOMNRF; solutionRoss(:,3)];
            fprintf(']');
        end

        if modelSwitch == 1 
            if harvestOMNRFTEST(1) < delta*harvestOMNRFTESTCoop(end)
                eOMNRF = 0;
            elseif (harvestOMNRFTEST(1) >= delta*harvestOMNRFTESTCoop(end)) && harvestOMNRFTEST(end) <= harvestOMNRFTESTCoop(end)
                eOMNRF = 1;
            elseif harvestOMNRFTEST(1) > harvestOMNRFTESTCoop(end)
                eOMNRF = 2;
            end
        end
       
        % TOTAL occurrences of all eOMNRF variables (ie. lets check how often it changes.)
        storeOc = [storeOc; eOMNRF];
        storeTime = [storeTime; T];
    end

    % Store Data in array per MC (k) and print to file.
    if modelSwitch == 0
        tempStore = [storeTime, storeCoop, storeBiomass];
        csvwrite(strcat('/Users/woof-redux/Documents/SESY_V2/data/',datestr(date),'_','VANILLA_run_', num2str(dirCount),'/','Vanilla_MC_',num2str(k),'.csv'),tempStore)
    elseif modelSwitch == 1
        tempStore = [storeTime, storeCoopSON, storeCoopOMNRF, storeBiomass];
        csvwrite(strcat('/Users/woof-redux/Documents/SESY_V2/data/',datestr(date),'_','ROSS_run_', num2str(dirCount),'/','Ross_MC_',num2str(k),'.csv'),tempStore)
    end
                
    % Plotting Scheme for Biomass and cooperator population(s).
    % Biomass may stay the same as in the above loop the if-else statement covers the change of initial state.
    figure(figure1);
    biomassPlot = plot(storeTime, storeBiomass);
    xlim([0, length(Tspan)]);
    title('x(t) versus t for 100 simulations');
    xlabel('Time (In Years)');
    ylabel('Population Biomass');
    hold on;

    if modelSwitch == 0
        figure(figure2); % Figure 2 for model type 0.
        coopPlot = plot(storeTime, storeCoop);
        xlim([0, length(Tspan)]);
        title('SON C(t) versus t for 100 simulations');
        xlabel('Time (In Years)');
        ylabel('Population of Cooperators');
        hold on;
    elseif modelSwitch == 1
        figure(figure2); % Figure 2 for model type 1.
        coopPlotOMNRF = plot(storeTime, storeCoopOMNRF);
        xlim([0, length(Tspan)]);
        title('OMNRF C(t) versus t for 100 simulations');
        xlabel('Time (In Years)');
        ylabel('Population of Cooperators');
        hold on;
        figure(figure3); % Figure 3 for model type 1.
        coopPlotOMNRF = plot(storeTime, storeCoopSON);
        xlim([0, length(Tspan)]);
        title('SON Cooperator population for 100 simulations of a time span of 100 years');
        xlabel('Time (In Years)');
        ylabel('Population of Cooperators');
        hold on;  
    end
end

fprintf('\n');
toc %End of time tracking, end of full simulation.



