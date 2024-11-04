%% Clear workspace and command window
clear; clc;
load('LH.mat')
load('LH2.mat')
load('LNG.mat')
load('cbuy.mat')
load('csell.mat')
load('cbuy2.mat')
load('csell2.mat')
load('RNW.mat')
load('CO2_factor.mat')
load('LE.mat')

%% Select the day to model
SOCmax = 2000; % kWh
SOCmin = 0; % kWh
rate = SOCmax/2; % kW
period = 1; % Periods of 1 hour
T = 24; % 24 hours
day2model = 3; % 1 = a day of June, 2 = a day of October, 3 = a day of December
Tday = 24/period;

%% Optimisation Setup
UB = rate * ones(T, 1);
LB = -rate * ones(T, 1);
A = []; B = []; Aeq = []; Beq = [];

x = cell(1, 1);
fval = cell(1, 1);
options = optimoptions('gamultiobj', 'MaxGenerations', 100000, 'OutputFcns', @outputfunction);

for i = 1:10
    [x{i, :}, fval{i, :}] = gamultiobj(@(x) objectives(x, T, period, cbuy2(:, day2model), csell2(:, day2model), LE(:, day2model), LH(:, day2model), RNW(:, day2model), CO2_factor(:, day2model), SOCmax, SOCmin), T, [], [], [], [], LB, UB, @(x) const(x, T, period, SOCmax, SOCmin), options);
end

EP_f = fval;
alpha_c2 = cell2mat(x);
EP_f2 = cell2mat(fval); % All Pareto solutions concatenated

%% "Net Flow Method" to find the best-compromised solution

% Normalize the data using min-max normalization
numAlternatives = size(EP_f2, 1);
normalizedData = (EP_f2 - min(EP_f2)) ./ (max(EP_f2) - min(EP_f2));

% Define weights for each objective
weights = [0.5, 0.5]; % Equal weights for 2 objectives

% Calculate Positive Flow and Negative Flow
positiveFlow = zeros(numAlternatives, 1);
negativeFlow = zeros(numAlternatives, 1);

for i = 1:numAlternatives
    for j = 1:numAlternatives
        if i ~= j
            % Positive Flow: how well alternative i performs compared to j
            positiveFlow(i) = positiveFlow(i) + sum(weights .* (normalizedData(i, :) > normalizedData(j, :)));
            
            % Negative Flow: how much worse alternative i performs compared to j
            negativeFlow(i) = negativeFlow(i) + sum(weights .* (normalizedData(i, :) < normalizedData(j, :)));
        end
    end
end

% Calculate the Net Flow for each alternative
netFlow = positiveFlow - negativeFlow;

% Rank the alternatives based on Net Flow (higher is better)
[~, BestAlternativeIndex] = max(netFlow);

% Best alternative values (corresponding to the best net flow score)
Best = alpha_c2(BestAlternativeIndex, :);

% Rank all the alternatives
[~, rankedAlternatives] = sort(netFlow, 'descend');

clear numAlternatives numObj normalizedData weights positiveFlow negativeFlow netFlow

%% Identify the Best Grid (Cost) and Best CO2 Solutions

% Best Cost solution (minimum grid consumption)
[~, minimo1_idx] = min(EP_f2(:,1), [], 'linear');  % Index for the best grid solution
x_Bills = alpha_c2(minimo1_idx, :);  % Best decision variables for cost

% Best CO2 solution (minimum CO2 emissions)
[~, minimo2_idx] = min(EP_f2(:,2), [], 'linear');  % Index for the best CO2 solution
x_CO2 = alpha_c2(minimo2_idx, :);  % Best decision variables for CO2

%% Analyze and Compare the Selected Solutions
choices = zeros(3, T);

% Store the best solutions: Cost, CO2, and Compromised (from NFM)
choices(1, :) = x_Bills;   % Best Cost solution
choices(2, :) = x_CO2;     % Best CO2 solution
choices(3, :) = Best;      % Best Compromised solution

% Array to hold the total cost and CO2 for each solution
totalResults = zeros(3, 2);

for i = 1:3
    % Analyze the cost and CO2 emissions for each solution
    [totalCost, tCO2Gen, ~, ~, ~, ~, ~] = analysis(choices(i, :), T, period, cbuy2, csell2, LE, LH, RNW, CO2_factor);
    
    totalResults(i, 1) = totalCost;  % Store total cost
    totalResults(i, 2) = tCO2Gen;    % Store total CO2 emissions
end

% Output total results for all 3 solutions (Best Cost, Best CO2, Best Compromise)
disp('Results for Best Cost, Best CO2, and Best Compromised Solutions:');
disp(totalResults);

%% Calculate Key Variables for Plotting

opt = 3; % 1 = Bills, 2 = CO2, 3 = Best
choice = choices(opt, :);
Pbat = choice'; % Battery power usage over time

Pgrid_Electricity = zeros(T, 1); % Power drawn from the grid
surplus = zeros(T, 1); % Surplus energy (after loads are satisfied)
purchaseG_Electricity = zeros(T, 1); % Electricity purchased from the grid
costElectricity = zeros(T, 1); % Cost of electricity for each time period
revenueElectricity = zeros(T, 1); % Revenue from selling surplus electricity

Pgrid_Hydrogen = zeros(T, 1); % Electricity used for hydrogen production
H2_heating = zeros(T, 1); % Hydrogen used for heating in kWh
NG_heating = zeros(T, 1); % Natural gas used for heating in kWh
costHeating = zeros(T, 1); % Cost of heating for each time period
CO2Heating = 0;
CO2grid = 0;
CO2total = zeros(T, 1); % CO2 emissions from heating

% Battery storage state (SOC)
SOC = zeros(T+1, 1);
SOC(1, 1) = (SOCmax + SOCmin) / 2; % Initial battery state of charge (SOC)

% Define constants
electrolyzer_efficiency = 0.8; % Efficiency of hydrogen production
gas_boiler_efficiency = 0.9; % Efficiency of natural gas boiler
co2_ng_per_kWh = 185; % CO2 emissions from natural gas in gCO2/kWh

for i = 1:T
    % Calculate the electricity demand and deficit/surplus after PV
    Pdeficit = LE(i, day2model) - RNW(i, day2model); % Load minus renewable (PV)
    
    if Pdeficit > 0
        % Not enough PV, use battery first if available
        if SOC(i) > SOCmin % Use battery if SOC is above minimum
            Pbat(i) = min(Pbat(i), Pdeficit); % Use battery to cover part of the deficit
            SOC(i+1) = SOC(i) - Pbat(i) * period; % Update battery state of charge
        else
            Pbat(i) = 0; % No battery usage if SOC is too low
        end
        
        % Any remaining deficit is covered by the grid
        Pgrid_Electricity(i, 1) = Pdeficit - Pbat(i);
        purchaseG_Electricity(i, 1) = Pgrid_Electricity(i, 1) * cbuy2(i, day2model) * period;
        CO2grid = purchaseG_Electricity(i, 1) * CO2_factor(i, day2model);
    else
        % Surplus PV generation: charge battery first if capacity is available
        surplus(i, 1) = -Pdeficit;
        
        if SOC(i) < SOCmax
            Pcharge = min(surplus(i, 1), SOCmax - SOC(i)); % Charge battery with surplus
            SOC(i+1) = SOC(i) + Pcharge * period; % Update SOC
            surplus(i, 1) = surplus(i, 1) - Pcharge; % Reduce surplus by the amount charged
        else
            SOC(i+1) = SOC(i); % Battery is full, SOC remains constant
        end
        
        % If surplus remains after charging, sell it
        if surplus(i, 1) > 0
            revenueElectricity(i, 1) = surplus(i, 1) * csell2(i, day2model) * period;
        else
            revenueElectricity(i, 1) = 0;
        end
    end
    
    % Heating demand (if hydrogen is insufficient)
    remaining_heating_demand = max(0, LH(i, day2model) - (H2_heating(i, 1) / electrolyzer_efficiency));
    
    if remaining_heating_demand > 0
        % Use natural gas for heating
        NG_heating(i, 1) = remaining_heating_demand / gas_boiler_efficiency;
        costHeating(i, 1) = NG_heating(i, 1) * cbuy2(i, day2model) * period;
        CO2Heating = NG_heating(i, 1) * co2_ng_per_kWh * period;
    else
        costHeating(i, 1) = 0;
        NG_heating(i, 1) = 0;
    end
    
    CO2total(i,1) = (CO2Heating + CO2grid)/1000; % In KG
    
    % Store electricity cost
    costElectricity(i, 1) = purchaseG_Electricity(i, 1) - revenueElectricity(i, 1);
end

% Total cost and CO2 emissions
totalElectricityCost = sum(costElectricity);
totalHeatingCost = sum(costHeating);
totalRevenue = sum(revenueElectricity);
for i=1:T
    tGridCost(i,1) = purchaseG_Electricity(i,1) + costHeating(i,1) - revenueElectricity(i,1);
end
totalCO2 = sum(CO2Heating); % Total CO2 emissions from heating

% Total cost and CO2 emissions for the system
totalCost = totalElectricityCost + totalHeatingCost - totalRevenue;
tCO2 = sum(CO2total); % Combined CO2 from electricity and heating

%% Non-optimized scenario
% Using only grid
non_optim_solution = zeros(24, 2);  

for t = 1:T
    % Calculate grid usage for electricity demand minus renewable energy
    GE(t, 1) = LE(t, day2model) - RNW(t, day2model);  
    if GE(t,1) > 0
        PG_E(t,1) = GE(t,1); 
    else
        PG_E(t,1) = 0;
    end
    
    % Heating demand fully supplied by grid (no hydrogen production)
    PG_H(t,1) = LH(t,day2model) / gas_boiler_efficiency; 
    
    % Total grid cost for electricity and heating
    Pgrid_Tot(t, 1) = (PG_E(t, 1) * cbuy(t, day2model) * period) + ...
                      (PG_H(t, 1) * cbuy(t, day2model) * period);
    
    % CO2 emissions from grid usage (both electricity and heating)
    totGrid(t,1) = LE(t,day2model)+LH(t,day2model);
    if totGrid(t,1) > 0
        CO2_E(t,1) = (CO2_factor(t,day2model) * totGrid(t,1)); 
    else
        CO2_E(t,1) = 0;
    end
    CO2_H(t,1) =PG_H(t, 1)*co2_ng_per_kWh;
    CO2_Grid(t, 1) = (CO2_E(t,1)+CO2_H(t,1))* (period / 1000);

    % Store results in the non-optimized solution matrix
    non_optim_solution(t, 1) = Pgrid_Tot(t, 1);  % Total grid cost
    non_optim_solution(t, 2) = CO2_Grid(t, 1);  % CO2 emissions
end

%% Optimised and non optimised results
% Plot non-optimized result
figure
dts=1;
index = 0:period:(24*dts)-period;
if day2model==1 
    dayofyear = "June"; 
elseif day2model==2 
    dayofyear = "October";
elseif day2model==3 
    dayofyear = "December";
end
index = 0:period:(24*dts)-period;
yyaxis left
stairs(index, non_optim_solution(:, 1),'-^b','LineWidth',1);
ylabel('Cost (£/kWh)','fontsize',15,'Color','b');
hold on;
yyaxis right
stairs(index, non_optim_solution(:, 2),'-xr', 'LineWidth',1);
ylabel('CO2 emissions (kg)','fontsize',15,'Color','r');
title({"Non-Optimised results of a day in "+dayofyear+":", "Grid usage: "+sum(non_optim_solution(:, 1))+"£      CO2 emissions: "+sum(non_optim_solution(:, 2))+"KG"},'fontweight','bold','fontsize',12);
xlabel('Time (Hours)','fontsize',15,'Color','k');
grid on;
legend({'Grid cost', 'CO2 emissions'}, 'Location', 'northeast', 'FontSize', 12);

% Plot optimised results
figure
hold on;
stairs(index,tGridCost, '-vb', 'LineWidth',1);
title({"Optimised results of a day in "+dayofyear+":","Grid usage: "+totalCost+"£    CO2 Emissions: "+tCO2+"KG"},'fontweight','bold','fontsize',12);
xlabel('Time (Hours)','fontsize',15,'Color','k');
yyaxis left
ylabel('Cost (£/kWh)','fontsize',15,'Color','b');  
yyaxis right
ylabel('CO2 emissions (kg)','fontsize',15,'Color','r')
hold on;
stairs(index,CO2total, '-or', 'LineWidth',1);
legend({'Grid cost', 'CO2 emissions'}, 'Location', 'northeast', 'FontSize', 12);

%% Comparison plot
dts=1;
index = 0:period:(24*dts)-period;
reCO2 = (CO2total./1000);
if day2model==1 
    dayofyear = "June"; 
elseif day2model==2 
    dayofyear = "October";
elseif day2model==3 
    dayofyear = "December";
end
totalLoad = LE + LH;
figure
title("Results of a day in "+dayofyear,'fontweight','bold','fontsize',15);
yyaxis left
P1 = stairs(index, non_optim_solution(:, 1),'--^b','LineWidth',2);
ax = gca;
ax.FontSize = 15;
ax.FontWeight = 'bold';
hold on
ylabel('Cost (£/kWh) / Power (kW)','fontweight','bold','fontsize',15,'Color','b');  
G2 = stairs(index,totalLoad(:,day2model), '-squareg', 'LineWidth',2);
G3 = stairs(index,RNW(:,day2model),'--pentagramg', 'LineWidth',2);
P2 = stairs(index,tGridCost(:,1), '-vb', 'LineWidth',2);
yyaxis right
P3 = plot(index, non_optim_solution(:, 2),'--xr', 'LineWidth',2); 
ylabel('CO2 emissions (kg)','fontweight','bold','fontsize',15,'Color','r')
P4 = plot(index,CO2total, '-or', 'LineWidth',2);
hold off
legend([G2, G3, P1, P2, P3, P4], ...
    'Total power load (kW)', ...
    'Solar power (kW)', ...
    'Non-optimised grid cost (£)', ...
    'Optimised grid cost (£)', ...
    'Non-optimised CO2 emissions (Kg)', ...
    'Optimised CO2 emissions (Kg)')
xlabel('Time (Hours)','fontweight','bold','fontsize',15,'Color','k');

%% Objectives
function [objectives] = objectives(x, T, period, cbuy, csell, LE, LH, Ppv, CEI, SOCmax, SOCmin)
    % Initialize variables
    Pbat = x(1:T)'; % Battery power usage over time
    Pgrid_Electricity = zeros(T, 1); % Power drawn from the grid
    surplus = zeros(T, 1); % Surplus energy (after loads are satisfied)
    purchaseG_Electricity = zeros(T, 1); % Electricity purchased from the grid
    revenueElectricity = zeros(T, 1); % Revenue from selling surplus electricity
    
    H2_heating = zeros(T, 1); % Hydrogen used for heating in kWh
    NG_heating = zeros(T, 1); % Natural gas used for heating in kWh
    costHeating = zeros(T, 1); % Cost of heating for each time period
    CO2total = zeros(T, 1);
    CO2grid = 0;
    
    % Battery storage state (SOC)
    SOC = zeros(T+1, 1);
    SOC(1, 1) = (SOCmax + SOCmin) / 2; % Initial battery state of charge (SOC)
    
    % Define constants
    electrolyzer_efficiency = 0.8; % Efficiency of hydrogen production
    gas_boiler_efficiency = 0.9; % Efficiency of natural gas boiler
    co2_ng_per_kWh = 185; % CO2 emissions from natural gas in gCO2/kWh
    
    for i = 1:T
        % Calculate the electricity demand and deficit/surplus after PV
        Pdeficit = LE(i, 1) - Ppv(i, 1); % Load minus renewable (PV)
        
        if Pdeficit > 0
            % Not enough PV, use battery first if available
            if SOC(i) > SOCmin % Use battery if SOC is above minimum
                Pbat(i) = min(Pbat(i), Pdeficit); % Use battery to cover part of the deficit
                SOC(i+1) = SOC(i) - Pbat(i) * period; % Update battery state of charge
            else
                Pbat(i) = 0; % No battery usage if SOC is too low
            end
            
            % Any remaining deficit is covered by the grid
            Pgrid_Electricity(i, 1) = Pdeficit - Pbat(i);
            purchaseG_Electricity(i, 1) = Pgrid_Electricity(i, 1) * cbuy(i, 1) * period;
            CO2grid = purchaseG_Electricity(i, 1) * CEI(i, 1);
        else
            % Surplus PV generation: charge battery first if capacity is available
            surplus(i, 1) = -Pdeficit;
            
            if SOC(i) < SOCmax
                Pcharge = min(surplus(i, 1), SOCmax - SOC(i)); % Charge battery with surplus
                SOC(i+1) = SOC(i) + Pcharge * period; % Update SOC
                surplus(i, 1) = surplus(i, 1) - Pcharge; % Reduce surplus by the amount charged
            else
                SOC(i+1) = SOC(i); % Battery is full, SOC remains constant
            end
            
            % If surplus remains after charging, sell it
            if surplus(i, 1) > 0
                revenueElectricity(i, 1) = surplus(i, 1) * csell(i, 1) * period;
            else
                revenueElectricity(i, 1) = 0;
            end
        end
        
        % Remaining heating demand (if hydrogen is insufficient)
        remaining_heating_demand = max(0, LH(i, 1) - (H2_heating(i, 1) / electrolyzer_efficiency));
        
        % If hydrogen is insufficient, use natural gas for heating
        if remaining_heating_demand > 0
            NG_heating(i, 1) = remaining_heating_demand / gas_boiler_efficiency;
            costHeating(i, 1) = NG_heating(i, 1) * cbuy(i, 1) * period; % Cost of natural gas for heating
            CO2Heating = NG_heating(i, 1) * co2_ng_per_kWh * period; % CO2 emissions from natural gas
        else
            NG_heating(i, 1) = 0; % No natural gas used if heating demand is met by hydrogen
            costHeating(i, 1) = 0; % No cost incurred
            CO2Heating = 0; % No CO2 emissions
        end
        CO2total(i, 1) = (CO2grid + CO2Heating)/1000; % In KG
    end
    
    % Total cost and CO2 emissions
    totalElectricityCost = sum(purchaseG_Electricity); % Total cost for electricity
    totalHeatingCost = sum(costHeating); % Total cost for heating
    totalRevenue = sum(revenueElectricity); % Total revenue from selling surplus electricity
    totalCost = totalElectricityCost + totalHeatingCost - totalRevenue; % Net total cost
    totalCO2 = sum(CO2total); % Total CO2 emissions
    
    % Define the cost function for the optimization solver
    objectives = [totalCost, totalCO2]; % First objective: cost, second objective: CO2 emissions
end

%% Nonlinear constraints
function [c,ceq] = const(x,T,period,SOCmax,SOCmin)
    Pbat = x(1:T)';
    Estor=zeros(T+1,1);
    Estor(1,1) = 0;
    cBAT = zeros(T,1);
    
    % Battery storage doesn't exceed max capacity
    for i=2:T+1
        Estor(i) = Estor(i-1,1) + (Pbat(i-1,1)*period);
    end
        
    C2 = (SOCmax + SOCmin)/2;
    
    for i=1:T
        cBAT(i,1)=abs(Estor(i+1,1)-C2)-C2;
    end
    
    c= cBAT;
    ceq= [];
end