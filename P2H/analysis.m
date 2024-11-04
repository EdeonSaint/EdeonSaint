function [totalCost, tCO2Gen, Pgrid_Electricity, Pgrid_Hydrogen, H2_heating, NG_heating, Estor] = analysis(x, T, period, cbuy, csell, Pload, heating_demand, Ppv, CEI)
    % Grid cost function
    % Initialize variables
    % Grid cost function
    % Initialize variables
    Pbat = x(1:T)'; % Battery power usage over time

    Pgrid_Electricity = zeros(T,1); % Power drawn from the grid
    surplus = zeros(T+1,1); % Surplus energy (if any)
    purchaseG_Electricity = zeros(T,1); % Electricity purchased from the grid
    costElectricity = zeros(T,1); % Cost of electricity for each time period
    revenueElectricity = zeros(T,1); % Revenue from selling electricity

    Pgrid_Hydrogen = zeros(T,1); % Electricity used for hydrogen production
    H2_heating = zeros(T,1); % Hydrogen used for heating in kWh
    NG_heating = zeros(T,1); % Natural gas used for heating in kWh
    costHeating = zeros(T,1); % Cost of heating for each time period
    CO2Heating = zeros(T,1); % CO2 emissions from heating

    % Define constants
    electrolyzer_efficiency = 0.7; % Efficiency of hydrogen production
    gas_boiler_efficiency = 0.9; % Efficiency of natural gas boiler
    co2_ng_per_kWh = 185; % CO2 emissions from natural gas in kgCO2/kWh

    for i=1:T
        % Calculate the electricity demand and surplus
        Pgrid_Electricity(i,1) = Pload(i,1) - Ppv(i,1) + Pbat(i,1); 

        if Pgrid_Electricity(i,1) < 0
            surplus(i+1,1) = -Pgrid_Electricity(i,1);
            purchaseG_Electricity(i,1) = 0; 
        else
            surplus(i+1,1) = 0;
            purchaseG_Electricity(i,1) = Pgrid_Electricity(i,1);
        end

        costElectricity(i,1) = (cbuy(i,1) * purchaseG_Electricity(i,1)) * period;

        % Allocate surplus energy (implicitly through optimization)
        % Surplus is first used for hydrogen production
        Pload_Hydrogen(i,1) = surplus(i+1,1); % All surplus goes to hydrogen production
        H2_heating(i,1) = Pload_Hydrogen(i,1) * electrolyzer_efficiency; % kWh from hydrogen

        % Calculate revenue from any remaining surplus (if it's better to sell it)
        %revenueElectricity(i,1) = surplus(i+1,1) * csell(i,1) * period;

        % Calculate remaining heating demand
        remaining_heating_demand = max(0, heating_demand(i,1) - H2_heating(i,1));
        
        % Use natural gas if hydrogen is insufficient
        NG_heating(i,1) = remaining_heating_demand / gas_boiler_efficiency;

        % Calculate costs for heating
        costHeating(i,1) = (H2_heating(i,1) * cbuy(i,1)) + (NG_heating(i,1) * 0.0521);

        % Calculate CO2 emissions for heating
        CO2Heating(i,1) = (H2_heating(i,1) * CEI(i,1)) + (NG_heating(i,1) * co2_ng_per_kWh);

        % Update total grid power required and CO2 emissions
        totalGrid(i,1) = Pgrid_Electricity(i,1) + Pload_Hydrogen(i,1);
        if totalGrid(i,1) > 0
            CO2Gen(i,1) = (CEI(i,1) * totalGrid(i,1)) * (period / 1000); % gCO2
        else
            CO2Gen(i,1) = 0;
        end
    end

    % Sum the total costs, revenues, and CO2 emissions
    tCostElectricity = sum(costElectricity);
    tRevenueElectricity = sum(revenueElectricity);
    tCostHeating = sum(costHeating);
    tCO2Gen = sum(CO2Gen) + sum(CO2Heating);

    % Adjust cost function by considering the revenue from selling electricity
    totalCost = tCostElectricity + tCostHeating - tRevenueElectricity;

    Estor=zeros(T+1,1);
    Estor(1,1) = 0;

    for i=2:T+1
        Estor(i) = Estor(i-1,1) + (Pbat(i-1,1)*period);
    end
    % Define the cost function for the optimization solver
    f = [totalCost, tCO2Gen];
end


