% 1. Load the necessary data and structure
load('CP_base_solar')
load('CP_base_wind')
load('CP_costs')
% the CP_costs structure must already have arrays of lists for the
% categories of locations:
% average: middle 45-55% of solar and wind no curtailment LCOA
% topsolar: top 10% of solar and average 45-55% of wind (no curtailment
% LCOA)
% topwind: top 10% of wind and average 45-55% of solar (no curtailment
% LCOA)
% top: top 15% of solar and wind no curtailment LCOA individually
PReq=100;
panelCapital = CP_base_solar.panelCapital; % $/kWp
turbineCapital = CP_base_wind.turbineCapital;
electrolyserCapital = CP_base_solar.electrolyserCapital; % $/GW fed in
ASUCapital = CP_base_solar.ASUCapital; %$/GW fed in
H2compCapital = CP_base_solar.H2compCapital; %$/GW fed in
N2compCapital = CP_base_solar.N2compCapital; %$/GW fed in
HBCapital = CP_base_solar.HBCapital; %dollars for 230 t/day plant
H2storageCapital = CP_base_solar.H2storageCapital; % $/kg
N2storageCapital = CP_base_solar.N2storageCapital; %$/kg
BattstorageCapital = CP_base_solar.BattstorageCapital; %$/kWh
BattpowerCapital = CP_base_solar.BattpowerCapital; %$/GW fed in
panel_OnM = CP_base_solar.panel_OnM; % Dollar per MW per year
turbine_OnM = CP_base_wind.turbine_OnM_perc;
elect_OnM = CP_base_solar.elect_OnM;
HB_OnM = CP_base_solar.HB_OnM; %dollar per year
discountRate = CP_base_solar.discountRate;
opYear = CP_base_solar.opYear; % years
annual=(CP_base_solar.discountRate/(1-(1+CP_base_solar.discountRate)^(-1*CP_base_solar.opYear)));

% 2.1.a Run the optimization for average locations with solar and wind
% combined
CP_costs.average_results = struct; % initialize structure
for i = 1:10
    loc = CP_costs.average_loc(i,1); % extract the location from the list
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc; % get the indeces using the location list that is synced between solar and wind
    panelPowers = CP_base_solar.panelPower(ind); % extract the power profiles of solar energy
    turbinePowers = CP_base_wind.turbinePower(ind); % extract the power profiles of wind energy
    years = CP_base_solar.year(ind); % extract the years
    panelskWp = CP_curt.panelkWp(ind); % extract the optimized panel size with curtailment
    turbineskWp = CP_curt.turbinekWp(ind); % extract the optimized turbine size with curtailment
    clear panel_save turbine_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6 % cycle through the years
        panelPower = panelPowers{y}; % extract the solar power profile for the year
        turbinePower = turbinePowers{y}; % extract the wind power profile for the year
        % combined power profile based on the optimized panel and turbine
        % size
        combPower_perkWp = (panelPower*panelskWp(y) + turbinePower*turbineskWp(y))/(panelskWp(y) + turbineskWp(y));
        % calculate the combined solar and wind kWp size for the case of no
        % curtailment
        combkWpGuess = 1e6; % kWp
        x0 = combkWpGuess;
        inventorySum = @(x) sum((combPower_perkWp*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
        x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
        combkWp_nocurt = x; %store panel kWp 
        % calculate the panel and turbine size individually for the case of
        % no curtailment
        panelkWp_nocurt = combkWp_nocurt * (panelskWp(y)/(turbineskWp(y)+panelskWp(y)));
        turbinekWp_nocurt = combkWp_nocurt * (turbineskWp(y)/(turbineskWp(y)+panelskWp(y)));

        N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
        % run the optimization for 0-99% curtailment, with a constant ratio
        % of turbine size to solar size
        parfor j=1:100
            E_supply = ((turbinekWp_nocurt/(1-(j-1)/100))*turbinePower/1E9) + (panelkWp_nocurt/(1-(j-1)/100))*panelPower/1E9; %GW
            % new panel and turbine size
            panelkWp_new = (panelkWp_nocurt/(1-(j-1)/100));
            turbinekWp_new = (turbinekWp_nocurt/(1-(j-1)/100));
            combProblem=optimproblem('ObjectiveSense','minimize');
            % Create variables
            H2_size=optimvar('H2_size',1,'LowerBound',1); % size of H2 generation relative to the HB
            Batt_size = optimvar('Batt_size',1); %size of the battery power intake capacity relative to HB
            H2_buffer = optimvar('H2_buffer',1,'LowerBound',0); % H2 buffer storage size in GJ
            Batt_buffer = optimvar('Batt_buffer',1,'LowerBound',0); % Battery buffer storage in GJ
            H2_initial = optimvar('H2_initial',1,'LowerBound',0); % initial level of the H2 buffer in GJ
            Batt_initial = optimvar('Batt_initial',1,'LowerBound',0); % initial level of the battery in GJ
            H2_store = optimvar('H2_store',8760,'LowerBound',0); % level of the H2 buffer at each hour over the year in GJ
            Batt_store = optimvar('Batt_store',8760,'LowerBound',0); % level of the battery at each hour over the year in GJ
            H2_extract = optimvar('H2_extract',8760,'LowerBound',0); % amount of energy fed to H2 production at each hour over the year in GW
            Elec_extract = optimvar('Elec_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)
            % Create constraints
            powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be equal to the energy supply
            Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
            Max_Batt = Elec_extract - (27.5/640)*PReq*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity
            H2Balance = optimconstr(8760); %initialize H2 balance of the year
            BattBalance = optimconstr(8760); % initialize battery balance over the year
            %H2 balance at all hours except the first
            H2Balance(2:end) = H2_store(2:end) == H2_store(1:end-1) +...
                (H2_extract(2:end) - PReq*(612.5/640000))*3600; % in GJ
            %H2 balance at the first hour
            H2Balance(1) = H2_store(1) == H2_initial +...
                (H2_extract(1) - PReq*(612.5/640000))*3600; 
            H2finalBalance = H2_store(end,:) == H2_initial; % final energy inventory must equal the initial inventory
            % battery balance at all hours except the first
            BattBalance(2:end) = Batt_store(2:end) == Batt_store(1:end-1) +...
                (Elec_extract(2:end) - PReq*(27.5/640000))*3600; % in GJ
            %battery balance at the first hour
            BattBalance(1) = Batt_store(1) == Batt_initial +...
                (Elec_extract(1)  - PReq*(27.5/640000))*3600; 
            BattfinalBalance = Batt_store(end,:) == Batt_initial; % final energy inventory must equal the initial inventory
            boundedH2 = H2_store <= H2_buffer; % H2 buffer level must be less than the H2 buffer size
            boundedH2Initial = H2_initial <= H2_buffer;
            boundedBatt = Batt_store <= Batt_buffer; % battery buffer level must be less than the battery size
            boundedBattInitial = Batt_initial <= Batt_buffer;
            % define annualized capital cost
            Capex=annual*((panelkWp_new*panelCapital) ...
                +(turbinekWp_new*turbineCapital) ...
                +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
                +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
                +(Batt_size*PReq*BattpowerCapital/1000) ...
                +HBCapital);
            % define operating cost 
            Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(panel_OnM*panelkWp_new/1000) ...
                +(turbine_OnM*turbinekWp_new*turbineCapital*annual) ...
                +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
                +(HBCapital*HB_OnM);
            combProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
                    
            % add constraints to the problem
            combProblem.Constraints.boundedH2 = boundedH2;
            combProblem.Constraints.boundedH2Initial = boundedH2Initial;
            combProblem.Constraints.boundedBatt = boundedBatt;
            combProblem.Constraints.boundedBattInitial = boundedBattInitial;
            combProblem.Constraints.powerBalance = powerBalance;
            combProblem.Constraints.H2Balance = H2Balance;
            combProblem.Constraints.BattBalance = BattBalance;
            combProblem.Constraints.H2finalBalance = H2finalBalance;
            combProblem.Constraints.BattfinalBalance = BattfinalBalance;
            combProblem.Constraints.Max_H2 = Max_H2;
            combProblem.Constraints.Max_Batt = Max_Batt;
            [sol,fval,exitflag,output] = solve(combProblem); % execute solver
            if exitflag==1 
                panelkWp(j,1)=panelkWp_new;
                turbinekWp(j,1)=turbinekWp_new;
                H2SizeGW(j,1)=sol.H2_size*PReq/1000;%GW
                N2SizeGW(j,1)=N2_size*PReq/1000; %GW
                BattSizeGW(j,1)=sol.Batt_size*PReq/1000; %GW
                H2bufferkg(j,1)=sol.H2_buffer*(1/0.2045); %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=sol.Batt_buffer*277.78;
                LCOA(j,1)=fval;
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000);
            
            else
                panelkWp(j,1)=0;
                turbinekWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1}=0;
                error_index(j,1)=i;
            end
       end
       panel_save(:,y) = panelkWp;
       turbine_save(:,y) = turbinekWp;
       H2Size_save(:,y) = H2SizeGW;
       N2Size_save(:,y) = N2SizeGW;
       BattSize_save(:,y) = BattSizeGW;
       H2buffer_save(:,y) = H2bufferkg;
       N2buffer_save(:,y) = N2bufferkg;
       battery_save(:,y) = batterykWh;
       LCOA_save(:,y) = LCOA;
       batteryUse_save{:,y} = batteryUse;
    end
    CP_costs.average_results.panel{i} = panel_save;
    CP_costs.average_results.turbine{i} = turbine_save;
    CP_costs.average_results.H2Size{i} = H2Size_save;
    CP_costs.average_results.N2Size{i} = N2Size_save;
    CP_costs.average_results.BattSize{i} = BattSize_save;
    CP_costs.average_results.H2buffer{i} = H2buffer_save;
    CP_costs.average_results.N2buffer{i} = N2buffer_save;
    CP_costs.average_results.batterykWh{i} = battery_save;
    CP_costs.average_results.LCOA{i} = LCOA_save;
    CP_costs.average_results.batteryUse{i} = batteryUse_save;
end

% 2.1.b Run the optimization for 0-99% curtailment for average locations using only solar energy
CP_costs.average_solaronly = struct;
for i = 1:10
    loc = CP_costs.average_loc(i,1);
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc;
    panelPowers = CP_base_solar.panelPower(ind);
    years = CP_base_solar.year(ind);
    clear panel_save turbine_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6
        panelPower = panelPowers{y};
        solarkWpGuess = 1e6; % kWp
        x0 = solarkWpGuess;
        inventorySum = @(x) sum((panelPower*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
        x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
        panelkWp_nocurt = x; %store panel kWp 
        N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
        parfor j=1:100
            E_supply = (panelkWp_nocurt/(1-(j-1)/100))*panelPower/1E9; %GW
            
            panelkWp_new = (panelkWp_nocurt/(1-(j-1)/100));
            solarProblem=optimproblem('ObjectiveSense','minimize');
            % Create variables
            H2_size=optimvar('H2_size',1,'LowerBound',1); % size of H2 generation relative to the HB
            Batt_size = optimvar('Batt_size',1); %size of the battery power intake capacity relative to HB
            H2_buffer = optimvar('H2_buffer',1,'LowerBound',0); % H2 buffer storage size in GJ
            Batt_buffer = optimvar('Batt_buffer',1,'LowerBound',0); % Battery buffer storage in GJ
            H2_initial = optimvar('H2_initial',1,'LowerBound',0); % initial level of the H2 buffer in GJ
            Batt_initial = optimvar('Batt_initial',1,'LowerBound',0); % initial level of the battery in GJ
            H2_store = optimvar('H2_store',8760,'LowerBound',0); % level of the H2 buffer at each hour over the year in GJ
            Batt_store = optimvar('Batt_store',8760,'LowerBound',0); % level of the battery at each hour over the year in GJ
            H2_extract = optimvar('H2_extract',8760,'LowerBound',0); % amount of energy fed to H2 production at each hour over the year in GW
            Elec_extract = optimvar('Elec_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)
            % Create constraints
            powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be equal to the energy supply
            Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
            Max_Batt = Elec_extract - (27.5/640)*PReq*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity
            H2Balance = optimconstr(8760); %initialize H2 balance of the year
            BattBalance = optimconstr(8760); % initialize battery balance over the year
            %H2 balance at all hours except the first
            H2Balance(2:end) = H2_store(2:end) == H2_store(1:end-1) +...
                (H2_extract(2:end) - PReq*(612.5/640000))*3600; % in GJ
            %H2 balance at the first hour
            H2Balance(1) = H2_store(1) == H2_initial +...
                (H2_extract(1) - PReq*(612.5/640000))*3600; 
            H2finalBalance = H2_store(end,:) == H2_initial; % final energy inventory must equal the initial inventory
            % battery balance at all hours except the first
            BattBalance(2:end) = Batt_store(2:end) == Batt_store(1:end-1) +...
                (Elec_extract(2:end) - PReq*(27.5/640000))*3600; % in GJ
            %battery balance at the first hour
            BattBalance(1) = Batt_store(1) == Batt_initial +...
                (Elec_extract(1)  - PReq*(27.5/640000))*3600; 
            BattfinalBalance = Batt_store(end,:) == Batt_initial; % final energy inventory must equal the initial inventory
            boundedH2 = H2_store <= H2_buffer; % H2 buffer level must be less than the H2 buffer size
            boundedH2Initial = H2_initial <= H2_buffer;
            boundedBatt = Batt_store <= Batt_buffer; % battery buffer level must be less than the battery size
            boundedBattInitial = Batt_initial <= Batt_buffer; 
            % define annualized capital cost
            Capex=annual*((panelkWp_new*panelCapital) ...
                +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
                +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
                +(Batt_size*PReq*BattpowerCapital/1000) ...
                +HBCapital);
            % define operating cost 
            Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(panel_OnM*panelkWp_new/1000) ...
                +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
                +(HBCapital*HB_OnM);
            solarProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
                    
            % add constraints to the problem
            solarProblem.Constraints.boundedH2 = boundedH2;
            solarProblem.Constraints.boundedH2Initial = boundedH2Initial;
            solarProblem.Constraints.boundedBatt = boundedBatt;
            solarProblem.Constraints.boundedBattInitial = boundedBattInitial;
            solarProblem.Constraints.powerBalance = powerBalance;
            solarProblem.Constraints.H2Balance = H2Balance;
            solarProblem.Constraints.BattBalance = BattBalance;
            solarProblem.Constraints.H2finalBalance = H2finalBalance;
            solarProblem.Constraints.BattfinalBalance = BattfinalBalance;
            solarProblem.Constraints.Max_H2 = Max_H2;
            solarProblem.Constraints.Max_Batt = Max_Batt;
            [sol,fval,exitflag,output] = solve(solarProblem); % execute solver
            if exitflag==1 
                panelkWp(j,1)=panelkWp_new;
                H2SizeGW(j,1)=sol.H2_size*PReq/1000;%GW
                N2SizeGW(j,1)=N2_size*PReq/1000; %GW
                BattSizeGW(j,1)=sol.Batt_size*PReq/1000; %GW
                H2bufferkg(j,1)=sol.H2_buffer*(1/0.2045); %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=sol.Batt_buffer*277.78;
                LCOA(j,1)=fval;
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000);
            
            else
                panelkWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1}=0;
                error_index(j,1)=i;
            end
       end
       panel_save(:,y) = panelkWp;
       H2Size_save(:,y) = H2SizeGW;
       N2Size_save(:,y) = N2SizeGW;
       BattSize_save(:,y) = BattSizeGW;
       H2buffer_save(:,y) = H2bufferkg;
       N2buffer_save(:,y) = N2bufferkg;
       battery_save(:,y) = batterykWh;
       LCOA_save(:,y) = LCOA;
       batteryUse_save{:,y} = batteryUse;
    end
    CP_costs.average_solaronly.panel{i} = panel_save;
    CP_costs.average_solaronly.H2Size{i} = H2Size_save;
    CP_costs.average_solaronly.N2Size{i} = N2Size_save;
    CP_costs.average_solaronly.BattSize{i} = BattSize_save;
    CP_costs.average_solaronly.H2buffer{i} = H2buffer_save;
    CP_costs.average_solaronly.N2buffer{i} = N2buffer_save;
    CP_costs.average_solaronly.batterykWh{i} = battery_save;
    CP_costs.average_solaronly.LCOA{i} = LCOA_save;
end

% 2.i.c Run the optimization for 0-99% curtailment for average locations
% using only wind energy
CP_costs.average_windonly = struct;
for i = 1:10
    loc = CP_costs.average_loc(i,1);
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc;
    turbinePowers = CP_base_wind.turbinePower(ind);
    years = CP_base_solar.year(ind);
    clear panel_save turbine_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6
        turbinePower = turbinePowers{y};
        turbinekWpGuess = 1e6; % kWp
        x0 = turbinekWpGuess;
        inventorySum = @(x) sum((turbinePower*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
        x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
        turbinekWp_nocurt = x; %store panel kWp 
        N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
        parfor j=1:100
            E_supply = ((turbinekWp_nocurt/(1-(j-1)/100))*turbinePower/1E9); %GW 
            turbinekWp_new = (turbinekWp_nocurt/(1-(j-1)/100));
            windProblem=optimproblem('ObjectiveSense','minimize');
            % Create variables
            H2_size=optimvar('H2_size',1,'LowerBound',1); % size of H2 generation relative to the HB
            Batt_size = optimvar('Batt_size',1); %size of the battery power intake capacity relative to HB
            H2_buffer = optimvar('H2_buffer',1,'LowerBound',0); % H2 buffer storage size in GJ
            Batt_buffer = optimvar('Batt_buffer',1,'LowerBound',0); % Battery buffer storage in GJ
            H2_initial = optimvar('H2_initial',1,'LowerBound',0); % initial level of the H2 buffer in GJ
            Batt_initial = optimvar('Batt_initial',1,'LowerBound',0); % initial level of the battery in GJ
            H2_store = optimvar('H2_store',8760,'LowerBound',0); % level of the H2 buffer at each hour over the year in GJ
            Batt_store = optimvar('Batt_store',8760,'LowerBound',0); % level of the battery at each hour over the year in GJ
            H2_extract = optimvar('H2_extract',8760,'LowerBound',0); % amount of energy fed to H2 production at each hour over the year in GW
            Elec_extract = optimvar('Elec_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)
            % Create constraints
            powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be equal to the energy supply
            Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
            Max_Batt = Elec_extract - (27.5/640)*PReq*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity
            H2Balance = optimconstr(8760); %initialize H2 balance of the year
            BattBalance = optimconstr(8760); % initialize battery balance over the year
            %H2 balance at all hours except the first
            H2Balance(2:end) = H2_store(2:end) == H2_store(1:end-1) +...
                (H2_extract(2:end) - PReq*(612.5/640000))*3600; % in GJ
            %H2 balance at the first hour
            H2Balance(1) = H2_store(1) == H2_initial +...
                (H2_extract(1) - PReq*(612.5/640000))*3600; 
            H2finalBalance = H2_store(end,:) == H2_initial; % final energy inventory must equal the initial inventory
            % battery balance at all hours except the first
            BattBalance(2:end) = Batt_store(2:end) == Batt_store(1:end-1) +...
                (Elec_extract(2:end) - PReq*(27.5/640000))*3600; % in GJ
            %battery balance at the first hour
            BattBalance(1) = Batt_store(1) == Batt_initial +...
                (Elec_extract(1)  - PReq*(27.5/640000))*3600; 
            BattfinalBalance = Batt_store(end,:) == Batt_initial; % final energy inventory must equal the initial inventory
            boundedH2 = H2_store <= H2_buffer; % H2 buffer level must be less than the H2 buffer size
            boundedH2Initial = H2_initial <= H2_buffer;
            boundedBatt = Batt_store <= Batt_buffer; % battery buffer level must be less than the battery size
            boundedBattInitial = Batt_initial <= Batt_buffer; 
            % define annualized capital cost
            Capex=annual*((turbinekWp_new*turbineCapital) ...
                +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
                +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
                +(Batt_size*PReq*BattpowerCapital/1000) ...
                +HBCapital);
            % define operating cost 
            Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(turbine_OnM*turbinekWp_new*turbineCapital*annual) ...
                +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
                +(HBCapital*HB_OnM);
            windProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
                    
            % add constraints to the problem
            windProblem.Constraints.boundedH2 = boundedH2;
            windProblem.Constraints.boundedH2Initial = boundedH2Initial;
            windProblem.Constraints.boundedBatt = boundedBatt;
            windProblem.Constraints.boundedBattInitial = boundedBattInitial;
            windProblem.Constraints.powerBalance = powerBalance;
            windProblem.Constraints.H2Balance = H2Balance;
            windProblem.Constraints.BattBalance = BattBalance;
            windProblem.Constraints.H2finalBalance = H2finalBalance;
            windProblem.Constraints.BattfinalBalance = BattfinalBalance;
            windProblem.Constraints.Max_H2 = Max_H2;
            windProblem.Constraints.Max_Batt = Max_Batt;
            [sol,fval,exitflag,output] = solve(windProblem); % execute solver
            if exitflag==1 
                turbinekWp(j,1)=turbinekWp_new;
                H2SizeGW(j,1)=sol.H2_size*PReq/1000;%GW
                N2SizeGW(j,1)=N2_size*PReq/1000; %GW
                BattSizeGW(j,1)=sol.Batt_size*PReq/1000; %GW
                H2bufferkg(j,1)=sol.H2_buffer*(1/0.2045); %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=sol.Batt_buffer*277.78;
                LCOA(j,1)=fval;
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000);
            
            else
                turbinekWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1}=0;
                error_index(j,1)=i;
            end
        end
       turbine_save(:,y) = turbinekWp;
       H2Size_save(:,y) = H2SizeGW;
       N2Size_save(:,y) = N2SizeGW;
       BattSize_save(:,y) = BattSizeGW;
       H2buffer_save(:,y) = H2bufferkg;
       N2buffer_save(:,y) = N2bufferkg;
       battery_save(:,y) = batterykWh;
       LCOA_save(:,y) = LCOA;
       batteryUse_save{:,y} = batteryUse;
    end
    CP_costs.average_windonly.turbine{i} = turbine_save;
    CP_costs.average_windonly.H2Size{i} = H2Size_save;
    CP_costs.average_windonly.N2Size{i} = N2Size_save;
    CP_costs.average_windonly.BattSize{i} = BattSize_save;
    CP_costs.average_windonly.H2buffer{i} = H2buffer_save;
    CP_costs.average_windonly.N2buffer{i} = N2buffer_save;
    CP_costs.average_windonly.batterykWh{i} = battery_save;
    CP_costs.average_windonly.LCOA{i} = LCOA_save;
end
save('CP_costs','CP_costs')

% 2.ii.a Run the optimization for top solar locations with solar and wind
% combined
CP_costs.topsolar_results = struct; % initialize structure
for i = 1:10
    loc = CP_costs.topsolar_loc(i,1); % extract the location from the list
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc; % get the indeces using the location list that is synced between solar and wind
    panelPowers = CP_base_solar.panelPower(ind); % extract the power profiles of solar energy
    turbinePowers = CP_base_wind.turbinePower(ind); % extract the power profiles of wind energy
    years = CP_base_solar.year(ind); % extract the years
    panelskWp = CP_curt.panelkWp(ind); % extract the optimized panel size with curtailment
    turbineskWp = CP_curt.turbinekWp(ind); % extract the optimized turbine size with curtailment
    clear panel_save turbine_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6 % cycle through the years
        panelPower = panelPowers{y}; % extract the solar power profile for the year
        turbinePower = turbinePowers{y}; % extract the wind power profile for the year
        % combined power profile based on the optimized panel and turbine
        % size
        combPower_perkWp = (panelPower*panelskWp(y) + turbinePower*turbineskWp(y))/(panelskWp(y) + turbineskWp(y));
        % calculate the combined solar and wind kWp size for the case of no
        % curtailment
        combkWpGuess = 1e6; % kWp
        x0 = combkWpGuess;
        inventorySum = @(x) sum((combPower_perkWp*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
        x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
        combkWp_nocurt = x; %store panel kWp 
        % calculate the panel and turbine size individually for the case of
        % no curtailment
        panelkWp_nocurt = combkWp_nocurt * (panelskWp(y)/(turbineskWp(y)+panelskWp(y)));
        turbinekWp_nocurt = combkWp_nocurt * (turbineskWp(y)/(turbineskWp(y)+panelskWp(y)));

        N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
        % run the optimization for 0-99% curtailment, with a constant ratio
        % of turbine size to solar size
        parfor j=1:100
            E_supply = ((turbinekWp_nocurt/(1-(j-1)/100))*turbinePower/1E9) + (panelkWp_nocurt/(1-(j-1)/100))*panelPower/1E9; %GW
            % new panel and turbine size
            panelkWp_new = (panelkWp_nocurt/(1-(j-1)/100));
            turbinekWp_new = (turbinekWp_nocurt/(1-(j-1)/100));
            combProblem=optimproblem('ObjectiveSense','minimize');
            % Create variables
            H2_size=optimvar('H2_size',1,'LowerBound',1); % size of H2 generation relative to the HB
            Batt_size = optimvar('Batt_size',1); %size of the battery power intake capacity relative to HB
            H2_buffer = optimvar('H2_buffer',1,'LowerBound',0); % H2 buffer storage size in GJ
            Batt_buffer = optimvar('Batt_buffer',1,'LowerBound',0); % Battery buffer storage in GJ
            H2_initial = optimvar('H2_initial',1,'LowerBound',0); % initial level of the H2 buffer in GJ
            Batt_initial = optimvar('Batt_initial',1,'LowerBound',0); % initial level of the battery in GJ
            H2_store = optimvar('H2_store',8760,'LowerBound',0); % level of the H2 buffer at each hour over the year in GJ
            Batt_store = optimvar('Batt_store',8760,'LowerBound',0); % level of the battery at each hour over the year in GJ
            H2_extract = optimvar('H2_extract',8760,'LowerBound',0); % amount of energy fed to H2 production at each hour over the year in GW
            Elec_extract = optimvar('Elec_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)
            % Create constraints
            powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be equal to the energy supply
            Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
            Max_Batt = Elec_extract - (27.5/640)*PReq*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity
            H2Balance = optimconstr(8760); %initialize H2 balance of the year
            BattBalance = optimconstr(8760); % initialize battery balance over the year
            %H2 balance at all hours except the first
            H2Balance(2:end) = H2_store(2:end) == H2_store(1:end-1) +...
                (H2_extract(2:end) - PReq*(612.5/640000))*3600; % in GJ
            %H2 balance at the first hour
            H2Balance(1) = H2_store(1) == H2_initial +...
                (H2_extract(1) - PReq*(612.5/640000))*3600; 
            H2finalBalance = H2_store(end,:) == H2_initial; % final energy inventory must equal the initial inventory
            % battery balance at all hours except the first
            BattBalance(2:end) = Batt_store(2:end) == Batt_store(1:end-1) +...
                (Elec_extract(2:end) - PReq*(27.5/640000))*3600; % in GJ
            %battery balance at the first hour
            BattBalance(1) = Batt_store(1) == Batt_initial +...
                (Elec_extract(1)  - PReq*(27.5/640000))*3600; 
            BattfinalBalance = Batt_store(end,:) == Batt_initial; % final energy inventory must equal the initial inventory
            boundedH2 = H2_store <= H2_buffer; % H2 buffer level must be less than the H2 buffer size
            boundedH2Initial = H2_initial <= H2_buffer;
            boundedBatt = Batt_store <= Batt_buffer; % battery buffer level must be less than the battery size
            boundedBattInitial = Batt_initial <= Batt_buffer;
            % define annualized capital cost
            Capex=annual*((panelkWp_new*panelCapital) ...
                +(turbinekWp_new*turbineCapital) ...
                +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
                +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
                +(Batt_size*PReq*BattpowerCapital/1000) ...
                +HBCapital);
            % define operating cost 
            Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(panel_OnM*panelkWp_new/1000) ...
                +(turbine_OnM*turbinekWp_new*turbineCapital*annual) ...
                +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
                +(HBCapital*HB_OnM);
            combProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
                    
            % add constraints to the problem
            combProblem.Constraints.boundedH2 = boundedH2;
            combProblem.Constraints.boundedH2Initial = boundedH2Initial;
            combProblem.Constraints.boundedBatt = boundedBatt;
            combProblem.Constraints.boundedBattInitial = boundedBattInitial;
            combProblem.Constraints.powerBalance = powerBalance;
            combProblem.Constraints.H2Balance = H2Balance;
            combProblem.Constraints.BattBalance = BattBalance;
            combProblem.Constraints.H2finalBalance = H2finalBalance;
            combProblem.Constraints.BattfinalBalance = BattfinalBalance;
            combProblem.Constraints.Max_H2 = Max_H2;
            combProblem.Constraints.Max_Batt = Max_Batt;
            [sol,fval,exitflag,output] = solve(combProblem); % execute solver
            if exitflag==1 
                panelkWp(j,1)=panelkWp_new;
                turbinekWp(j,1)=turbinekWp_new;
                H2SizeGW(j,1)=sol.H2_size*PReq/1000;%GW
                N2SizeGW(j,1)=N2_size*PReq/1000; %GW
                BattSizeGW(j,1)=sol.Batt_size*PReq/1000; %GW
                H2bufferkg(j,1)=sol.H2_buffer*(1/0.2045); %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=sol.Batt_buffer*277.78;
                LCOA(j,1)=fval;
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000);
            
            else
                panelkWp(j,1)=0;
                turbinekWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1}=0;
                error_index(j,1)=i;
            end
       end
       panel_save(:,y) = panelkWp;
       turbine_save(:,y) = turbinekWp;
       H2Size_save(:,y) = H2SizeGW;
       N2Size_save(:,y) = N2SizeGW;
       BattSize_save(:,y) = BattSizeGW;
       H2buffer_save(:,y) = H2bufferkg;
       N2buffer_save(:,y) = N2bufferkg;
       battery_save(:,y) = batterykWh;
       LCOA_save(:,y) = LCOA;
       batteryUse_save{:,y} = batteryUse;
    end
    CP_costs.topsolar_results.panel{i} = panel_save;
    CP_costs.topsolar_results.turbine{i} = turbine_save;
    CP_costs.topsolar_results.H2Size{i} = H2Size_save;
    CP_costs.topsolar_results.N2Size{i} = N2Size_save;
    CP_costs.topsolar_results.BattSize{i} = BattSize_save;
    CP_costs.topsolar_results.H2buffer{i} = H2buffer_save;
    CP_costs.topsolar_results.N2buffer{i} = N2buffer_save;
    CP_costs.topsolar_results.batterykWh{i} = battery_save;
    CP_costs.topsolar_results.LCOA{i} = LCOA_save;
    CP_costs.topsolar_results.batteryUse{i} = batteryUse_save;
end

% 2.ii.b Run the optimization for 0-99% curtailment for top solar locations using only solar energy
CP_costs.topsolar_solaronly = struct;
for i = 1:10
    loc = CP_costs.topsolar_loc(i,1);
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc;
    panelPowers = CP_base_solar.panelPower(ind);
    years = CP_base_solar.year(ind);
    clear panel_save turbine_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6
        panelPower = panelPowers{y};
        solarkWpGuess = 1e6; % kWp
        x0 = solarkWpGuess;
        inventorySum = @(x) sum((panelPower*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
        x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
        panelkWp_nocurt = x; %store panel kWp 
        N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
        parfor j=1:100
            E_supply = (panelkWp_nocurt/(1-(j-1)/100))*panelPower/1E9; %GW
            
            panelkWp_new = (panelkWp_nocurt/(1-(j-1)/100));
            solarProblem=optimproblem('ObjectiveSense','minimize');
            % Create variables
            H2_size=optimvar('H2_size',1,'LowerBound',1); % size of H2 generation relative to the HB
            Batt_size = optimvar('Batt_size',1); %size of the battery power intake capacity relative to HB
            H2_buffer = optimvar('H2_buffer',1,'LowerBound',0); % H2 buffer storage size in GJ
            Batt_buffer = optimvar('Batt_buffer',1,'LowerBound',0); % Battery buffer storage in GJ
            H2_initial = optimvar('H2_initial',1,'LowerBound',0); % initial level of the H2 buffer in GJ
            Batt_initial = optimvar('Batt_initial',1,'LowerBound',0); % initial level of the battery in GJ
            H2_store = optimvar('H2_store',8760,'LowerBound',0); % level of the H2 buffer at each hour over the year in GJ
            Batt_store = optimvar('Batt_store',8760,'LowerBound',0); % level of the battery at each hour over the year in GJ
            H2_extract = optimvar('H2_extract',8760,'LowerBound',0); % amount of energy fed to H2 production at each hour over the year in GW
            Elec_extract = optimvar('Elec_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)
            % Create constraints
            powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be equal to the energy supply
            Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
            Max_Batt = Elec_extract - (27.5/640)*PReq*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity
            H2Balance = optimconstr(8760); %initialize H2 balance of the year
            BattBalance = optimconstr(8760); % initialize battery balance over the year
            %H2 balance at all hours except the first
            H2Balance(2:end) = H2_store(2:end) == H2_store(1:end-1) +...
                (H2_extract(2:end) - PReq*(612.5/640000))*3600; % in GJ
            %H2 balance at the first hour
            H2Balance(1) = H2_store(1) == H2_initial +...
                (H2_extract(1) - PReq*(612.5/640000))*3600; 
            H2finalBalance = H2_store(end,:) == H2_initial; % final energy inventory must equal the initial inventory
            % battery balance at all hours except the first
            BattBalance(2:end) = Batt_store(2:end) == Batt_store(1:end-1) +...
                (Elec_extract(2:end) - PReq*(27.5/640000))*3600; % in GJ
            %battery balance at the first hour
            BattBalance(1) = Batt_store(1) == Batt_initial +...
                (Elec_extract(1)  - PReq*(27.5/640000))*3600; 
            BattfinalBalance = Batt_store(end,:) == Batt_initial; % final energy inventory must equal the initial inventory
            boundedH2 = H2_store <= H2_buffer; % H2 buffer level must be less than the H2 buffer size
            boundedH2Initial = H2_initial <= H2_buffer;
            boundedBatt = Batt_store <= Batt_buffer; % battery buffer level must be less than the battery size
            boundedBattInitial = Batt_initial <= Batt_buffer; 
            % define annualized capital cost
            Capex=annual*((panelkWp_new*panelCapital) ...
                +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
                +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
                +(Batt_size*PReq*BattpowerCapital/1000) ...
                +HBCapital);
            % define operating cost 
            Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(panel_OnM*panelkWp_new/1000) ...
                +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
                +(HBCapital*HB_OnM);
            solarProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
                    
            % add constraints to the problem
            solarProblem.Constraints.boundedH2 = boundedH2;
            solarProblem.Constraints.boundedH2Initial = boundedH2Initial;
            solarProblem.Constraints.boundedBatt = boundedBatt;
            solarProblem.Constraints.boundedBattInitial = boundedBattInitial;
            solarProblem.Constraints.powerBalance = powerBalance;
            solarProblem.Constraints.H2Balance = H2Balance;
            solarProblem.Constraints.BattBalance = BattBalance;
            solarProblem.Constraints.H2finalBalance = H2finalBalance;
            solarProblem.Constraints.BattfinalBalance = BattfinalBalance;
            solarProblem.Constraints.Max_H2 = Max_H2;
            solarProblem.Constraints.Max_Batt = Max_Batt;
            [sol,fval,exitflag,output] = solve(solarProblem); % execute solver
            if exitflag==1 
                panelkWp(j,1)=panelkWp_new;
                H2SizeGW(j,1)=sol.H2_size*PReq/1000;%GW
                N2SizeGW(j,1)=N2_size*PReq/1000; %GW
                BattSizeGW(j,1)=sol.Batt_size*PReq/1000; %GW
                H2bufferkg(j,1)=sol.H2_buffer*(1/0.2045); %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=sol.Batt_buffer*277.78;
                LCOA(j,1)=fval;
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000);
            
            else
                panelkWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1}=0;
                error_index(j,1)=i;
            end
       end
       panel_save(:,y) = panelkWp;
       H2Size_save(:,y) = H2SizeGW;
       N2Size_save(:,y) = N2SizeGW;
       BattSize_save(:,y) = BattSizeGW;
       H2buffer_save(:,y) = H2bufferkg;
       N2buffer_save(:,y) = N2bufferkg;
       battery_save(:,y) = batterykWh;
       LCOA_save(:,y) = LCOA;
       batteryUse_save{:,y} = batteryUse;
    end
    CP_costs.topsolar_solaronly.panel{i} = panel_save;
    CP_costs.topsolar_solaronly.H2Size{i} = H2Size_save;
    CP_costs.topsolar_solaronly.N2Size{i} = N2Size_save;
    CP_costs.topsolar_solaronly.BattSize{i} = BattSize_save;
    CP_costs.topsolar_solaronly.H2buffer{i} = H2buffer_save;
    CP_costs.topsolar_solaronly.N2buffer{i} = N2buffer_save;
    CP_costs.topsolar_solaronly.batterykWh{i} = battery_save;
    CP_costs.topsolar_solaronly.LCOA{i} = LCOA_save;
end

% 2.ii.c Run the optimization for 0-99% curtailment for top solar locations
% using only wind energy
CP_costs.topsolar_windonly = struct;
for i = 1:10
    loc = CP_costs.topsolar_loc(i,1);
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc;
    turbinePowers = CP_base_wind.turbinePower(ind);
    years = CP_base_solar.year(ind);
    clear panel_save turbine_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6
        turbinePower = turbinePowers{y};
        turbinekWpGuess = 1e6; % kWp
        x0 = turbinekWpGuess;
        inventorySum = @(x) sum((turbinePower*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
        x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
        turbinekWp_nocurt = x; %store panel kWp 
        N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
        parfor j=1:100
            E_supply = ((turbinekWp_nocurt/(1-(j-1)/100))*turbinePower/1E9); %GW 
            turbinekWp_new = (turbinekWp_nocurt/(1-(j-1)/100));
            windProblem=optimproblem('ObjectiveSense','minimize');
            % Create variables
            H2_size=optimvar('H2_size',1,'LowerBound',1); % size of H2 generation relative to the HB
            Batt_size = optimvar('Batt_size',1); %size of the battery power intake capacity relative to HB
            H2_buffer = optimvar('H2_buffer',1,'LowerBound',0); % H2 buffer storage size in GJ
            Batt_buffer = optimvar('Batt_buffer',1,'LowerBound',0); % Battery buffer storage in GJ
            H2_initial = optimvar('H2_initial',1,'LowerBound',0); % initial level of the H2 buffer in GJ
            Batt_initial = optimvar('Batt_initial',1,'LowerBound',0); % initial level of the battery in GJ
            H2_store = optimvar('H2_store',8760,'LowerBound',0); % level of the H2 buffer at each hour over the year in GJ
            Batt_store = optimvar('Batt_store',8760,'LowerBound',0); % level of the battery at each hour over the year in GJ
            H2_extract = optimvar('H2_extract',8760,'LowerBound',0); % amount of energy fed to H2 production at each hour over the year in GW
            Elec_extract = optimvar('Elec_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)
            % Create constraints
            powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be equal to the energy supply
            Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
            Max_Batt = Elec_extract - (27.5/640)*PReq*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity
            H2Balance = optimconstr(8760); %initialize H2 balance of the year
            BattBalance = optimconstr(8760); % initialize battery balance over the year
            %H2 balance at all hours except the first
            H2Balance(2:end) = H2_store(2:end) == H2_store(1:end-1) +...
                (H2_extract(2:end) - PReq*(612.5/640000))*3600; % in GJ
            %H2 balance at the first hour
            H2Balance(1) = H2_store(1) == H2_initial +...
                (H2_extract(1) - PReq*(612.5/640000))*3600; 
            H2finalBalance = H2_store(end,:) == H2_initial; % final energy inventory must equal the initial inventory
            % battery balance at all hours except the first
            BattBalance(2:end) = Batt_store(2:end) == Batt_store(1:end-1) +...
                (Elec_extract(2:end) - PReq*(27.5/640000))*3600; % in GJ
            %battery balance at the first hour
            BattBalance(1) = Batt_store(1) == Batt_initial +...
                (Elec_extract(1)  - PReq*(27.5/640000))*3600; 
            BattfinalBalance = Batt_store(end,:) == Batt_initial; % final energy inventory must equal the initial inventory
            boundedH2 = H2_store <= H2_buffer; % H2 buffer level must be less than the H2 buffer size
            boundedH2Initial = H2_initial <= H2_buffer;
            boundedBatt = Batt_store <= Batt_buffer; % battery buffer level must be less than the battery size
            boundedBattInitial = Batt_initial <= Batt_buffer; 
            % define annualized capital cost
            Capex=annual*((turbinekWp_new*turbineCapital) ...
                +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
                +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
                +(Batt_size*PReq*BattpowerCapital/1000) ...
                +HBCapital);
            % define operating cost 
            Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(turbine_OnM*turbinekWp_new*turbineCapital*annual) ...
                +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
                +(HBCapital*HB_OnM);
            windProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
                    
            % add constraints to the problem
            windProblem.Constraints.boundedH2 = boundedH2;
            windProblem.Constraints.boundedH2Initial = boundedH2Initial;
            windProblem.Constraints.boundedBatt = boundedBatt;
            windProblem.Constraints.boundedBattInitial = boundedBattInitial;
            windProblem.Constraints.powerBalance = powerBalance;
            windProblem.Constraints.H2Balance = H2Balance;
            windProblem.Constraints.BattBalance = BattBalance;
            windProblem.Constraints.H2finalBalance = H2finalBalance;
            windProblem.Constraints.BattfinalBalance = BattfinalBalance;
            windProblem.Constraints.Max_H2 = Max_H2;
            windProblem.Constraints.Max_Batt = Max_Batt;
            [sol,fval,exitflag,output] = solve(windProblem); % execute solver
            if exitflag==1 
                turbinekWp(j,1)=turbinekWp_new;
                H2SizeGW(j,1)=sol.H2_size*PReq/1000;%GW
                N2SizeGW(j,1)=N2_size*PReq/1000; %GW
                BattSizeGW(j,1)=sol.Batt_size*PReq/1000; %GW
                H2bufferkg(j,1)=sol.H2_buffer*(1/0.2045); %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=sol.Batt_buffer*277.78;
                LCOA(j,1)=fval;
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000);
            
            else
                turbinekWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1}=0;
                error_index(j,1)=i;
            end
        end
       turbine_save(:,y) = turbinekWp;
       H2Size_save(:,y) = H2SizeGW;
       N2Size_save(:,y) = N2SizeGW;
       BattSize_save(:,y) = BattSizeGW;
       H2buffer_save(:,y) = H2bufferkg;
       N2buffer_save(:,y) = N2bufferkg;
       battery_save(:,y) = batterykWh;
       LCOA_save(:,y) = LCOA;
       batteryUse_save{:,y} = batteryUse;
    end
    CP_costs.topsolar_windonly.turbine{i} = turbine_save;
    CP_costs.topsolar_windonly.H2Size{i} = H2Size_save;
    CP_costs.topsolar_windonly.N2Size{i} = N2Size_save;
    CP_costs.topsolar_windonly.BattSize{i} = BattSize_save;
    CP_costs.topsolar_windonly.H2buffer{i} = H2buffer_save;
    CP_costs.topsolar_windonly.N2buffer{i} = N2buffer_save;
    CP_costs.topsolar_windonly.batterykWh{i} = battery_save;
    CP_costs.topsolar_windonly.LCOA{i} = LCOA_save;
end
save('CP_costs','CP_costs')

% 2.iii.a Run the optimization for top wind locations with solar and wind
% combined
CP_costs.topwind_results = struct; % initialize structure
for i = 1:10
    loc = CP_costs.topwind_loc(i,1); % extract the location from the list
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc; % get the indeces using the location list that is synced between solar and wind
    panelPowers = CP_base_solar.panelPower(ind); % extract the power profiles of solar energy
    turbinePowers = CP_base_wind.turbinePower(ind); % extract the power profiles of wind energy
    years = CP_base_solar.year(ind); % extract the years
    panelskWp = CP_curt.panelkWp(ind); % extract the optimized panel size with curtailment
    turbineskWp = CP_curt.turbinekWp(ind); % extract the optimized turbine size with curtailment
    clear panel_save turbine_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6 % cycle through the years
        panelPower = panelPowers{y}; % extract the solar power profile for the year
        turbinePower = turbinePowers{y}; % extract the wind power profile for the year
        % combined power profile based on the optimized panel and turbine
        % size
        combPower_perkWp = (panelPower*panelskWp(y) + turbinePower*turbineskWp(y))/(panelskWp(y) + turbineskWp(y));
        % calculate the combined solar and wind kWp size for the case of no
        % curtailment
        combkWpGuess = 1e6; % kWp
        x0 = combkWpGuess;
        inventorySum = @(x) sum((combPower_perkWp*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
        x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
        combkWp_nocurt = x; %store panel kWp 
        % calculate the panel and turbine size individually for the case of
        % no curtailment
        panelkWp_nocurt = combkWp_nocurt * (panelskWp(y)/(turbineskWp(y)+panelskWp(y)));
        turbinekWp_nocurt = combkWp_nocurt * (turbineskWp(y)/(turbineskWp(y)+panelskWp(y)));

        N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
        % run the optimization for 0-99% curtailment, with a constant ratio
        % of turbine size to solar size
        parfor j=1:100
            E_supply = ((turbinekWp_nocurt/(1-(j-1)/100))*turbinePower/1E9) + (panelkWp_nocurt/(1-(j-1)/100))*panelPower/1E9; %GW
            % new panel and turbine size
            panelkWp_new = (panelkWp_nocurt/(1-(j-1)/100));
            turbinekWp_new = (turbinekWp_nocurt/(1-(j-1)/100));
            combProblem=optimproblem('ObjectiveSense','minimize');
            % Create variables
            H2_size=optimvar('H2_size',1,'LowerBound',1); % size of H2 generation relative to the HB
            Batt_size = optimvar('Batt_size',1); %size of the battery power intake capacity relative to HB
            H2_buffer = optimvar('H2_buffer',1,'LowerBound',0); % H2 buffer storage size in GJ
            Batt_buffer = optimvar('Batt_buffer',1,'LowerBound',0); % Battery buffer storage in GJ
            H2_initial = optimvar('H2_initial',1,'LowerBound',0); % initial level of the H2 buffer in GJ
            Batt_initial = optimvar('Batt_initial',1,'LowerBound',0); % initial level of the battery in GJ
            H2_store = optimvar('H2_store',8760,'LowerBound',0); % level of the H2 buffer at each hour over the year in GJ
            Batt_store = optimvar('Batt_store',8760,'LowerBound',0); % level of the battery at each hour over the year in GJ
            H2_extract = optimvar('H2_extract',8760,'LowerBound',0); % amount of energy fed to H2 production at each hour over the year in GW
            Elec_extract = optimvar('Elec_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)
            % Create constraints
            powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be equal to the energy supply
            Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
            Max_Batt = Elec_extract - (27.5/640)*PReq*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity
            H2Balance = optimconstr(8760); %initialize H2 balance of the year
            BattBalance = optimconstr(8760); % initialize battery balance over the year
            %H2 balance at all hours except the first
            H2Balance(2:end) = H2_store(2:end) == H2_store(1:end-1) +...
                (H2_extract(2:end) - PReq*(612.5/640000))*3600; % in GJ
            %H2 balance at the first hour
            H2Balance(1) = H2_store(1) == H2_initial +...
                (H2_extract(1) - PReq*(612.5/640000))*3600; 
            H2finalBalance = H2_store(end,:) == H2_initial; % final energy inventory must equal the initial inventory
            % battery balance at all hours except the first
            BattBalance(2:end) = Batt_store(2:end) == Batt_store(1:end-1) +...
                (Elec_extract(2:end) - PReq*(27.5/640000))*3600; % in GJ
            %battery balance at the first hour
            BattBalance(1) = Batt_store(1) == Batt_initial +...
                (Elec_extract(1)  - PReq*(27.5/640000))*3600; 
            BattfinalBalance = Batt_store(end,:) == Batt_initial; % final energy inventory must equal the initial inventory
            boundedH2 = H2_store <= H2_buffer; % H2 buffer level must be less than the H2 buffer size
            boundedH2Initial = H2_initial <= H2_buffer;
            boundedBatt = Batt_store <= Batt_buffer; % battery buffer level must be less than the battery size
            boundedBattInitial = Batt_initial <= Batt_buffer;
            % define annualized capital cost
            Capex=annual*((panelkWp_new*panelCapital) ...
                +(turbinekWp_new*turbineCapital) ...
                +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
                +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
                +(Batt_size*PReq*BattpowerCapital/1000) ...
                +HBCapital);
            % define operating cost 
            Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(panel_OnM*panelkWp_new/1000) ...
                +(turbine_OnM*turbinekWp_new*turbineCapital*annual) ...
                +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
                +(HBCapital*HB_OnM);
            combProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
                    
            % add constraints to the problem
            combProblem.Constraints.boundedH2 = boundedH2;
            combProblem.Constraints.boundedH2Initial = boundedH2Initial;
            combProblem.Constraints.boundedBatt = boundedBatt;
            combProblem.Constraints.boundedBattInitial = boundedBattInitial;
            combProblem.Constraints.powerBalance = powerBalance;
            combProblem.Constraints.H2Balance = H2Balance;
            combProblem.Constraints.BattBalance = BattBalance;
            combProblem.Constraints.H2finalBalance = H2finalBalance;
            combProblem.Constraints.BattfinalBalance = BattfinalBalance;
            combProblem.Constraints.Max_H2 = Max_H2;
            combProblem.Constraints.Max_Batt = Max_Batt;
            [sol,fval,exitflag,output] = solve(combProblem); % execute solver
            if exitflag==1 
                panelkWp(j,1)=panelkWp_new;
                turbinekWp(j,1)=turbinekWp_new;
                H2SizeGW(j,1)=sol.H2_size*PReq/1000;%GW
                N2SizeGW(j,1)=N2_size*PReq/1000; %GW
                BattSizeGW(j,1)=sol.Batt_size*PReq/1000; %GW
                H2bufferkg(j,1)=sol.H2_buffer*(1/0.2045); %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=sol.Batt_buffer*277.78;
                LCOA(j,1)=fval;
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000);
            
            else
                panelkWp(j,1)=0;
                turbinekWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1}=0;
                error_index(j,1)=i;
            end
       end
       panel_save(:,y) = panelkWp;
       turbine_save(:,y) = turbinekWp;
       H2Size_save(:,y) = H2SizeGW;
       N2Size_save(:,y) = N2SizeGW;
       BattSize_save(:,y) = BattSizeGW;
       H2buffer_save(:,y) = H2bufferkg;
       N2buffer_save(:,y) = N2bufferkg;
       battery_save(:,y) = batterykWh;
       LCOA_save(:,y) = LCOA;
       batteryUse_save{:,y} = batteryUse;
    end
    CP_costs.topwind_results.panel{i} = panel_save;
    CP_costs.topwind_results.turbine{i} = turbine_save;
    CP_costs.topwind_results.H2Size{i} = H2Size_save;
    CP_costs.topwind_results.N2Size{i} = N2Size_save;
    CP_costs.topwind_results.BattSize{i} = BattSize_save;
    CP_costs.topwind_results.H2buffer{i} = H2buffer_save;
    CP_costs.topwind_results.N2buffer{i} = N2buffer_save;
    CP_costs.topwind_results.batterykWh{i} = battery_save;
    CP_costs.topwind_results.LCOA{i} = LCOA_save;
    CP_costs.topwind_results.batteryUse{i} = batteryUse_save;
end

% 2.iii.b Run the optimization for 0-99% curtailment for top wind locations using only solar energy
CP_costs.topwind_solaronly = struct;
for i = 1:10
    loc = CP_costs.topwind_loc(i,1);
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc;
    panelPowers = CP_base_solar.panelPower(ind);
    years = CP_base_solar.year(ind);
    clear panel_save turbine_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6
        panelPower = panelPowers{y};
        solarkWpGuess = 1e6; % kWp
        x0 = solarkWpGuess;
        inventorySum = @(x) sum((panelPower*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
        x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
        panelkWp_nocurt = x; %store panel kWp 
        N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
        parfor j=1:100
            E_supply = (panelkWp_nocurt/(1-(j-1)/100))*panelPower/1E9; %GW
            
            panelkWp_new = (panelkWp_nocurt/(1-(j-1)/100));
            solarProblem=optimproblem('ObjectiveSense','minimize');
            % Create variables
            H2_size=optimvar('H2_size',1,'LowerBound',1); % size of H2 generation relative to the HB
            Batt_size = optimvar('Batt_size',1); %size of the battery power intake capacity relative to HB
            H2_buffer = optimvar('H2_buffer',1,'LowerBound',0); % H2 buffer storage size in GJ
            Batt_buffer = optimvar('Batt_buffer',1,'LowerBound',0); % Battery buffer storage in GJ
            H2_initial = optimvar('H2_initial',1,'LowerBound',0); % initial level of the H2 buffer in GJ
            Batt_initial = optimvar('Batt_initial',1,'LowerBound',0); % initial level of the battery in GJ
            H2_store = optimvar('H2_store',8760,'LowerBound',0); % level of the H2 buffer at each hour over the year in GJ
            Batt_store = optimvar('Batt_store',8760,'LowerBound',0); % level of the battery at each hour over the year in GJ
            H2_extract = optimvar('H2_extract',8760,'LowerBound',0); % amount of energy fed to H2 production at each hour over the year in GW
            Elec_extract = optimvar('Elec_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)
            % Create constraints
            powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be equal to the energy supply
            Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
            Max_Batt = Elec_extract - (27.5/640)*PReq*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity
            H2Balance = optimconstr(8760); %initialize H2 balance of the year
            BattBalance = optimconstr(8760); % initialize battery balance over the year
            %H2 balance at all hours except the first
            H2Balance(2:end) = H2_store(2:end) == H2_store(1:end-1) +...
                (H2_extract(2:end) - PReq*(612.5/640000))*3600; % in GJ
            %H2 balance at the first hour
            H2Balance(1) = H2_store(1) == H2_initial +...
                (H2_extract(1) - PReq*(612.5/640000))*3600; 
            H2finalBalance = H2_store(end,:) == H2_initial; % final energy inventory must equal the initial inventory
            % battery balance at all hours except the first
            BattBalance(2:end) = Batt_store(2:end) == Batt_store(1:end-1) +...
                (Elec_extract(2:end) - PReq*(27.5/640000))*3600; % in GJ
            %battery balance at the first hour
            BattBalance(1) = Batt_store(1) == Batt_initial +...
                (Elec_extract(1)  - PReq*(27.5/640000))*3600; 
            BattfinalBalance = Batt_store(end,:) == Batt_initial; % final energy inventory must equal the initial inventory
            boundedH2 = H2_store <= H2_buffer; % H2 buffer level must be less than the H2 buffer size
            boundedH2Initial = H2_initial <= H2_buffer;
            boundedBatt = Batt_store <= Batt_buffer; % battery buffer level must be less than the battery size
            boundedBattInitial = Batt_initial <= Batt_buffer; 
            % define annualized capital cost
            Capex=annual*((panelkWp_new*panelCapital) ...
                +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
                +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
                +(Batt_size*PReq*BattpowerCapital/1000) ...
                +HBCapital);
            % define operating cost 
            Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(panel_OnM*panelkWp_new/1000) ...
                +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
                +(HBCapital*HB_OnM);
            solarProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
                    
            % add constraints to the problem
            solarProblem.Constraints.boundedH2 = boundedH2;
            solarProblem.Constraints.boundedH2Initial = boundedH2Initial;
            solarProblem.Constraints.boundedBatt = boundedBatt;
            solarProblem.Constraints.boundedBattInitial = boundedBattInitial;
            solarProblem.Constraints.powerBalance = powerBalance;
            solarProblem.Constraints.H2Balance = H2Balance;
            solarProblem.Constraints.BattBalance = BattBalance;
            solarProblem.Constraints.H2finalBalance = H2finalBalance;
            solarProblem.Constraints.BattfinalBalance = BattfinalBalance;
            solarProblem.Constraints.Max_H2 = Max_H2;
            solarProblem.Constraints.Max_Batt = Max_Batt;
            [sol,fval,exitflag,output] = solve(solarProblem); % execute solver
            if exitflag==1 
                panelkWp(j,1)=panelkWp_new;
                H2SizeGW(j,1)=sol.H2_size*PReq/1000;%GW
                N2SizeGW(j,1)=N2_size*PReq/1000; %GW
                BattSizeGW(j,1)=sol.Batt_size*PReq/1000; %GW
                H2bufferkg(j,1)=sol.H2_buffer*(1/0.2045); %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=sol.Batt_buffer*277.78;
                LCOA(j,1)=fval;
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000);
            
            else
                panelkWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1}=0;
                error_index(j,1)=i;
            end
       end
       panel_save(:,y) = panelkWp;
       H2Size_save(:,y) = H2SizeGW;
       N2Size_save(:,y) = N2SizeGW;
       BattSize_save(:,y) = BattSizeGW;
       H2buffer_save(:,y) = H2bufferkg;
       N2buffer_save(:,y) = N2bufferkg;
       battery_save(:,y) = batterykWh;
       LCOA_save(:,y) = LCOA;
       batteryUse_save{:,y} = batteryUse;
    end
    CP_costs.topwind_solaronly.panel{i} = panel_save;
    CP_costs.topwind_solaronly.H2Size{i} = H2Size_save;
    CP_costs.topwind_solaronly.N2Size{i} = N2Size_save;
    CP_costs.topwind_solaronly.BattSize{i} = BattSize_save;
    CP_costs.topwind_solaronly.H2buffer{i} = H2buffer_save;
    CP_costs.topwind_solaronly.N2buffer{i} = N2buffer_save;
    CP_costs.topwind_solaronly.batterykWh{i} = battery_save;
    CP_costs.topwind_solaronly.LCOA{i} = LCOA_save;
end

% 2.iii.c Run the optimization for 0-99% curtailment for top wind locations
% using only wind energy
CP_costs.topwind_windonly = struct;
for i = 1:10
    loc = CP_costs.topwind_loc(i,1);
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc;
    turbinePowers = CP_base_wind.turbinePower(ind);
    years = CP_base_solar.year(ind);
    clear panel_save turbine_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6
        turbinePower = turbinePowers{y};
        turbinekWpGuess = 1e6; % kWp
        x0 = turbinekWpGuess;
        inventorySum = @(x) sum((turbinePower*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
        x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
        turbinekWp_nocurt = x; %store panel kWp 
        N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
        parfor j=1:100
            E_supply = ((turbinekWp_nocurt/(1-(j-1)/100))*turbinePower/1E9); %GW 
            turbinekWp_new = (turbinekWp_nocurt/(1-(j-1)/100));
            windProblem=optimproblem('ObjectiveSense','minimize');
            % Create variables
            H2_size=optimvar('H2_size',1,'LowerBound',1); % size of H2 generation relative to the HB
            Batt_size = optimvar('Batt_size',1); %size of the battery power intake capacity relative to HB
            H2_buffer = optimvar('H2_buffer',1,'LowerBound',0); % H2 buffer storage size in GJ
            Batt_buffer = optimvar('Batt_buffer',1,'LowerBound',0); % Battery buffer storage in GJ
            H2_initial = optimvar('H2_initial',1,'LowerBound',0); % initial level of the H2 buffer in GJ
            Batt_initial = optimvar('Batt_initial',1,'LowerBound',0); % initial level of the battery in GJ
            H2_store = optimvar('H2_store',8760,'LowerBound',0); % level of the H2 buffer at each hour over the year in GJ
            Batt_store = optimvar('Batt_store',8760,'LowerBound',0); % level of the battery at each hour over the year in GJ
            H2_extract = optimvar('H2_extract',8760,'LowerBound',0); % amount of energy fed to H2 production at each hour over the year in GW
            Elec_extract = optimvar('Elec_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)
            % Create constraints
            powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be equal to the energy supply
            Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
            Max_Batt = Elec_extract - (27.5/640)*PReq*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity
            H2Balance = optimconstr(8760); %initialize H2 balance of the year
            BattBalance = optimconstr(8760); % initialize battery balance over the year
            %H2 balance at all hours except the first
            H2Balance(2:end) = H2_store(2:end) == H2_store(1:end-1) +...
                (H2_extract(2:end) - PReq*(612.5/640000))*3600; % in GJ
            %H2 balance at the first hour
            H2Balance(1) = H2_store(1) == H2_initial +...
                (H2_extract(1) - PReq*(612.5/640000))*3600; 
            H2finalBalance = H2_store(end,:) == H2_initial; % final energy inventory must equal the initial inventory
            % battery balance at all hours except the first
            BattBalance(2:end) = Batt_store(2:end) == Batt_store(1:end-1) +...
                (Elec_extract(2:end) - PReq*(27.5/640000))*3600; % in GJ
            %battery balance at the first hour
            BattBalance(1) = Batt_store(1) == Batt_initial +...
                (Elec_extract(1)  - PReq*(27.5/640000))*3600; 
            BattfinalBalance = Batt_store(end,:) == Batt_initial; % final energy inventory must equal the initial inventory
            boundedH2 = H2_store <= H2_buffer; % H2 buffer level must be less than the H2 buffer size
            boundedH2Initial = H2_initial <= H2_buffer;
            boundedBatt = Batt_store <= Batt_buffer; % battery buffer level must be less than the battery size
            boundedBattInitial = Batt_initial <= Batt_buffer; 
            % define annualized capital cost
            Capex=annual*((turbinekWp_new*turbineCapital) ...
                +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
                +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
                +(Batt_size*PReq*BattpowerCapital/1000) ...
                +HBCapital);
            % define operating cost 
            Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(turbine_OnM*turbinekWp_new*turbineCapital*annual) ...
                +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
                +(HBCapital*HB_OnM);
            windProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
                    
            % add constraints to the problem
            windProblem.Constraints.boundedH2 = boundedH2;
            windProblem.Constraints.boundedH2Initial = boundedH2Initial;
            windProblem.Constraints.boundedBatt = boundedBatt;
            windProblem.Constraints.boundedBattInitial = boundedBattInitial;
            windProblem.Constraints.powerBalance = powerBalance;
            windProblem.Constraints.H2Balance = H2Balance;
            windProblem.Constraints.BattBalance = BattBalance;
            windProblem.Constraints.H2finalBalance = H2finalBalance;
            windProblem.Constraints.BattfinalBalance = BattfinalBalance;
            windProblem.Constraints.Max_H2 = Max_H2;
            windProblem.Constraints.Max_Batt = Max_Batt;
            [sol,fval,exitflag,output] = solve(windProblem); % execute solver
            if exitflag==1 
                turbinekWp(j,1)=turbinekWp_new;
                H2SizeGW(j,1)=sol.H2_size*PReq/1000;%GW
                N2SizeGW(j,1)=N2_size*PReq/1000; %GW
                BattSizeGW(j,1)=sol.Batt_size*PReq/1000; %GW
                H2bufferkg(j,1)=sol.H2_buffer*(1/0.2045); %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=sol.Batt_buffer*277.78;
                LCOA(j,1)=fval;
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000);
            
            else
                turbinekWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1}=0;
                error_index(j,1)=i;
            end
        end
       turbine_save(:,y) = turbinekWp;
       H2Size_save(:,y) = H2SizeGW;
       N2Size_save(:,y) = N2SizeGW;
       BattSize_save(:,y) = BattSizeGW;
       H2buffer_save(:,y) = H2bufferkg;
       N2buffer_save(:,y) = N2bufferkg;
       battery_save(:,y) = batterykWh;
       LCOA_save(:,y) = LCOA;
       batteryUse_save{:,y} = batteryUse;
    end
    CP_costs.topwind_windonly.turbine{i} = turbine_save;
    CP_costs.topwind_windonly.H2Size{i} = H2Size_save;
    CP_costs.topwind_windonly.N2Size{i} = N2Size_save;
    CP_costs.topwind_windonly.BattSize{i} = BattSize_save;
    CP_costs.topwind_windonly.H2buffer{i} = H2buffer_save;
    CP_costs.topwind_windonly.N2buffer{i} = N2buffer_save;
    CP_costs.topwind_windonly.batterykWh{i} = battery_save;
    CP_costs.topwind_windonly.LCOA{i} = LCOA_save;
end
save('CP_costs','CP_costs')

% 2.iv.a Run the optimization for top locations with solar and wind
% combined
CP_costs.top_results = struct; % initialize structure
for i = 1:10
    loc = CP_costs.top_loc(i,1); % extract the location from the list
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc; % get the indeces using the location list that is synced between solar and wind
    panelPowers = CP_base_solar.panelPower(ind); % extract the power profiles of solar energy
    turbinePowers = CP_base_wind.turbinePower(ind); % extract the power profiles of wind energy
    years = CP_base_solar.year(ind); % extract the years
    panelskWp = CP_curt.panelkWp(ind); % extract the optimized panel size with curtailment
    turbineskWp = CP_curt.turbinekWp(ind); % extract the optimized turbine size with curtailment
    clear panel_save turbine_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6 % cycle through the years
        panelPower = panelPowers{y}; % extract the solar power profile for the year
        turbinePower = turbinePowers{y}; % extract the wind power profile for the year
        % combined power profile based on the optimized panel and turbine
        % size
        combPower_perkWp = (panelPower*panelskWp(y) + turbinePower*turbineskWp(y))/(panelskWp(y) + turbineskWp(y));
        % calculate the combined solar and wind kWp size for the case of no
        % curtailment
        combkWpGuess = 1e6; % kWp
        x0 = combkWpGuess;
        inventorySum = @(x) sum((combPower_perkWp*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
        x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
        combkWp_nocurt = x; %store panel kWp 
        % calculate the panel and turbine size individually for the case of
        % no curtailment
        panelkWp_nocurt = combkWp_nocurt * (panelskWp(y)/(turbineskWp(y)+panelskWp(y)));
        turbinekWp_nocurt = combkWp_nocurt * (turbineskWp(y)/(turbineskWp(y)+panelskWp(y)));

        N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
        % run the optimization for 0-99% curtailment, with a constant ratio
        % of turbine size to solar size
        parfor j=1:100
            E_supply = ((turbinekWp_nocurt/(1-(j-1)/100))*turbinePower/1E9) + (panelkWp_nocurt/(1-(j-1)/100))*panelPower/1E9; %GW
            % new panel and turbine size
            panelkWp_new = (panelkWp_nocurt/(1-(j-1)/100));
            turbinekWp_new = (turbinekWp_nocurt/(1-(j-1)/100));
            combProblem=optimproblem('ObjectiveSense','minimize');
            % Create variables
            H2_size=optimvar('H2_size',1,'LowerBound',1); % size of H2 generation relative to the HB
            Batt_size = optimvar('Batt_size',1); %size of the battery power intake capacity relative to HB
            H2_buffer = optimvar('H2_buffer',1,'LowerBound',0); % H2 buffer storage size in GJ
            Batt_buffer = optimvar('Batt_buffer',1,'LowerBound',0); % Battery buffer storage in GJ
            H2_initial = optimvar('H2_initial',1,'LowerBound',0); % initial level of the H2 buffer in GJ
            Batt_initial = optimvar('Batt_initial',1,'LowerBound',0); % initial level of the battery in GJ
            H2_store = optimvar('H2_store',8760,'LowerBound',0); % level of the H2 buffer at each hour over the year in GJ
            Batt_store = optimvar('Batt_store',8760,'LowerBound',0); % level of the battery at each hour over the year in GJ
            H2_extract = optimvar('H2_extract',8760,'LowerBound',0); % amount of energy fed to H2 production at each hour over the year in GW
            Elec_extract = optimvar('Elec_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)
            % Create constraints
            powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be equal to the energy supply
            Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
            Max_Batt = Elec_extract - (27.5/640)*PReq*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity
            H2Balance = optimconstr(8760); %initialize H2 balance of the year
            BattBalance = optimconstr(8760); % initialize battery balance over the year
            %H2 balance at all hours except the first
            H2Balance(2:end) = H2_store(2:end) == H2_store(1:end-1) +...
                (H2_extract(2:end) - PReq*(612.5/640000))*3600; % in GJ
            %H2 balance at the first hour
            H2Balance(1) = H2_store(1) == H2_initial +...
                (H2_extract(1) - PReq*(612.5/640000))*3600; 
            H2finalBalance = H2_store(end,:) == H2_initial; % final energy inventory must equal the initial inventory
            % battery balance at all hours except the first
            BattBalance(2:end) = Batt_store(2:end) == Batt_store(1:end-1) +...
                (Elec_extract(2:end) - PReq*(27.5/640000))*3600; % in GJ
            %battery balance at the first hour
            BattBalance(1) = Batt_store(1) == Batt_initial +...
                (Elec_extract(1)  - PReq*(27.5/640000))*3600; 
            BattfinalBalance = Batt_store(end,:) == Batt_initial; % final energy inventory must equal the initial inventory
            boundedH2 = H2_store <= H2_buffer; % H2 buffer level must be less than the H2 buffer size
            boundedH2Initial = H2_initial <= H2_buffer;
            boundedBatt = Batt_store <= Batt_buffer; % battery buffer level must be less than the battery size
            boundedBattInitial = Batt_initial <= Batt_buffer;
            % define annualized capital cost
            Capex=annual*((panelkWp_new*panelCapital) ...
                +(turbinekWp_new*turbineCapital) ...
                +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
                +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
                +(Batt_size*PReq*BattpowerCapital/1000) ...
                +HBCapital);
            % define operating cost 
            Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(panel_OnM*panelkWp_new/1000) ...
                +(turbine_OnM*turbinekWp_new*turbineCapital*annual) ...
                +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
                +(HBCapital*HB_OnM);
            combProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
                    
            % add constraints to the problem
            combProblem.Constraints.boundedH2 = boundedH2;
            combProblem.Constraints.boundedH2Initial = boundedH2Initial;
            combProblem.Constraints.boundedBatt = boundedBatt;
            combProblem.Constraints.boundedBattInitial = boundedBattInitial;
            combProblem.Constraints.powerBalance = powerBalance;
            combProblem.Constraints.H2Balance = H2Balance;
            combProblem.Constraints.BattBalance = BattBalance;
            combProblem.Constraints.H2finalBalance = H2finalBalance;
            combProblem.Constraints.BattfinalBalance = BattfinalBalance;
            combProblem.Constraints.Max_H2 = Max_H2;
            combProblem.Constraints.Max_Batt = Max_Batt;
            [sol,fval,exitflag,output] = solve(combProblem); % execute solver
            if exitflag==1 
                panelkWp(j,1)=panelkWp_new;
                turbinekWp(j,1)=turbinekWp_new;
                H2SizeGW(j,1)=sol.H2_size*PReq/1000;%GW
                N2SizeGW(j,1)=N2_size*PReq/1000; %GW
                BattSizeGW(j,1)=sol.Batt_size*PReq/1000; %GW
                H2bufferkg(j,1)=sol.H2_buffer*(1/0.2045); %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=sol.Batt_buffer*277.78;
                LCOA(j,1)=fval;
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000);
            
            else
                panelkWp(j,1)=0;
                turbinekWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1}=0;
                error_index(j,1)=i;
            end
       end
       panel_save(:,y) = panelkWp;
       turbine_save(:,y) = turbinekWp;
       H2Size_save(:,y) = H2SizeGW;
       N2Size_save(:,y) = N2SizeGW;
       BattSize_save(:,y) = BattSizeGW;
       H2buffer_save(:,y) = H2bufferkg;
       N2buffer_save(:,y) = N2bufferkg;
       battery_save(:,y) = batterykWh;
       LCOA_save(:,y) = LCOA;
       batteryUse_save{:,y} = batteryUse;
    end
    CP_costs.top_results.panel{i} = panel_save;
    CP_costs.top_results.turbine{i} = turbine_save;
    CP_costs.top_results.H2Size{i} = H2Size_save;
    CP_costs.top_results.N2Size{i} = N2Size_save;
    CP_costs.top_results.BattSize{i} = BattSize_save;
    CP_costs.top_results.H2buffer{i} = H2buffer_save;
    CP_costs.top_results.N2buffer{i} = N2buffer_save;
    CP_costs.top_results.batterykWh{i} = battery_save;
    CP_costs.top_results.LCOA{i} = LCOA_save;
    CP_costs.top_results.batteryUse{i} = batteryUse_save;
end

% 2.iv.b Run the optimization for 0-99% curtailment for top locations using only solar energy
CP_costs.top_solaronly = struct;
for i = 1:10
    loc = CP_costs.top_loc(i,1);
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc;
    panelPowers = CP_base_solar.panelPower(ind);
    years = CP_base_solar.year(ind);
    clear panel_save turbine_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6
        panelPower = panelPowers{y};
        solarkWpGuess = 1e6; % kWp
        x0 = solarkWpGuess;
        inventorySum = @(x) sum((panelPower*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
        x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
        panelkWp_nocurt = x; %store panel kWp 
        N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
        parfor j=1:100
            E_supply = (panelkWp_nocurt/(1-(j-1)/100))*panelPower/1E9; %GW
            
            panelkWp_new = (panelkWp_nocurt/(1-(j-1)/100));
            solarProblem=optimproblem('ObjectiveSense','minimize');
            % Create variables
            H2_size=optimvar('H2_size',1,'LowerBound',1); % size of H2 generation relative to the HB
            Batt_size = optimvar('Batt_size',1); %size of the battery power intake capacity relative to HB
            H2_buffer = optimvar('H2_buffer',1,'LowerBound',0); % H2 buffer storage size in GJ
            Batt_buffer = optimvar('Batt_buffer',1,'LowerBound',0); % Battery buffer storage in GJ
            H2_initial = optimvar('H2_initial',1,'LowerBound',0); % initial level of the H2 buffer in GJ
            Batt_initial = optimvar('Batt_initial',1,'LowerBound',0); % initial level of the battery in GJ
            H2_store = optimvar('H2_store',8760,'LowerBound',0); % level of the H2 buffer at each hour over the year in GJ
            Batt_store = optimvar('Batt_store',8760,'LowerBound',0); % level of the battery at each hour over the year in GJ
            H2_extract = optimvar('H2_extract',8760,'LowerBound',0); % amount of energy fed to H2 production at each hour over the year in GW
            Elec_extract = optimvar('Elec_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)
            % Create constraints
            powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be equal to the energy supply
            Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
            Max_Batt = Elec_extract - (27.5/640)*PReq*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity
            H2Balance = optimconstr(8760); %initialize H2 balance of the year
            BattBalance = optimconstr(8760); % initialize battery balance over the year
            %H2 balance at all hours except the first
            H2Balance(2:end) = H2_store(2:end) == H2_store(1:end-1) +...
                (H2_extract(2:end) - PReq*(612.5/640000))*3600; % in GJ
            %H2 balance at the first hour
            H2Balance(1) = H2_store(1) == H2_initial +...
                (H2_extract(1) - PReq*(612.5/640000))*3600; 
            H2finalBalance = H2_store(end,:) == H2_initial; % final energy inventory must equal the initial inventory
            % battery balance at all hours except the first
            BattBalance(2:end) = Batt_store(2:end) == Batt_store(1:end-1) +...
                (Elec_extract(2:end) - PReq*(27.5/640000))*3600; % in GJ
            %battery balance at the first hour
            BattBalance(1) = Batt_store(1) == Batt_initial +...
                (Elec_extract(1)  - PReq*(27.5/640000))*3600; 
            BattfinalBalance = Batt_store(end,:) == Batt_initial; % final energy inventory must equal the initial inventory
            boundedH2 = H2_store <= H2_buffer; % H2 buffer level must be less than the H2 buffer size
            boundedH2Initial = H2_initial <= H2_buffer;
            boundedBatt = Batt_store <= Batt_buffer; % battery buffer level must be less than the battery size
            boundedBattInitial = Batt_initial <= Batt_buffer; 
            % define annualized capital cost
            Capex=annual*((panelkWp_new*panelCapital) ...
                +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
                +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
                +(Batt_size*PReq*BattpowerCapital/1000) ...
                +HBCapital);
            % define operating cost 
            Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(panel_OnM*panelkWp_new/1000) ...
                +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
                +(HBCapital*HB_OnM);
            solarProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
                    
            % add constraints to the problem
            solarProblem.Constraints.boundedH2 = boundedH2;
            solarProblem.Constraints.boundedH2Initial = boundedH2Initial;
            solarProblem.Constraints.boundedBatt = boundedBatt;
            solarProblem.Constraints.boundedBattInitial = boundedBattInitial;
            solarProblem.Constraints.powerBalance = powerBalance;
            solarProblem.Constraints.H2Balance = H2Balance;
            solarProblem.Constraints.BattBalance = BattBalance;
            solarProblem.Constraints.H2finalBalance = H2finalBalance;
            solarProblem.Constraints.BattfinalBalance = BattfinalBalance;
            solarProblem.Constraints.Max_H2 = Max_H2;
            solarProblem.Constraints.Max_Batt = Max_Batt;
            [sol,fval,exitflag,output] = solve(solarProblem); % execute solver
            if exitflag==1 
                panelkWp(j,1)=panelkWp_new;
                H2SizeGW(j,1)=sol.H2_size*PReq/1000;%GW
                N2SizeGW(j,1)=N2_size*PReq/1000; %GW
                BattSizeGW(j,1)=sol.Batt_size*PReq/1000; %GW
                H2bufferkg(j,1)=sol.H2_buffer*(1/0.2045); %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=sol.Batt_buffer*277.78;
                LCOA(j,1)=fval;
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000);
            
            else
                panelkWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1}=0;
                error_index(j,1)=i;
            end
       end
       panel_save(:,y) = panelkWp;
       H2Size_save(:,y) = H2SizeGW;
       N2Size_save(:,y) = N2SizeGW;
       BattSize_save(:,y) = BattSizeGW;
       H2buffer_save(:,y) = H2bufferkg;
       N2buffer_save(:,y) = N2bufferkg;
       battery_save(:,y) = batterykWh;
       LCOA_save(:,y) = LCOA;
       batteryUse_save{:,y} = batteryUse;
    end
    CP_costs.top_solaronly.panel{i} = panel_save;
    CP_costs.top_solaronly.H2Size{i} = H2Size_save;
    CP_costs.top_solaronly.N2Size{i} = N2Size_save;
    CP_costs.top_solaronly.BattSize{i} = BattSize_save;
    CP_costs.top_solaronly.H2buffer{i} = H2buffer_save;
    CP_costs.top_solaronly.N2buffer{i} = N2buffer_save;
    CP_costs.top_solaronly.batterykWh{i} = battery_save;
    CP_costs.top_solaronly.LCOA{i} = LCOA_save;
end

% 2.iv.c Run the optimization for 0-99% curtailment for top locations
% using only wind energy
CP_costs.top_windonly = struct;
for i = 1:10
    loc = CP_costs.top_loc(i,1);
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc;
    turbinePowers = CP_base_wind.turbinePower(ind);
    years = CP_base_solar.year(ind);
    clear panel_save turbine_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6
        turbinePower = turbinePowers{y};
        turbinekWpGuess = 1e6; % kWp
        x0 = turbinekWpGuess;
        inventorySum = @(x) sum((turbinePower*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
        x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
        turbinekWp_nocurt = x; %store panel kWp 
        N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
        parfor j=1:100
            E_supply = ((turbinekWp_nocurt/(1-(j-1)/100))*turbinePower/1E9); %GW 
            turbinekWp_new = (turbinekWp_nocurt/(1-(j-1)/100));
            windProblem=optimproblem('ObjectiveSense','minimize');
            % Create variables
            H2_size=optimvar('H2_size',1,'LowerBound',1); % size of H2 generation relative to the HB
            Batt_size = optimvar('Batt_size',1); %size of the battery power intake capacity relative to HB
            H2_buffer = optimvar('H2_buffer',1,'LowerBound',0); % H2 buffer storage size in GJ
            Batt_buffer = optimvar('Batt_buffer',1,'LowerBound',0); % Battery buffer storage in GJ
            H2_initial = optimvar('H2_initial',1,'LowerBound',0); % initial level of the H2 buffer in GJ
            Batt_initial = optimvar('Batt_initial',1,'LowerBound',0); % initial level of the battery in GJ
            H2_store = optimvar('H2_store',8760,'LowerBound',0); % level of the H2 buffer at each hour over the year in GJ
            Batt_store = optimvar('Batt_store',8760,'LowerBound',0); % level of the battery at each hour over the year in GJ
            H2_extract = optimvar('H2_extract',8760,'LowerBound',0); % amount of energy fed to H2 production at each hour over the year in GW
            Elec_extract = optimvar('Elec_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)
            % Create constraints
            powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be equal to the energy supply
            Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
            Max_Batt = Elec_extract - (27.5/640)*PReq*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity
            H2Balance = optimconstr(8760); %initialize H2 balance of the year
            BattBalance = optimconstr(8760); % initialize battery balance over the year
            %H2 balance at all hours except the first
            H2Balance(2:end) = H2_store(2:end) == H2_store(1:end-1) +...
                (H2_extract(2:end) - PReq*(612.5/640000))*3600; % in GJ
            %H2 balance at the first hour
            H2Balance(1) = H2_store(1) == H2_initial +...
                (H2_extract(1) - PReq*(612.5/640000))*3600; 
            H2finalBalance = H2_store(end,:) == H2_initial; % final energy inventory must equal the initial inventory
            % battery balance at all hours except the first
            BattBalance(2:end) = Batt_store(2:end) == Batt_store(1:end-1) +...
                (Elec_extract(2:end) - PReq*(27.5/640000))*3600; % in GJ
            %battery balance at the first hour
            BattBalance(1) = Batt_store(1) == Batt_initial +...
                (Elec_extract(1)  - PReq*(27.5/640000))*3600; 
            BattfinalBalance = Batt_store(end,:) == Batt_initial; % final energy inventory must equal the initial inventory
            boundedH2 = H2_store <= H2_buffer; % H2 buffer level must be less than the H2 buffer size
            boundedH2Initial = H2_initial <= H2_buffer;
            boundedBatt = Batt_store <= Batt_buffer; % battery buffer level must be less than the battery size
            boundedBattInitial = Batt_initial <= Batt_buffer; 
            % define annualized capital cost
            Capex=annual*((turbinekWp_new*turbineCapital) ...
                +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
                +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
                +(Batt_size*PReq*BattpowerCapital/1000) ...
                +HBCapital);
            % define operating cost 
            Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
                +(turbine_OnM*turbinekWp_new*turbineCapital*annual) ...
                +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
                +(HBCapital*HB_OnM);
            windProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
                    
            % add constraints to the problem
            windProblem.Constraints.boundedH2 = boundedH2;
            windProblem.Constraints.boundedH2Initial = boundedH2Initial;
            windProblem.Constraints.boundedBatt = boundedBatt;
            windProblem.Constraints.boundedBattInitial = boundedBattInitial;
            windProblem.Constraints.powerBalance = powerBalance;
            windProblem.Constraints.H2Balance = H2Balance;
            windProblem.Constraints.BattBalance = BattBalance;
            windProblem.Constraints.H2finalBalance = H2finalBalance;
            windProblem.Constraints.BattfinalBalance = BattfinalBalance;
            windProblem.Constraints.Max_H2 = Max_H2;
            windProblem.Constraints.Max_Batt = Max_Batt;
            [sol,fval,exitflag,output] = solve(windProblem); % execute solver
            if exitflag==1 
                turbinekWp(j,1)=turbinekWp_new;
                H2SizeGW(j,1)=sol.H2_size*PReq/1000;%GW
                N2SizeGW(j,1)=N2_size*PReq/1000; %GW
                BattSizeGW(j,1)=sol.Batt_size*PReq/1000; %GW
                H2bufferkg(j,1)=sol.H2_buffer*(1/0.2045); %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=sol.Batt_buffer*277.78;
                LCOA(j,1)=fval;
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000);
            
            else
                turbinekWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1}=0;
                error_index(j,1)=i;
            end
        end
       turbine_save(:,y) = turbinekWp;
       H2Size_save(:,y) = H2SizeGW;
       N2Size_save(:,y) = N2SizeGW;
       BattSize_save(:,y) = BattSizeGW;
       H2buffer_save(:,y) = H2bufferkg;
       N2buffer_save(:,y) = N2bufferkg;
       battery_save(:,y) = batterykWh;
       LCOA_save(:,y) = LCOA;
       batteryUse_save{:,y} = batteryUse;
    end
    CP_costs.top_windonly.turbine{i} = turbine_save;
    CP_costs.top_windonly.H2Size{i} = H2Size_save;
    CP_costs.top_windonly.N2Size{i} = N2Size_save;
    CP_costs.top_windonly.BattSize{i} = BattSize_save;
    CP_costs.top_windonly.H2buffer{i} = H2buffer_save;
    CP_costs.top_windonly.N2buffer{i} = N2buffer_save;
    CP_costs.top_windonly.batterykWh{i} = battery_save;
    CP_costs.top_windonly.LCOA{i} = LCOA_save;
end
save('CP_costs','CP_costs')


% 3.i.a Calculate the utilization costs for top locations combined solar and
% wind
clearvars -except CP_costs panelCapital panel_OnM turbineCapital turbine_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
results = CP_costs.top_results;
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_turbineCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    panels = flip(results.panel{i},1);
    turbines = flip(results.turbine{i},1);
    H2Sizes = flip(results.H2Size{i},1);
    N2Sizes = flip(results.N2Size{i},1);
    BattSizes = flip(results.BattSize{i},1);
    H2buffers = flip(results.H2buffer{i},1);
    N2buffers = flip(results.N2buffer{i},1);
    batteries = flip(results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel and turbine size for full energy
        %utilization
        E_cost = (results.panel{i}(1,y) * (annual * panelCapital + (1/1000)*panel_OnM) ...
            +(results.turbine{i}(1,y)) * (1 + turbine_OnM) * annual * turbineCapital ) /(876000);
        % find the cost of each process component for each percent energy
        % utilization
        panelCost = panels(:,y) * (annual * panelCapital + (1/1000)*panel_OnM);
        turbineCost = turbines(:,y) * (1 + turbine_OnM) * annual * turbineCapital;
        H2SizeCost = H2Sizes(:,y) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM);
        BattSizeCost = BattSizes(:,y) * BattpowerCapital * annual;
        H2bufferCost = H2buffers(:,y) * H2storageCapital * annual;
        BatteryCost = batteries(:,y) * BattstorageCapital * annual;
        % find the change in the cost of each process component after
        % rescaling by the percent energy utilization such that the solar
        % panel and wind turbine size is constant rather than the total production - as was
        % the case during the process optimization for each percent
        % curtailment
        delta_panelCost(:,y) = (panelCost(2:end,1).*(S(2:end,1)/100) - panelCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_turbineCost(:,y) = (turbineCost(2:end,1).*(S(2:end,1)/100) - turbineCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2SizeCost(:,y) = (H2SizeCost(2:end,1).*(S(2:end,1)/100) - H2SizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BattSizeCost(:,y) = (BattSizeCost(2:end,1).*(S(2:end,1)/100) - BattSizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2bufferCost(:,y) = (H2bufferCost(2:end,1).*(S(2:end,1)/100) - H2bufferCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BatteryCost(:,y) = (BatteryCost(2:end,1).*(S(2:end,1)/100) - BatteryCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        %save the energy cost for the location and year in the CP_value
        %structure
        CP_costs.top_Ecost(i,y) = E_cost;
    end
    % add the change in cost for the first step as well
    delta_panelCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_panelCost];
    delta_turbineCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_turbineCost];
    delta_H2SizeCost = [H2Sizes(1,:) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/100) * (1/8760) ; delta_H2SizeCost];
    delta_BattSizeCost = [BattSizes(1,:) * BattpowerCapital * annual * (1/100) * (1/8760) ; delta_BattSizeCost];
    delta_H2bufferCost = [H2buffers(1,:) * H2storageCapital * annual * (1/100) * (1/8760) ; delta_H2bufferCost];
    delta_BatteryCost = [batteries(1,:) * BattstorageCapital * annual * (1/100) * (1/8760) ; delta_BatteryCost];
    % find the change in cost for the steady-state energy utilization
    % (constant array)
    delta_H2SizeCost_SS = ones(100,6) * (612.5/640) * (1/1000) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/8760);
    delta_ASUSizeCost_SS = ones(100,6) * (22.5/640) * (1/1000) * (ASUCapital+N2compCapital) * (annual + HB_OnM) * (1/8760);
    delta_HBSizeCost_SS = ones(100,6) * (1/100) * (HBCapital) * (annual + HB_OnM) * (1/8760);
    %create table of the added cost for each of the process components as a
    %function of energy utilization. add to the CP_value structure
    costs = table(mean(delta_panelCost,2), ...
        mean(delta_turbineCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["panel_cost","turbine_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    CP_costs.top_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    CP_costs.top_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    CP_costs.top_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end
% 3.i.b Calculate the utilization costs for top locations with only solar
clearvars -except CP_costs panelCapital panel_OnM turbineCapital turbine_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
results = CP_costs.top_solaronly;
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_turbineCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    panels = flip(results.panel{i},1);
    H2Sizes = flip(results.H2Size{i},1);
    N2Sizes = flip(results.N2Size{i},1);
    BattSizes = flip(results.BattSize{i},1);
    H2buffers = flip(results.H2buffer{i},1);
    N2buffers = flip(results.N2buffer{i},1);
    batteries = flip(results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel and turbine size for full energy
        %utilization
        E_cost = (results.panel{i}(1,y) * (annual * panelCapital + (1/1000)*panel_OnM)) /(876000);
        % find the cost of each process component for each percent energy
        % utilization
        panelCost = panels(:,y) * (annual * panelCapital + (1/1000)*panel_OnM);
        H2SizeCost = H2Sizes(:,y) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM);
        BattSizeCost = BattSizes(:,y) * BattpowerCapital * annual;
        H2bufferCost = H2buffers(:,y) * H2storageCapital * annual;
        BatteryCost = batteries(:,y) * BattstorageCapital * annual;
        % find the change in the cost of each process component after
        % rescaling by the percent energy utilization such that the solar
        % panel size is constant rather than the total production - as was
        % the case during the process optimization for each percent
        % curtailment
        delta_panelCost(:,y) = (panelCost(2:end,1).*(S(2:end,1)/100) - panelCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2SizeCost(:,y) = (H2SizeCost(2:end,1).*(S(2:end,1)/100) - H2SizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BattSizeCost(:,y) = (BattSizeCost(2:end,1).*(S(2:end,1)/100) - BattSizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2bufferCost(:,y) = (H2bufferCost(2:end,1).*(S(2:end,1)/100) - H2bufferCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BatteryCost(:,y) = (BatteryCost(2:end,1).*(S(2:end,1)/100) - BatteryCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        %save the energy cost for the location and year in the CP_value
        %structure
        CP_costs.top_solaronly_Ecost(i,y) = E_cost;
    end
    % add the change in cost for the first step as well
    delta_panelCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_panelCost];
    delta_H2SizeCost = [H2Sizes(1,:) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/100) * (1/8760) ; delta_H2SizeCost];
    delta_BattSizeCost = [BattSizes(1,:) * BattpowerCapital * annual * (1/100) * (1/8760) ; delta_BattSizeCost];
    delta_H2bufferCost = [H2buffers(1,:) * H2storageCapital * annual * (1/100) * (1/8760) ; delta_H2bufferCost];
    delta_BatteryCost = [batteries(1,:) * BattstorageCapital * annual * (1/100) * (1/8760) ; delta_BatteryCost];
    % find the change in cost for the steady-state energy utilization
    % (constant array)
    delta_H2SizeCost_SS = ones(100,6) * (612.5/640) * (1/1000) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/8760);
    delta_ASUSizeCost_SS = ones(100,6) * (22.5/640) * (1/1000) * (ASUCapital+N2compCapital) * (annual + HB_OnM) * (1/8760);
    delta_HBSizeCost_SS = ones(100,6) * (1/100) * (HBCapital) * (annual + HB_OnM) * (1/8760);
    %create table of the added cost for each of the process components as a
    %function of energy utilization. add to the CP_value structure
    costs = table(mean(delta_panelCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["panel_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    CP_costs.top_solaronly_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    CP_costs.top_solaronly_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    CP_costs.top_solaronly_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end
% 3.i.c Calculate the cost and value for top locations with only wind
clearvars -except CP_costs panelCapital panel_OnM turbineCapital turbine_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
results = CP_costs.top_windonly;
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_turbineCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    turbines = flip(results.turbine{i},1);
    H2Sizes = flip(results.H2Size{i},1);
    N2Sizes = flip(results.N2Size{i},1);
    BattSizes = flip(results.BattSize{i},1);
    H2buffers = flip(results.H2buffer{i},1);
    N2buffers = flip(results.N2buffer{i},1);
    batteries = flip(results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel and turbine size for full energy
        %utilization
        E_cost = ((results.turbine{i}(1,y)) * (1 + turbine_OnM) * annual * turbineCapital ) /(876000);
        % find the cost of each process component for each percent energy
        % utilization
        turbineCost = turbines(:,y) * (1 + turbine_OnM) * annual * turbineCapital;
        H2SizeCost = H2Sizes(:,y) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM);
        BattSizeCost = BattSizes(:,y) * BattpowerCapital * annual;
        H2bufferCost = H2buffers(:,y) * H2storageCapital * annual;
        BatteryCost = batteries(:,y) * BattstorageCapital * annual;
        % find the change in the cost of each process component after
        % rescaling by the percent energy utilization such that the wind
        % turbine size is constant rather than the total production - as was
        % the case during the process optimization for each percent
        % curtailment
        delta_turbineCost(:,y) = (turbineCost(2:end,1).*(S(2:end,1)/100) - turbineCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2SizeCost(:,y) = (H2SizeCost(2:end,1).*(S(2:end,1)/100) - H2SizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BattSizeCost(:,y) = (BattSizeCost(2:end,1).*(S(2:end,1)/100) - BattSizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2bufferCost(:,y) = (H2bufferCost(2:end,1).*(S(2:end,1)/100) - H2bufferCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BatteryCost(:,y) = (BatteryCost(2:end,1).*(S(2:end,1)/100) - BatteryCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        %save the energy cost for the location and year in the CP_value
        %structure
        CP_costs.top_windonly_Ecost(i,y) = E_cost;
    end
    % add the change in cost for the first step as well
    delta_turbineCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_turbineCost];
    delta_H2SizeCost = [H2Sizes(1,:) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/100) * (1/8760) ; delta_H2SizeCost];
    delta_BattSizeCost = [BattSizes(1,:) * BattpowerCapital * annual * (1/100) * (1/8760) ; delta_BattSizeCost];
    delta_H2bufferCost = [H2buffers(1,:) * H2storageCapital * annual * (1/100) * (1/8760) ; delta_H2bufferCost];
    delta_BatteryCost = [batteries(1,:) * BattstorageCapital * annual * (1/100) * (1/8760) ; delta_BatteryCost];
    % find the change in cost for the steady-state energy utilization
    % (constant array)
    delta_H2SizeCost_SS = ones(100,6) * (612.5/640) * (1/1000) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/8760);
    delta_ASUSizeCost_SS = ones(100,6) * (22.5/640) * (1/1000) * (ASUCapital+N2compCapital) * (annual + HB_OnM) * (1/8760);
    delta_HBSizeCost_SS = ones(100,6) * (1/100) * (HBCapital) * (annual + HB_OnM) * (1/8760);
    %create table of the added cost for each of the process components as a
    %function of energy utilization. add to the CP_value structure
    costs = table(mean(delta_turbineCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["turbine_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    CP_costs.top_windonly_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    CP_costs.top_windonly_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    CP_costs.top_windonly_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end


% 3.ii.a Calculate the utilziation costs for top solar locations combined solar and
% wind
clearvars -except CP_costs panelCapital panel_OnM turbineCapital turbine_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
results = CP_costs.topsolar_results;
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_turbineCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    panels = flip(results.panel{i},1);
    turbines = flip(results.turbine{i},1);
    H2Sizes = flip(results.H2Size{i},1);
    N2Sizes = flip(results.N2Size{i},1);
    BattSizes = flip(results.BattSize{i},1);
    H2buffers = flip(results.H2buffer{i},1);
    N2buffers = flip(results.N2buffer{i},1);
    batteries = flip(results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel and turbine size for full energy
        %utilization
        E_cost = (results.panel{i}(1,y) * (annual * panelCapital + (1/1000)*panel_OnM) ...
            +(results.turbine{i}(1,y)) * (1 + turbine_OnM) * annual * turbineCapital ) /(876000);
        % find the cost of each process component for each percent energy
        % utilization
        panelCost = panels(:,y) * (annual * panelCapital + (1/1000)*panel_OnM);
        turbineCost = turbines(:,y) * (1 + turbine_OnM) * annual * turbineCapital;
        H2SizeCost = H2Sizes(:,y) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM);
        BattSizeCost = BattSizes(:,y) * BattpowerCapital * annual;
        H2bufferCost = H2buffers(:,y) * H2storageCapital * annual;
        BatteryCost = batteries(:,y) * BattstorageCapital * annual;
        % find the change in the cost of each process component after
        % rescaling by the percent energy utilization such that the solar
        % panel and wind turbine size is constant rather than the total production - as was
        % the case during the process optimization for each percent
        % curtailment
        delta_panelCost(:,y) = (panelCost(2:end,1).*(S(2:end,1)/100) - panelCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_turbineCost(:,y) = (turbineCost(2:end,1).*(S(2:end,1)/100) - turbineCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2SizeCost(:,y) = (H2SizeCost(2:end,1).*(S(2:end,1)/100) - H2SizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BattSizeCost(:,y) = (BattSizeCost(2:end,1).*(S(2:end,1)/100) - BattSizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2bufferCost(:,y) = (H2bufferCost(2:end,1).*(S(2:end,1)/100) - H2bufferCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BatteryCost(:,y) = (BatteryCost(2:end,1).*(S(2:end,1)/100) - BatteryCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        %save the energy cost for the location and year in the CP_value
        %structure
        CP_costs.topsolar_Ecost(i,y) = E_cost;
    end
    % add the change in cost for the first step as well
    delta_panelCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_panelCost];
    delta_turbineCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_turbineCost];
    delta_H2SizeCost = [H2Sizes(1,:) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/100) * (1/8760) ; delta_H2SizeCost];
    delta_BattSizeCost = [BattSizes(1,:) * BattpowerCapital * annual * (1/100) * (1/8760) ; delta_BattSizeCost];
    delta_H2bufferCost = [H2buffers(1,:) * H2storageCapital * annual * (1/100) * (1/8760) ; delta_H2bufferCost];
    delta_BatteryCost = [batteries(1,:) * BattstorageCapital * annual * (1/100) * (1/8760) ; delta_BatteryCost];
    % find the change in cost for the steady-state energy utilization
    % (constant array)
    delta_H2SizeCost_SS = ones(100,6) * (612.5/640) * (1/1000) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/8760);
    delta_ASUSizeCost_SS = ones(100,6) * (22.5/640) * (1/1000) * (ASUCapital+N2compCapital) * (annual + HB_OnM) * (1/8760);
    delta_HBSizeCost_SS = ones(100,6) * (1/100) * (HBCapital) * (annual + HB_OnM) * (1/8760);
    %create table of the added cost for each of the process components as a
    %function of energy utilization. add to the CP_value structure
    costs = table(mean(delta_panelCost,2), ...
        mean(delta_turbineCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["panel_cost","turbine_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    CP_costs.topsolar_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    CP_costs.topsolar_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    CP_costs.topsolar_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end
% 3.ii.b Calculate the utilziation costs for top solar locations with only solar
clearvars -except CP_costs panelCapital panel_OnM turbineCapital turbine_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
results = CP_costs.topsolar_solaronly;
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_turbineCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    panels = flip(results.panel{i},1);
    H2Sizes = flip(results.H2Size{i},1);
    N2Sizes = flip(results.N2Size{i},1);
    BattSizes = flip(results.BattSize{i},1);
    H2buffers = flip(results.H2buffer{i},1);
    N2buffers = flip(results.N2buffer{i},1);
    batteries = flip(results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel and turbine size for full energy
        %utilization
        E_cost = (results.panel{i}(1,y) * (annual * panelCapital + (1/1000)*panel_OnM)) /(876000);
        % find the cost of each process component for each percent energy
        % utilization
        panelCost = panels(:,y) * (annual * panelCapital + (1/1000)*panel_OnM);
        H2SizeCost = H2Sizes(:,y) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM);
        BattSizeCost = BattSizes(:,y) * BattpowerCapital * annual;
        H2bufferCost = H2buffers(:,y) * H2storageCapital * annual;
        BatteryCost = batteries(:,y) * BattstorageCapital * annual;
        % find the change in the cost of each process component after
        % rescaling by the percent energy utilization such that the solar
        % panel size is constant rather than the total production - as was
        % the case during the process optimization for each percent
        % curtailment
        delta_panelCost(:,y) = (panelCost(2:end,1).*(S(2:end,1)/100) - panelCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2SizeCost(:,y) = (H2SizeCost(2:end,1).*(S(2:end,1)/100) - H2SizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BattSizeCost(:,y) = (BattSizeCost(2:end,1).*(S(2:end,1)/100) - BattSizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2bufferCost(:,y) = (H2bufferCost(2:end,1).*(S(2:end,1)/100) - H2bufferCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BatteryCost(:,y) = (BatteryCost(2:end,1).*(S(2:end,1)/100) - BatteryCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        %save the energy cost for the location and year in the CP_value
        %structure
        CP_costs.topsolar_solaronly_Ecost(i,y) = E_cost;
    end
    % add the change in cost for the first step as well
    delta_panelCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_panelCost];
    delta_H2SizeCost = [H2Sizes(1,:) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/100) * (1/8760) ; delta_H2SizeCost];
    delta_BattSizeCost = [BattSizes(1,:) * BattpowerCapital * annual * (1/100) * (1/8760) ; delta_BattSizeCost];
    delta_H2bufferCost = [H2buffers(1,:) * H2storageCapital * annual * (1/100) * (1/8760) ; delta_H2bufferCost];
    delta_BatteryCost = [batteries(1,:) * BattstorageCapital * annual * (1/100) * (1/8760) ; delta_BatteryCost];
    % find the change in cost for the steady-state energy utilization
    % (constant array)
    delta_H2SizeCost_SS = ones(100,6) * (612.5/640) * (1/1000) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/8760);
    delta_ASUSizeCost_SS = ones(100,6) * (22.5/640) * (1/1000) * (ASUCapital+N2compCapital) * (annual + HB_OnM) * (1/8760);
    delta_HBSizeCost_SS = ones(100,6) * (1/100) * (HBCapital) * (annual + HB_OnM) * (1/8760);
    %create table of the added cost for each of the process components as a
    %function of energy utilization. add to the CP_value structure
    costs = table(mean(delta_panelCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["panel_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    CP_costs.topsolar_solaronly_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    CP_costs.topsolar_solaronly_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    CP_costs.topsolar_solaronly_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end
% 3.ii.c Calculate the utilization costs for top solar locations with only wind
clearvars -except CP_costs panelCapital panel_OnM turbineCapital turbine_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
results = CP_costs.topsolar_windonly;
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_turbineCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    turbines = flip(results.turbine{i},1);
    H2Sizes = flip(results.H2Size{i},1);
    N2Sizes = flip(results.N2Size{i},1);
    BattSizes = flip(results.BattSize{i},1);
    H2buffers = flip(results.H2buffer{i},1);
    N2buffers = flip(results.N2buffer{i},1);
    batteries = flip(results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel and turbine size for full energy
        %utilization
        E_cost = ((results.turbine{i}(1,y)) * (1 + turbine_OnM) * annual * turbineCapital ) /(876000);
        % find the cost of each process component for each percent energy
        % utilization
        turbineCost = turbines(:,y) * (1 + turbine_OnM) * annual * turbineCapital;
        H2SizeCost = H2Sizes(:,y) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM);
        BattSizeCost = BattSizes(:,y) * BattpowerCapital * annual;
        H2bufferCost = H2buffers(:,y) * H2storageCapital * annual;
        BatteryCost = batteries(:,y) * BattstorageCapital * annual;
        % find the change in the cost of each process component after
        % rescaling by the percent energy utilization such that the wind
        % turbine size is constant rather than the total production - as was
        % the case during the process optimization for each percent
        % curtailment
        delta_turbineCost(:,y) = (turbineCost(2:end,1).*(S(2:end,1)/100) - turbineCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2SizeCost(:,y) = (H2SizeCost(2:end,1).*(S(2:end,1)/100) - H2SizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BattSizeCost(:,y) = (BattSizeCost(2:end,1).*(S(2:end,1)/100) - BattSizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2bufferCost(:,y) = (H2bufferCost(2:end,1).*(S(2:end,1)/100) - H2bufferCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BatteryCost(:,y) = (BatteryCost(2:end,1).*(S(2:end,1)/100) - BatteryCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        %save the energy cost for the location and year in the CP_value
        %structure
        CP_costs.topsolar_windonly_Ecost(i,y) = E_cost;
    end
    % add the change in cost for the first step as well
    delta_turbineCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_turbineCost];
    delta_H2SizeCost = [H2Sizes(1,:) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/100) * (1/8760) ; delta_H2SizeCost];
    delta_BattSizeCost = [BattSizes(1,:) * BattpowerCapital * annual * (1/100) * (1/8760) ; delta_BattSizeCost];
    delta_H2bufferCost = [H2buffers(1,:) * H2storageCapital * annual * (1/100) * (1/8760) ; delta_H2bufferCost];
    delta_BatteryCost = [batteries(1,:) * BattstorageCapital * annual * (1/100) * (1/8760) ; delta_BatteryCost];
    % find the change in cost for the steady-state energy utilization
    % (constant array)
    delta_H2SizeCost_SS = ones(100,6) * (612.5/640) * (1/1000) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/8760);
    delta_ASUSizeCost_SS = ones(100,6) * (22.5/640) * (1/1000) * (ASUCapital+N2compCapital) * (annual + HB_OnM) * (1/8760);
    delta_HBSizeCost_SS = ones(100,6) * (1/100) * (HBCapital) * (annual + HB_OnM) * (1/8760);
    %create table of the added cost for each of the process components as a
    %function of energy utilization. add to the CP_value structure
    costs = table(mean(delta_turbineCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["turbine_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    CP_costs.topsolar_windonly_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    CP_costs.topsolar_windonly_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    CP_costs.topsolar_windonly_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end

% 3.iii.a Calculate the utilization costs for top wind locations combined solar and
% wind
clearvars -except CP_costs panelCapital panel_OnM turbineCapital turbine_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
results = CP_costs.topwind_results;
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_turbineCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    panels = flip(results.panel{i},1);
    turbines = flip(results.turbine{i},1);
    H2Sizes = flip(results.H2Size{i},1);
    N2Sizes = flip(results.N2Size{i},1);
    BattSizes = flip(results.BattSize{i},1);
    H2buffers = flip(results.H2buffer{i},1);
    N2buffers = flip(results.N2buffer{i},1);
    batteries = flip(results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel and turbine size for full energy
        %utilization
        E_cost = (results.panel{i}(1,y) * (annual * panelCapital + (1/1000)*panel_OnM) ...
            +(results.turbine{i}(1,y)) * (1 + turbine_OnM) * annual * turbineCapital ) /(876000);
        % find the cost of each process component for each percent energy
        % utilization
        panelCost = panels(:,y) * (annual * panelCapital + (1/1000)*panel_OnM);
        turbineCost = turbines(:,y) * (1 + turbine_OnM) * annual * turbineCapital;
        H2SizeCost = H2Sizes(:,y) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM);
        BattSizeCost = BattSizes(:,y) * BattpowerCapital * annual;
        H2bufferCost = H2buffers(:,y) * H2storageCapital * annual;
        BatteryCost = batteries(:,y) * BattstorageCapital * annual;
        % find the change in the cost of each process component after
        % rescaling by the percent energy utilization such that the solar
        % panel and wind turbine size is constant rather than the total production - as was
        % the case during the process optimization for each percent
        % curtailment
        delta_panelCost(:,y) = (panelCost(2:end,1).*(S(2:end,1)/100) - panelCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_turbineCost(:,y) = (turbineCost(2:end,1).*(S(2:end,1)/100) - turbineCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2SizeCost(:,y) = (H2SizeCost(2:end,1).*(S(2:end,1)/100) - H2SizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BattSizeCost(:,y) = (BattSizeCost(2:end,1).*(S(2:end,1)/100) - BattSizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2bufferCost(:,y) = (H2bufferCost(2:end,1).*(S(2:end,1)/100) - H2bufferCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BatteryCost(:,y) = (BatteryCost(2:end,1).*(S(2:end,1)/100) - BatteryCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        %save the energy cost for the location and year in the CP_value
        %structure
        CP_costs.topwind_Ecost(i,y) = E_cost;
    end
    % add the change in cost for the first step as well
    delta_panelCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_panelCost];
    delta_turbineCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_turbineCost];
    delta_H2SizeCost = [H2Sizes(1,:) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/100) * (1/8760) ; delta_H2SizeCost];
    delta_BattSizeCost = [BattSizes(1,:) * BattpowerCapital * annual * (1/100) * (1/8760) ; delta_BattSizeCost];
    delta_H2bufferCost = [H2buffers(1,:) * H2storageCapital * annual * (1/100) * (1/8760) ; delta_H2bufferCost];
    delta_BatteryCost = [batteries(1,:) * BattstorageCapital * annual * (1/100) * (1/8760) ; delta_BatteryCost];
    % find the change in cost for the steady-state energy utilization
    % (constant array)
    delta_H2SizeCost_SS = ones(100,6) * (612.5/640) * (1/1000) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/8760);
    delta_ASUSizeCost_SS = ones(100,6) * (22.5/640) * (1/1000) * (ASUCapital+N2compCapital) * (annual + HB_OnM) * (1/8760);
    delta_HBSizeCost_SS = ones(100,6) * (1/100) * (HBCapital) * (annual + HB_OnM) * (1/8760);
    %create table of the added cost for each of the process components as a
    %function of energy utilization. add to the CP_value structure
    costs = table(mean(delta_panelCost,2), ...
        mean(delta_turbineCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["panel_cost","turbine_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    CP_costs.topwind_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    CP_costs.topwind_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    CP_costs.topwind_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end
% 3.iii.b Calculate the utilization costs for top solar locations with only solar
clearvars -except CP_costs panelCapital panel_OnM turbineCapital turbine_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
results = CP_costs.topwind_solaronly;
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_turbineCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    panels = flip(results.panel{i},1);
    H2Sizes = flip(results.H2Size{i},1);
    N2Sizes = flip(results.N2Size{i},1);
    BattSizes = flip(results.BattSize{i},1);
    H2buffers = flip(results.H2buffer{i},1);
    N2buffers = flip(results.N2buffer{i},1);
    batteries = flip(results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel and turbine size for full energy
        %utilization
        E_cost = (results.panel{i}(1,y) * (annual * panelCapital + (1/1000)*panel_OnM)) /(876000);
        % find the cost of each process component for each percent energy
        % utilization
        panelCost = panels(:,y) * (annual * panelCapital + (1/1000)*panel_OnM);
        H2SizeCost = H2Sizes(:,y) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM);
        BattSizeCost = BattSizes(:,y) * BattpowerCapital * annual;
        H2bufferCost = H2buffers(:,y) * H2storageCapital * annual;
        BatteryCost = batteries(:,y) * BattstorageCapital * annual;
        % find the change in the cost of each process component after
        % rescaling by the percent energy utilization such that the solar
        % panel size is constant rather than the total production - as was
        % the case during the process optimization for each percent
        % curtailment
        delta_panelCost(:,y) = (panelCost(2:end,1).*(S(2:end,1)/100) - panelCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2SizeCost(:,y) = (H2SizeCost(2:end,1).*(S(2:end,1)/100) - H2SizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BattSizeCost(:,y) = (BattSizeCost(2:end,1).*(S(2:end,1)/100) - BattSizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2bufferCost(:,y) = (H2bufferCost(2:end,1).*(S(2:end,1)/100) - H2bufferCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BatteryCost(:,y) = (BatteryCost(2:end,1).*(S(2:end,1)/100) - BatteryCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        %save the energy cost for the location and year in the CP_value
        %structure
        CP_costs.topwind_solaronly_Ecost(i,y) = E_cost;
    end
    % add the change in cost for the first step as well
    delta_panelCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_panelCost];
    delta_H2SizeCost = [H2Sizes(1,:) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/100) * (1/8760) ; delta_H2SizeCost];
    delta_BattSizeCost = [BattSizes(1,:) * BattpowerCapital * annual * (1/100) * (1/8760) ; delta_BattSizeCost];
    delta_H2bufferCost = [H2buffers(1,:) * H2storageCapital * annual * (1/100) * (1/8760) ; delta_H2bufferCost];
    delta_BatteryCost = [batteries(1,:) * BattstorageCapital * annual * (1/100) * (1/8760) ; delta_BatteryCost];
    % find the change in cost for the steady-state energy utilization
    % (constant array)
    delta_H2SizeCost_SS = ones(100,6) * (612.5/640) * (1/1000) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/8760);
    delta_ASUSizeCost_SS = ones(100,6) * (22.5/640) * (1/1000) * (ASUCapital+N2compCapital) * (annual + HB_OnM) * (1/8760);
    delta_HBSizeCost_SS = ones(100,6) * (1/100) * (HBCapital) * (annual + HB_OnM) * (1/8760);
    %create table of the added cost for each of the process components as a
    %function of energy utilization. add to the CP_value structure
    costs = table(mean(delta_panelCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["panel_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    CP_costs.topwind_solaronly_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    CP_costs.topwind_solaronly_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    CP_costs.topwind_solaronly_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end
% 3.iii.c Calculate the utilization costs for top wind locations with only wind
clearvars -except CP_costs panelCapital panel_OnM turbineCapital turbine_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
results = CP_costs.topwind_windonly;
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_turbineCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    turbines = flip(results.turbine{i},1);
    H2Sizes = flip(results.H2Size{i},1);
    N2Sizes = flip(results.N2Size{i},1);
    BattSizes = flip(results.BattSize{i},1);
    H2buffers = flip(results.H2buffer{i},1);
    N2buffers = flip(results.N2buffer{i},1);
    batteries = flip(results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel and turbine size for full energy
        %utilization
        E_cost = ((results.turbine{i}(1,y)) * (1 + turbine_OnM) * annual * turbineCapital ) /(876000);
        % find the cost of each process component for each percent energy
        % utilization
        turbineCost = turbines(:,y) * (1 + turbine_OnM) * annual * turbineCapital;
        H2SizeCost = H2Sizes(:,y) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM);
        BattSizeCost = BattSizes(:,y) * BattpowerCapital * annual;
        H2bufferCost = H2buffers(:,y) * H2storageCapital * annual;
        BatteryCost = batteries(:,y) * BattstorageCapital * annual;
        % find the change in the cost of each process component after
        % rescaling by the percent energy utilization such that the wind
        % turbine size is constant rather than the total production - as was
        % the case during the process optimization for each percent
        % curtailment
        delta_turbineCost(:,y) = (turbineCost(2:end,1).*(S(2:end,1)/100) - turbineCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2SizeCost(:,y) = (H2SizeCost(2:end,1).*(S(2:end,1)/100) - H2SizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BattSizeCost(:,y) = (BattSizeCost(2:end,1).*(S(2:end,1)/100) - BattSizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2bufferCost(:,y) = (H2bufferCost(2:end,1).*(S(2:end,1)/100) - H2bufferCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BatteryCost(:,y) = (BatteryCost(2:end,1).*(S(2:end,1)/100) - BatteryCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        %save the energy cost for the location and year in the CP_value
        %structure
        CP_costs.topwind_windonly_Ecost(i,y) = E_cost;
    end
    % add the change in cost for the first step as well
    delta_turbineCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_turbineCost];
    delta_H2SizeCost = [H2Sizes(1,:) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/100) * (1/8760) ; delta_H2SizeCost];
    delta_BattSizeCost = [BattSizes(1,:) * BattpowerCapital * annual * (1/100) * (1/8760) ; delta_BattSizeCost];
    delta_H2bufferCost = [H2buffers(1,:) * H2storageCapital * annual * (1/100) * (1/8760) ; delta_H2bufferCost];
    delta_BatteryCost = [batteries(1,:) * BattstorageCapital * annual * (1/100) * (1/8760) ; delta_BatteryCost];
    % find the change in cost for the steady-state energy utilization
    % (constant array)
    delta_H2SizeCost_SS = ones(100,6) * (612.5/640) * (1/1000) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/8760);
    delta_ASUSizeCost_SS = ones(100,6) * (22.5/640) * (1/1000) * (ASUCapital+N2compCapital) * (annual + HB_OnM) * (1/8760);
    delta_HBSizeCost_SS = ones(100,6) * (1/100) * (HBCapital) * (annual + HB_OnM) * (1/8760);
    %create table of the added cost for each of the process components as a
    %function of energy utilization. add to the CP_value structure
    costs = table(mean(delta_turbineCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["turbine_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    CP_costs.topwind_windonly_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    CP_costs.topwind_windonly_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    CP_costs.topwind_windonly_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end

% 3.iv.a Calculate the utilization costs for average locations combined solar and
% wind
clearvars -except CP_costs panelCapital panel_OnM turbineCapital turbine_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
results = CP_costs.average_results;
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_turbineCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    panels = flip(results.panel{i},1);
    turbines = flip(results.turbine{i},1);
    H2Sizes = flip(results.H2Size{i},1);
    N2Sizes = flip(results.N2Size{i},1);
    BattSizes = flip(results.BattSize{i},1);
    H2buffers = flip(results.H2buffer{i},1);
    N2buffers = flip(results.N2buffer{i},1);
    batteries = flip(results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel and turbine size for full energy
        %utilization
        E_cost = (results.panel{i}(1,y) * (annual * panelCapital + (1/1000)*panel_OnM) ...
            +(results.turbine{i}(1,y)) * (1 + turbine_OnM) * annual * turbineCapital ) /(876000);
        % find the cost of each process component for each percent energy
        % utilization
        panelCost = panels(:,y) * (annual * panelCapital + (1/1000)*panel_OnM);
        turbineCost = turbines(:,y) * (1 + turbine_OnM) * annual * turbineCapital;
        H2SizeCost = H2Sizes(:,y) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM);
        BattSizeCost = BattSizes(:,y) * BattpowerCapital * annual;
        H2bufferCost = H2buffers(:,y) * H2storageCapital * annual;
        BatteryCost = batteries(:,y) * BattstorageCapital * annual;
        % find the change in the cost of each process component after
        % rescaling by the percent energy utilization such that the solar
        % panel and wind turbine size is constant rather than the total production - as was
        % the case during the process optimization for each percent
        % curtailment
        delta_panelCost(:,y) = (panelCost(2:end,1).*(S(2:end,1)/100) - panelCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_turbineCost(:,y) = (turbineCost(2:end,1).*(S(2:end,1)/100) - turbineCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2SizeCost(:,y) = (H2SizeCost(2:end,1).*(S(2:end,1)/100) - H2SizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BattSizeCost(:,y) = (BattSizeCost(2:end,1).*(S(2:end,1)/100) - BattSizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2bufferCost(:,y) = (H2bufferCost(2:end,1).*(S(2:end,1)/100) - H2bufferCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BatteryCost(:,y) = (BatteryCost(2:end,1).*(S(2:end,1)/100) - BatteryCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        %save the energy cost for the location and year in the CP_value
        %structure
        CP_costs.average_Ecost(i,y) = E_cost;
    end
    % add the change in cost for the first step as well
    delta_panelCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_panelCost];
    delta_turbineCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_turbineCost];
    delta_H2SizeCost = [H2Sizes(1,:) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/100) * (1/8760) ; delta_H2SizeCost];
    delta_BattSizeCost = [BattSizes(1,:) * BattpowerCapital * annual * (1/100) * (1/8760) ; delta_BattSizeCost];
    delta_H2bufferCost = [H2buffers(1,:) * H2storageCapital * annual * (1/100) * (1/8760) ; delta_H2bufferCost];
    delta_BatteryCost = [batteries(1,:) * BattstorageCapital * annual * (1/100) * (1/8760) ; delta_BatteryCost];
    % find the change in cost for the steady-state energy utilization
    % (constant array)
    delta_H2SizeCost_SS = ones(100,6) * (612.5/640) * (1/1000) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/8760);
    delta_ASUSizeCost_SS = ones(100,6) * (22.5/640) * (1/1000) * (ASUCapital+N2compCapital) * (annual + HB_OnM) * (1/8760);
    delta_HBSizeCost_SS = ones(100,6) * (1/100) * (HBCapital) * (annual + HB_OnM) * (1/8760);
    %create table of the added cost for each of the process components as a
    %function of energy utilization. add to the CP_value structure
    costs = table(mean(delta_panelCost,2), ...
        mean(delta_turbineCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["panel_cost","turbine_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    CP_costs.average_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    CP_costs.average_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    CP_costs.average_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end
% 3.iv.b Calculate the utilization costs for average locations with only solar
clearvars -except CP_costs panelCapital panel_OnM turbineCapital turbine_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
results = CP_costs.average_solaronly;
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_turbineCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    panels = flip(results.panel{i},1);
    H2Sizes = flip(results.H2Size{i},1);
    N2Sizes = flip(results.N2Size{i},1);
    BattSizes = flip(results.BattSize{i},1);
    H2buffers = flip(results.H2buffer{i},1);
    N2buffers = flip(results.N2buffer{i},1);
    batteries = flip(results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel and turbine size for full energy
        %utilization
        E_cost = (results.panel{i}(1,y) * (annual * panelCapital + (1/1000)*panel_OnM)) /(876000);
        % find the cost of each process component for each percent energy
        % utilization
        panelCost = panels(:,y) * (annual * panelCapital + (1/1000)*panel_OnM);
        H2SizeCost = H2Sizes(:,y) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM);
        BattSizeCost = BattSizes(:,y) * BattpowerCapital * annual;
        H2bufferCost = H2buffers(:,y) * H2storageCapital * annual;
        BatteryCost = batteries(:,y) * BattstorageCapital * annual;
        % find the change in the cost of each process component after
        % rescaling by the percent energy utilization such that the solar
        % panel size is constant rather than the total production - as was
        % the case during the process optimization for each percent
        % curtailment
        delta_panelCost(:,y) = (panelCost(2:end,1).*(S(2:end,1)/100) - panelCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2SizeCost(:,y) = (H2SizeCost(2:end,1).*(S(2:end,1)/100) - H2SizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BattSizeCost(:,y) = (BattSizeCost(2:end,1).*(S(2:end,1)/100) - BattSizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2bufferCost(:,y) = (H2bufferCost(2:end,1).*(S(2:end,1)/100) - H2bufferCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BatteryCost(:,y) = (BatteryCost(2:end,1).*(S(2:end,1)/100) - BatteryCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        %save the energy cost for the location and year in the CP_value
        %structure
        CP_costs.average_solaronly_Ecost(i,y) = E_cost;
    end
    % add the change in cost for the first step as well
    delta_panelCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_panelCost];
    delta_H2SizeCost = [H2Sizes(1,:) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/100) * (1/8760) ; delta_H2SizeCost];
    delta_BattSizeCost = [BattSizes(1,:) * BattpowerCapital * annual * (1/100) * (1/8760) ; delta_BattSizeCost];
    delta_H2bufferCost = [H2buffers(1,:) * H2storageCapital * annual * (1/100) * (1/8760) ; delta_H2bufferCost];
    delta_BatteryCost = [batteries(1,:) * BattstorageCapital * annual * (1/100) * (1/8760) ; delta_BatteryCost];
    % find the change in cost for the steady-state energy utilization
    % (constant array)
    delta_H2SizeCost_SS = ones(100,6) * (612.5/640) * (1/1000) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/8760);
    delta_ASUSizeCost_SS = ones(100,6) * (22.5/640) * (1/1000) * (ASUCapital+N2compCapital) * (annual + HB_OnM) * (1/8760);
    delta_HBSizeCost_SS = ones(100,6) * (1/100) * (HBCapital) * (annual + HB_OnM) * (1/8760);
    %create table of the added cost for each of the process components as a
    %function of energy utilization. add to the CP_value structure
    costs = table(mean(delta_panelCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["panel_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    CP_costs.average_solaronly_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    CP_costs.average_solaronly_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    CP_costs.average_solaronly_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end
% 3.iv.c Calculate the utilization costs for average locations with only wind
clearvars -except CP_costs panelCapital panel_OnM turbineCapital turbine_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
results = CP_costs.average_windonly;
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_turbineCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    turbines = flip(results.turbine{i},1);
    H2Sizes = flip(results.H2Size{i},1);
    N2Sizes = flip(results.N2Size{i},1);
    BattSizes = flip(results.BattSize{i},1);
    H2buffers = flip(results.H2buffer{i},1);
    N2buffers = flip(results.N2buffer{i},1);
    batteries = flip(results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel and turbine size for full energy
        %utilization
        E_cost = ((results.turbine{i}(1,y)) * (1 + turbine_OnM) * annual * turbineCapital ) /(876000);
        % find the cost of each process component for each percent energy
        % utilization
        turbineCost = turbines(:,y) * (1 + turbine_OnM) * annual * turbineCapital;
        H2SizeCost = H2Sizes(:,y) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM);
        BattSizeCost = BattSizes(:,y) * BattpowerCapital * annual;
        H2bufferCost = H2buffers(:,y) * H2storageCapital * annual;
        BatteryCost = batteries(:,y) * BattstorageCapital * annual;
        % find the change in the cost of each process component after
        % rescaling by the percent energy utilization such that the wind
        % turbine size is constant rather than the total production - as was
        % the case during the process optimization for each percent
        % curtailment
        delta_turbineCost(:,y) = (turbineCost(2:end,1).*(S(2:end,1)/100) - turbineCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2SizeCost(:,y) = (H2SizeCost(2:end,1).*(S(2:end,1)/100) - H2SizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BattSizeCost(:,y) = (BattSizeCost(2:end,1).*(S(2:end,1)/100) - BattSizeCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_H2bufferCost(:,y) = (H2bufferCost(2:end,1).*(S(2:end,1)/100) - H2bufferCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        delta_BatteryCost(:,y) = (BatteryCost(2:end,1).*(S(2:end,1)/100) - BatteryCost(1:end-1,1).*(S(1:end-1,1)/100)) / 8760; %$ per MWh added
        %save the energy cost for the location and year in the CP_value
        %structure
        CP_costs.average_windonly_Ecost(i,y) = E_cost;
    end
    % add the change in cost for the first step as well
    delta_turbineCost = [[NaN,NaN,NaN,NaN,NaN,NaN]; delta_turbineCost];
    delta_H2SizeCost = [H2Sizes(1,:) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/100) * (1/8760) ; delta_H2SizeCost];
    delta_BattSizeCost = [BattSizes(1,:) * BattpowerCapital * annual * (1/100) * (1/8760) ; delta_BattSizeCost];
    delta_H2bufferCost = [H2buffers(1,:) * H2storageCapital * annual * (1/100) * (1/8760) ; delta_H2bufferCost];
    delta_BatteryCost = [batteries(1,:) * BattstorageCapital * annual * (1/100) * (1/8760) ; delta_BatteryCost];
    % find the change in cost for the steady-state energy utilization
    % (constant array)
    delta_H2SizeCost_SS = ones(100,6) * (612.5/640) * (1/1000) * (electrolyserCapital+H2compCapital) * (annual + elect_OnM) * (1/8760);
    delta_ASUSizeCost_SS = ones(100,6) * (22.5/640) * (1/1000) * (ASUCapital+N2compCapital) * (annual + HB_OnM) * (1/8760);
    delta_HBSizeCost_SS = ones(100,6) * (1/100) * (HBCapital) * (annual + HB_OnM) * (1/8760);
    %create table of the added cost for each of the process components as a
    %function of energy utilization. add to the CP_value structure
    costs = table(mean(delta_turbineCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["turbine_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    CP_costs.average_windonly_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    CP_costs.average_windonly_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    CP_costs.average_windonly_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end
save('CP_costs','CP_costs')
clearvars

% 4.i Calculate the amount of solar and wind utilization in each month for
% combined solar and wind at average locations
load('CP_base_wind')
load('CP_base_solar')
load('CP_curt')
load('CP_costs')
load('CP_nocurt')

S = linspace(1,100)'; %list is a scaling array for percent curtailment. 
% It is used to resize all of the process units for each level of curtailment 
% because the optimization assumes constant total production rather than 
% constant solar panel and wind turbine size
market_value = 1000 * (17/1000000)*(3600000/640); %$ per MWh. given a price of green ammonia as $1000/tonne
PReq = 100; %nomial size 100 MW
results = CP_costs.average_results;
for i = 1:10 % cycle through all 10 locations
    %load the process design data from the optimization of 0-99%
    %curtailment
    panels = flip(results.panel{i},1);
    turbines = flip(results.turbine{i},1);
    H2Sizes = flip(results.H2Size{i},1); % GW
    BattSizes = flip(results.BattSize{i},1); % GW
    H2buffers = flip(results.H2buffer{i},1) * (0.2045) * 0.000278; % GWh
    batteries = flip(results.batterykWh{i},1) * (1/1E6); % GWh
    battery_uses = (results.batteryUse{i}); % GW
    % extract the panel power and turbine power and the optimal amount of curtailment at the
    % location and year
    loc = CP_costs.average_loc(i,1);
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc;
    panelPowers = CP_base_solar.panelPower(ind);
    turbinePowers = CP_base_wind.turbinePower(ind);
    opt_curtailment = round(CP_curt.curtailment(ind)*100);
    
    % initialize a table of monthly results for each of the six years.
    % Results are organized for solar and wind separately. 
    month = ["Value_Type","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
    varTypes = ["string","double","double","double","double","double","double","double","double","double","double","double","double"];
    monthlyResults = table('Size', [77 numel(month)], 'VariableTypes',varTypes,'VariableNames', month);
    monthlyResults.Value_Type = [
        "Year 1 Total";"Year 2 Total";"Year 3 Total";"Year 4 Total";"Year 5 Total";"Year 6 Total";"Average Total"; ...
        "Year 1 Solar Energy";"Year 2 Solar Energy";"Year 3 Solar Energy";"Year 4 Solar Energy";"Year 5 Solar Energy";"Year 6 Solar Energy";"Average Solar Energy"; ...
        "Solar Year 1 Used";"Solar Year 2 Used";"Solar Year 3 Used";"Solar Year 4 Used";"Solar Year 5 Used";"Solar Year 6 Used";"Solar Average Used";...
        "Solar Year 1 Wasted";"Solar Year 2 Wasted";"Solar Year 3 Wasted";"Solar Year 4 Wasted";"Solar Year 5 Wasted";"Solar Year 6 Wasted";"Solar Average Wasted";
        "Solar Year 1 LAEC";"Solar Year 2 LAEC";"Solar Year 3 LAEC";"Solar Year 4 LAEC";"Solar Year 5 LAEC";"Solar Year 6 LAEC";"Solar Average LAEC"; ...
        "Solar Year 1 Xused";"Solar Year 2 Xused";"Solar Year 3 Xused";"Solar Year 4 Xused";"Solar Year 5 Xused";"Solar Year 6 Xused";"Solar Average Xused"; ...
        "Year 1 Wind Energy";"Year 2 Wind Energy";"Year 3 Wind Energy";"Year 4 Wind Energy";"Year 5 Wind Energy";"Year 6 Wind Energy";"Average Wind Energy"; ...
        "Wind Year 1 Used";"Wind Year 2 Used";"Wind Year 3 Used";"Wind Year 4 Used";"Wind Year 5 Used";"Wind Year 6 Used";"Wind Average Used";...
        "Wind Year 1 Wasted";"Wind Year 2 Wasted";"Wind Year 3 Wasted";"Wind Year 4 Wasted";"Wind Year 5 Wasted";"Wind Year 6 Wasted";"Wind Average Wasted";
        "Wind Year 1 LAEC";"Wind Year 2 LAEC";"Wind Year 3 LAEC";"Wind Year 4 LAEC";"Wind Year 5 LAEC";"Wind Year 6 LAEC";"Wind Average LAEC"; ...
        "Wind Year 1 Xused";"Wind Year 2 Xused";"Wind Year 3 Xused";"Wind Year 4 Xused";"Wind Year 5 Xused";"Wind Year 6 Xused";"Wind Average Xused"];

    for y = 1:6 %cycle through each of the years
        battery_uses_year = flip(battery_uses{1,y},1); %extract and flip the optimization result for the battery usage over a year (includes all 6 years at a location)
        %initialize arrays for the usage of energy in each monthly for each of the optimizations between 0-99% curtailment 
        % arrays are for solar and wind separately
        S_jan_used = zeros(100,1);
        S_feb_used = zeros(100,1);
        S_mar_used = zeros(100,1);
        S_apr_used = zeros(100,1);
        S_may_used = zeros(100,1);
        S_jun_used = zeros(100,1);
        S_jul_used = zeros(100,1);
        S_aug_used = zeros(100,1);
        S_sep_used = zeros(100,1);
        S_oct_used = zeros(100,1);
        S_nov_used = zeros(100,1);
        S_dec_used = zeros(100,1);
        W_jan_used = zeros(100,1);
        W_feb_used = zeros(100,1);
        W_mar_used = zeros(100,1);
        W_apr_used = zeros(100,1);
        W_may_used = zeros(100,1);
        W_jun_used = zeros(100,1);
        W_jul_used = zeros(100,1);
        W_aug_used = zeros(100,1);
        W_sep_used = zeros(100,1);
        W_oct_used = zeros(100,1);
        W_nov_used = zeros(100,1);
        W_dec_used = zeros(100,1);
         
        % Start of sub-section to find the optimal amount of solar and wind
        % curtailment respectively at the optimal total curtailment
        k = 100 - (opt_curtailment(y));
        buffCap = H2buffers(k,y) * S(k)/100; %set the buffer capacity based on the optimization result scale by the percent curtailment
        Curtailed = zeros(8760,1); % initialize array of curtailment for each hour over a year in GWh
        demand = single((S(k)/100) * 0.1 * (612.5/640) * ones(8760,1)); %set the hourly demand matrix in GW, scaled according to the percent curtailment
        power = single(panels(end,y) * panelPowers{y} * (1/1E9) + turbines(end,y) * turbinePowers{y} * (1/1E9)); %initialize the hourly power profile
        fract_solar = (panels(end,y) * panelPowers{y})./(panels(end,y)*panelPowers{y} + turbines(end,y)*turbinePowers{y});
        fract_solar(isnan(fract_solar)) = 0;
        Solar_E = [sum(power(1:744).*fract_solar(1:744)),...
            sum(power(745:1416).*fract_solar(745:1416)),...
            sum(power(1417:2160).*fract_solar(1417:2160)),...
            sum(power(2161:2880).*fract_solar(2161:2880)),...
            sum(power(2881:3624).*fract_solar(2881:3624)),...
            sum(power(3625:4344).*fract_solar(3625:4344)),...
            sum(power(4345:5088).*fract_solar(4345:5088)),...
            sum(power(5089:5832).*fract_solar(5089:5832)),...
            sum(power(5833:6552).*fract_solar(5833:6552)),...
            sum(power(6553:7296).*fract_solar(6553:7296)),...
            sum(power(7297:8016).*fract_solar(7297:8016)),...
            sum(power(8017:8760).*fract_solar(8017:8760))];
        Wind_E = [sum(power(1:744).*(1-fract_solar(1:744))),...
            sum(power(745:1416).*(1-fract_solar(745:1416))),...
            sum(power(1417:2160).*(1-fract_solar(1417:2160))),...
            sum(power(2161:2880).*(1-fract_solar(2161:2880))),...
            sum(power(2881:3624).*(1-fract_solar(2881:3624))),...
            sum(power(3625:4344).*(1-fract_solar(3625:4344))),...
            sum(power(4345:5088).*(1-fract_solar(4345:5088))),...
            sum(power(5089:5832).*(1-fract_solar(5089:5832))),...
            sum(power(5833:6552).*(1-fract_solar(5833:6552))),...
            sum(power(6553:7296).*(1-fract_solar(6553:7296))),...
            sum(power(7297:8016).*(1-fract_solar(7297:8016))),...
            sum(power(8017:8760).*(1-fract_solar(8017:8760)))];

        power_orig = power; %set the original power profile to use later
        battery_use = battery_uses_year{k,1}; %extract the particular hourly battery use over a year for a given location and year
        battery_use = (battery_use + (27.5/640)*PReq*(1/1000)) * S(k)/100; %baseline electrical requirement added to battery usage, making it electricity extraction, but keep name as battery
        battery_use(battery_use < 0 ) = 0; %if battery use is negative, set it to zero.
        power = power - battery_use; %power supply after removing the electricity extraction
        elect_curtail = power - H2Sizes(k,y)*S(k)/100; %calculate the amount of electricity curtailed due to it being above the rated capacity of the electrolyser
        elect_curtail(elect_curtail < 0) = 0; %remove entries less than zeros (keep only those where electricity must be curtailed)
        Curtailed = Curtailed + elect_curtail; %add the curtailed electricity to the curtailment array
        power(power > (H2Sizes(k,y)*S(k)/100)) = H2Sizes(k,y)*S(k)/100; %adjust the power supply according to the maximum power set by the electrolyser capacity
        cum_energy = [0;cumsum(power)]; % create a vector of cumulative energy supply for each hour over a year
        cum_demand = [0;cumsum(demand)]; % create a vector of the cumulative energy demand (as hydrogen) for each hour over a year
        Current_def = -1 * (cum_energy(1) - cum_energy(end)) + (cum_demand(1) - cum_demand(end)); % calculate the initial overall energy balance (should be positive)
        
        active_vector = ones(8760,1); %create a vector which describes which hours during the year are "active", meaning more energy during the hour can still be curtailed
        dummy_power = power; %initialize dummy array of power supply
        cum_dummy_power = [0;cumsum(dummy_power)]; % calculate cumulative dummy power

        while Current_def > 0.001 % continue to curtail energy until the overall energy balance over a year is clost to zero
            Avail_power = dummy_power .* active_vector; % calculate an array of the available power
            cum_Avail_power = [0;cumsum(Avail_power)]; % array of cumulative available power
            % create a matrix of the maxaimum percent curtailment which
            % can occur between the starting hour (row entry) and the
            % ending hour (column entry) before the energy deficit over
            % that period is equal to the hydrogen storage capacity
            max_curtail = ((-1 * triu(cum_dummy_power - cum_dummy_power') + triu(cum_demand - cum_demand')) + buffCap*ones(8761,8761))./(-1 * triu(cum_Avail_power - cum_Avail_power'));
                
            max_curtail(max_curtail<=0) = NaN; % remove entries less than zero
            curt = min(min(max_curtail)); % find the overall maximum amount of curtailment before the hydrogen storage capacity is reached
            % if the maximum curtailment is less than the amount of
            % curtailment which will cause the energy balance over the
            % year to close, record the starting and ending indexes.
            % Else, set the percent curtailment the amount which will
            % close the energy balance.    
            if curt < (Current_def/(-1*cum_Avail_power(1) + cum_Avail_power(end)))
                [r,c] = find(max_curtail == curt);
            else
                curt = Current_def/(-1+cum_Avail_power(1) + cum_Avail_power(end));
                r = [0];
                c = [0];
            end
            % add the hourly energy curtailed in this iteration to the
            % array of total curtailed energy
            Curtailed = Curtailed + Avail_power * curt;
            dummy_power = dummy_power - (Avail_power * curt); %update the dummy power array after removing curtailed energy
            cum_dummy_power = [0;cumsum(dummy_power)]; % update cumulative power array
            % if the percent curtailment is not set by closing the
            % yearly energy balance, then remove the active indexes
            % which reached their maximum level of curtailment (the
            % deficit over that period is equal to the hydrogen storage
            % capacity)   
            if r(1,1)~=0
                for n = 1:height(r)
                    active_vector(r(n,1):c(n,1)-1,1) = zeros((c(n,1)-r(n,1)),1);
                end
            end
            %update the overall energy balance over a year
            Current_def = -1 * (cum_dummy_power(1) - cum_dummy_power(end)) + (cum_demand(1) - cum_demand(end));
        end
        S_curt = round(100 * sum(Curtailed .* fract_solar)/sum(power_orig.*fract_solar));
        W_curt = round(100 * sum(Curtailed .* (1-fract_solar))/sum(power_orig.*(1-fract_solar)));
        % End of sub-section to find the optimal amount of solar and wind
        % curtailment respectively at the optimal total curtailment
        
        k = 1;
        parfor k = 1:99 %cycle through each of 1-99% curtailment
            disp(k)
            buffCap = H2buffers(k,y) * S(k)/100; %set the buffer capacity based on the optimization result scale by the percent curtailment
            Curtailed = zeros(8760,1); % initialize array of curtailment for each hour over a year in GWh
            demand = single((S(k)/100) * 0.1 * (612.5/640) * ones(8760,1)); %set the hourly demand matrix in GW, scaled according to the percent curtailment
            power = single(panels(end,y) * panelPowers{y} * (1/1E9) + turbines(end,y) * turbinePowers{y} * (1/1E9)); %initialize the hourly power profile
            fract_solar = (panels(end,y) * panelPowers{y})./(panels(end,y)*panelPowers{y} + turbines(end,y)*turbinePowers{y})
            fract_solar(isnan(fract_solar)) = 0;

            power_orig = power; %set the original power profile to use later
            battery_use = battery_uses_year{k,1}; %extract the particular hourly battery use over a year for a given location and year
            battery_use = (battery_use + (27.5/640)*PReq*(1/1000)) * S(k)/100; %baseline electrical requirement added to battery usage, making is electricity extraction, but keep name as battery
            battery_use(battery_use < 0 ) = 0; %if battery use is negative, set it to zero.
            power = power - battery_use; %power supply after removing the electricity extraction
            elect_curtail = power - H2Sizes(k,y)*S(k)/100; %calculate the amount of electricity curtailed due to it being above the rated capacity of the electrolyser
            elect_curtail(elect_curtail < 0) = 0; %remove entries less than zeros (keep only those where electricity must be curtailed)
            Curtailed = Curtailed + elect_curtail; %add the curtailed electricity to the curtailment array
            power(power > (H2Sizes(k,y)*S(k)/100)) = H2Sizes(k,y)*S(k)/100; %adjust the power supply according to the maximum power set by the electrolyser capacity
            cum_energy = [0;cumsum(power)]; % create a vector of cumulative energy supply for each hour over a year
            cum_demand = [0;cumsum(demand)]; % create a vector of the cumulative energy demand (as hydrogen) for each hour over a year
            Current_def = -1 * (cum_energy(1) - cum_energy(end)) + (cum_demand(1) - cum_demand(end)); % calculate the initial overall energy balance (should be positive)
            
            active_vector = ones(8760,1); %create a vector which describes which hours during the year are "active", meaning more energy during the hour can still be curtailed
            dummy_power = power; %initialize dummy array of power supply
            cum_dummy_power = [0;cumsum(dummy_power)]; % calculate cumulative dummy power

            while Current_def > 0.001 % continue to curtail energy until the overall energy balance over a year is close to zero
                Avail_power = dummy_power .* active_vector; % calculate an array of the available power
                cum_Avail_power = [0;cumsum(Avail_power)]; % array of cumulative available power
                % create a matrix of the maxaimum percent curtailment which
                % can occur between the starting hour (row entry) and the
                % ending hour (column entry) before the energy deficit over
                % that period is equal to the hydrogen storage capacity
                max_curtail = ((-1 * triu(cum_dummy_power - cum_dummy_power') + triu(cum_demand - cum_demand')) + buffCap*ones(8761,8761))./(-1 * triu(cum_Avail_power - cum_Avail_power'));
                
                max_curtail(max_curtail<=0) = NaN; % remove entries less than zero
                curt = min(min(max_curtail)); % find the overall maximum amount of curtailment before the hydrogen storage capacity is reached
                % if the maximum curtailment is less than the amount of
                % curtailment which will cause the energy balance over the
                % year to close, record the starting and ending indexes.
                % Else, set the percent curtailment the amount which will
                % close the energy balance. 
                if curt < (Current_def/(-1*cum_Avail_power(1) + cum_Avail_power(end)))
                    [r,c] = find(max_curtail == curt);
                else
                    curt = Current_def/(-1+cum_Avail_power(1) + cum_Avail_power(end));
                    r = [0];
                    c = [0];
                end
                % add the hourly energy curtailed in this iteration to the
                % array of total curtailed energy
                Curtailed = Curtailed + Avail_power * curt;
                dummy_power = dummy_power - (Avail_power * curt); %update the dummy power array after removing curtailed energy
                cum_dummy_power = [0;cumsum(dummy_power)]; % update cumulative power array
                % if the percent curtailment is not set by closing the
                % yearly energy balance, then remove the active indexes
                % which reached their maximum level of curtailment (the
                % deficit over that period is equal to the hydrogen storage
                % capacity)
                if r(1,1)~=0
                    for n = 1:height(r)
                        active_vector(r(n,1):c(n,1)-1,1) = zeros((c(n,1)-r(n,1)),1);
                    end
                end
                %update the overall energy balance over a year
                Current_def = -1 * (cum_dummy_power(1) - cum_dummy_power(end)) + (cum_demand(1) - cum_demand(end));
            end
            % save the percent of energy used in each month for the current
            % iteration of total percent curtailment
            S_jan_used(k+1,1) = sum((power_orig(1:744) - Curtailed(1:744)).*fract_solar(1:744))/sum(power_orig(1:744).*fract_solar(1:744));
            S_feb_used(k+1,1) = sum((power_orig(745:1416) - Curtailed(745:1416)).*fract_solar(745:1416))/sum(power_orig(745:1416).*fract_solar(745:1416));
            S_mar_used(k+1,1) = sum((power_orig(1417:2160) - Curtailed(1417:2160)).*fract_solar(1417:2160))/sum(power_orig(1417:2160).*fract_solar(1417:2160));
            S_apr_used(k+1,1) = sum((power_orig(2161:2880) - Curtailed(2161:2880)).*fract_solar(2161:2880))/sum(power_orig(2161:2880).*fract_solar(2161:2880));
            S_may_used(k+1,1) = sum((power_orig(2881:3624) - Curtailed(2881:3624)).*fract_solar(2881:3624))/sum(power_orig(2881:3624).*fract_solar(2881:3624));
            S_jun_used(k+1,1) = sum((power_orig(3625:4344) - Curtailed(3625:4344)).*fract_solar(3625:4344))/sum(power_orig(3625:4344).*fract_solar(3625:4344));
            S_jul_used(k+1,1) = sum((power_orig(4345:5088) - Curtailed(4345:5088)).*fract_solar(4345:5088))/sum(power_orig(4345:5088).*fract_solar(4345:5088));
            S_aug_used(k+1,1) = sum((power_orig(5089:5832) - Curtailed(5089:5832)).*fract_solar(5089:5832))/sum(power_orig(5089:5832).*fract_solar(5089:5832));
            S_sep_used(k+1,1) = sum((power_orig(5833:6552) - Curtailed(5833:6552)).*fract_solar(5833:6552))/sum(power_orig(5833:6552).*fract_solar(5833:6552));
            S_oct_used(k+1,1) = sum((power_orig(6553:7296) - Curtailed(6553:7296)).*fract_solar(6553:7296))/sum(power_orig(6553:7296).*fract_solar(6553:7296));
            S_nov_used(k+1,1) = sum((power_orig(7297:8016) - Curtailed(7297:8016)).*fract_solar(7297:8016))/sum(power_orig(7297:8016).*fract_solar(7297:8016));
            S_dec_used(k+1,1) = sum((power_orig(8017:8760) - Curtailed(8017:8760)).*fract_solar(8017:8760))/sum(power_orig(8017:8760).*fract_solar(8017:8760));
            
            W_jan_used(k+1,1) = sum((power_orig(1:744) - Curtailed(1:744)).*(1-fract_solar(1:744)))/sum(power_orig(1:744).*(1-fract_solar(1:744)));
            W_feb_used(k+1,1) = sum((power_orig(745:1416) - Curtailed(745:1416)).*(1-fract_solar(745:1416)))/sum(power_orig(745:1416).*(1-fract_solar(745:1416)));
            W_mar_used(k+1,1) = sum((power_orig(1417:2160) - Curtailed(1417:2160)).*(1-fract_solar(1417:2160)))/sum(power_orig(1417:2160).*(1-fract_solar(1417:2160)));
            W_apr_used(k+1,1) = sum((power_orig(2161:2880) - Curtailed(2161:2880)).*(1-fract_solar(2161:2880)))/sum(power_orig(2161:2880).*(1-fract_solar(2161:2880)));
            W_may_used(k+1,1) = sum((power_orig(2881:3624) - Curtailed(2881:3624)).*(1-fract_solar(2881:3624)))/sum(power_orig(2881:3624).*(1-fract_solar(2881:3624)));
            W_jun_used(k+1,1) = sum((power_orig(3625:4344) - Curtailed(3625:4344)).*(1-fract_solar(3625:4344)))/sum(power_orig(3625:4344).*(1-fract_solar(3625:4344)));
            W_jul_used(k+1,1) = sum((power_orig(4345:5088) - Curtailed(4345:5088)).*(1-fract_solar(4345:5088)))/sum(power_orig(4345:5088).*(1-fract_solar(4345:5088)));
            W_aug_used(k+1,1) = sum((power_orig(5089:5832) - Curtailed(5089:5832)).*(1-fract_solar(5089:5832)))/sum(power_orig(5089:5832).*(1-fract_solar(5089:5832)));
            W_sep_used(k+1,1) = sum((power_orig(5833:6552) - Curtailed(5833:6552)).*(1-fract_solar(5833:6552)))/sum(power_orig(5833:6552).*(1-fract_solar(5833:6552)));
            W_oct_used(k+1,1) = sum((power_orig(6553:7296) - Curtailed(6553:7296)).*(1-fract_solar(6553:7296)))/sum(power_orig(6553:7296).*(1-fract_solar(6553:7296)));
            W_nov_used(k+1,1) = sum((power_orig(7297:8016) - Curtailed(7297:8016)).*(1-fract_solar(7297:8016)))/sum(power_orig(7297:8016).*(1-fract_solar(7297:8016)));
            W_dec_used(k+1,1) = sum((power_orig(8017:8760) - Curtailed(8017:8760)).*(1-fract_solar(8017:8760)))/sum(power_orig(8017:8760).*(1-fract_solar(8017:8760)));
        end
        % create matrix of percent of energy used each month for each of
        % the total percent curtailments
        S_month_used = [S_jan_used,S_feb_used,S_mar_used,S_apr_used,S_may_used,S_jun_used,S_jul_used,S_aug_used,S_sep_used,S_oct_used,S_nov_used,S_dec_used];
        W_month_used = [W_jan_used,W_feb_used,W_mar_used,W_apr_used,W_may_used,W_jun_used,W_jul_used,W_aug_used,W_sep_used,W_oct_used,W_nov_used,W_dec_used];
        %difference between the curtailment in the month at a given curtailment level and the previous level
        % add the final difference in energy used for the case case of 0% curtailment
        S_usedDiff = [S_month_used(2:end,:) - S_month_used(1:end-1,:)]; %difference between the curtailment in the month at a given curtailment level and the previous level
        S_usedDiff = [S_usedDiff ; 1-sum(S_usedDiff)];
        W_usedDiff = [W_month_used(2:end,:) - W_month_used(1:end-1,:)]; %difference between the curtailment in the month at a given curtailment level and the previous level
        W_usedDiff = [W_usedDiff ; 1-sum(W_usedDiff)];
        % remove entries that at NaN due to there being no solar or wind
        % installed
        W_usedDiff(isnan(W_usedDiff)) = 0;
        S_usedDiff(isnan(S_usedDiff)) = 0;

        costStep = CP_costs.average_yearlyLCOU{i}(:,y); %cost of utilizing energy at each percent curtailment level
        LAEC = CP_costs.average_yearlyLAEC{i}(:,y); % total energy cost (except the energy supply cost) at each percent curtailment level
        S_E_cost = CP_costs.average_solaronly_Ecost(i,y); % solar energy supply cost in $/MWh
        W_E_cost = CP_costs.average_windonly_Ecost(i,y); % wind energy supply cost in $/MWh
        % total value of energy in the month is the sum of the difference in curtailment multiplied by the value of the energy wasted
        monthlyValueTotal = sum(S_usedDiff .* (market_value - costStep)) + sum(W_usedDiff .* (market_value - costStep)); %value of energy in the month is the sum of the difference in curtailment multiplied by the value of the energy wasted
        % value of the used energy in the monthly is value of the energy
        % used multiplied by the change in the amount of energy for each
        % percent curtailment, only up to the optimial percent curtailment,
        % normalized by the total energy used in the month
        S_monthlyValueUsed = sum(S_usedDiff(1:(100-opt_curtailment(y)),:) .* (market_value - costStep(1:(100-opt_curtailment(y)),:))) ./ sum(S_usedDiff(1:(100-opt_curtailment(y)),:));
        % similarly for the value of curtailed energy, only analysing the
        % percent curtailment between the optimal energy utilization and
        % 100% energy utilization
        S_monthlyValueWasted = sum(S_usedDiff((101-opt_curtailment(y)):end,:) .* (market_value - costStep((101-opt_curtailment(y)):end,:))) ./ sum(S_usedDiff((101-opt_curtailment(y)):end,:));
        % total cost of the energy used, including the energy supply cost
        % and the cost of utilizing energy
        S_monthlyCostUsed = sum(S_usedDiff(1:(100-opt_curtailment(y)),:) .* ((S_E_cost*(100/(100-S_curt))) + LAEC(1:(100-opt_curtailment(y)),:))) ./ sum(S_usedDiff(1:(100-opt_curtailment(y)),:));
        % total fraction of energy used in each month
        S_monthlyXUsed = sum(S_usedDiff(1:(100-opt_curtailment(y)),:));
        % repeat for the wind energy 
        W_monthlyValueUsed = sum(W_usedDiff(1:(100-opt_curtailment(y)),:) .* (market_value - costStep(1:(100-opt_curtailment(y)),:))) ./ sum(W_usedDiff(1:(100-opt_curtailment(y)),:));
        W_monthlyValueWasted = sum(W_usedDiff((101-opt_curtailment(y)):end,:) .* (market_value - costStep((101-opt_curtailment(y)):end,:))) ./ sum(W_usedDiff((101-opt_curtailment(y)):end,:));
        W_monthlyCostUsed = sum(W_usedDiff(1:(100-opt_curtailment(y)),:) .* ((W_E_cost*(100/(100-W_curt))) + LAEC(1:(100-opt_curtailment(y)),:))) ./ sum(W_usedDiff(1:(100-opt_curtailment(y)),:));
        W_monthlyXUsed = sum(W_usedDiff(1:(100-opt_curtailment(y)),:));
        % remove the NaN entries
        S_monthlyValueUsed(isnan(S_monthlyValueUsed)) = 0;
        S_monthlyValueWasted(isnan(S_monthlyValueWasted)) = 0;
        S_monthlyCostUsed(isnan(S_monthlyCostUsed)) = 0;
        W_monthlyValueUsed(isnan(W_monthlyValueUsed)) = 0;
        W_monthlyValueWasted(isnan(W_monthlyValueWasted)) = 0;
        W_monthlyCostUsed(isnan(W_monthlyCostUsed)) = 0;
        %save the above metrics in the monthlyResults table
        monthlyResults{y,2:end} = monthlyValueTotal; %add results for that year to the record table
        monthlyResults{y+7,2:end} = Solar_E;
        monthlyResults{y+14,2:end} = S_monthlyValueUsed;
        monthlyResults{y+21,2:end} = S_monthlyValueWasted;
        monthlyResults{y+28,2:end} = S_monthlyCostUsed;
        monthlyResults{y+35,2:end} = S_monthlyXUsed;
        monthlyResults{y+42,2:end} = Wind_E;
        monthlyResults{y+49,2:end} = W_monthlyValueUsed;
        monthlyResults{y+56,2:end} = W_monthlyValueWasted;
        monthlyResults{y+63,2:end} = W_monthlyCostUsed;
        monthlyResults{y+70,2:end} = W_monthlyXUsed;
    end
    % find the average of the metrics over the 6 years, removing the zero
    % entries (i.e. only consider the year as counting toward the average
    % if solar or wind were installed for that year)
    for m = 2:13
        monthlyResults{7,m} = mean(nonzeros(monthlyResults{1:6,m}));
        monthlyResults{14,m} = mean(nonzeros(monthlyResults{8:13,m}));
        monthlyResults{21,m} = mean(nonzeros(monthlyResults{15:20,m}));
        monthlyResults{28,m} = mean(nonzeros(monthlyResults{22:27,m}));
        monthlyResults{35,m} = mean(nonzeros(monthlyResults{29:34,m}));
        monthlyResults{42,m} = mean(nonzeros(monthlyResults{36:41,m}));
        monthlyResults{49,m} = mean(nonzeros(monthlyResults{43:48,m}));
        monthlyResults{56,m} = mean(nonzeros(monthlyResults{50:55,m}));
        monthlyResults{63,m} = mean(nonzeros(monthlyResults{57:62,m}));
        monthlyResults{70,m} = mean(nonzeros(monthlyResults{64:69,m}));
        monthlyResults{77,m} = mean(nonzeros(monthlyResults{71:76,m}));
    end
    CP_costs.average_monthlyValues{i} = monthlyResults;
end
save('CP_costs','CP_costs')
clearvars -except CP_costs CP_base_solar CP_base_wind CP_nocurt CP_curt

% 4.ii Calculate the amount of solar and wind utilization in each month for
% combined solar and wind at top solar locations
S = linspace(1,100)'; %list is a scaling array for percent curtailment. 
% It is used to resize all of the process units for each level of curtailment 
% because the optimization assumes constant total production rather than 
% constant solar panel and wind turbine size
market_value = 1000 * (17/1000000)*(3600000/640); %$ per MWh. given a price of green ammonia as $1000/tonne
PReq = 100; %nomial size 100 MW
results = CP_costs.topsolar_results;
for i = 1:10 % cycle through all 10 locations
    %load the process design data from the optimization of 0-99%
    %curtailment
    panels = flip(results.panel{i},1);
    turbines = flip(results.turbine{i},1);
    H2Sizes = flip(results.H2Size{i},1); % GW
    BattSizes = flip(results.BattSize{i},1); % GW
    H2buffers = flip(results.H2buffer{i},1) * (0.2045) * 0.000278; % GWh
    batteries = flip(results.batterykWh{i},1) * (1/1E6); % GWh
    battery_uses = (results.batteryUse{i}); % GW
    % extract the panel power and turbine power and the optimal amount of curtailment at the
    % location and year
    loc = CP_costs.topsolar_loc(i,1);
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc;
    panelPowers = CP_base_solar.panelPower(ind);
    turbinePowers = CP_base_wind.turbinePower(ind);
    opt_curtailment = round(CP_curt.curtailment(ind)*100);
    
    % initialize a table of monthly results for each of the six years.
    % Results are organized for solar and wind separately. 
    month = ["Value_Type","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
    varTypes = ["string","double","double","double","double","double","double","double","double","double","double","double","double"];
    monthlyResults = table('Size', [77 numel(month)], 'VariableTypes',varTypes,'VariableNames', month);
    monthlyResults.Value_Type = [
        "Year 1 Total";"Year 2 Total";"Year 3 Total";"Year 4 Total";"Year 5 Total";"Year 6 Total";"Average Total"; ...
        "Year 1 Solar Energy";"Year 2 Solar Energy";"Year 3 Solar Energy";"Year 4 Solar Energy";"Year 5 Solar Energy";"Year 6 Solar Energy";"Average Solar Energy"; ...
        "Solar Year 1 Used";"Solar Year 2 Used";"Solar Year 3 Used";"Solar Year 4 Used";"Solar Year 5 Used";"Solar Year 6 Used";"Solar Average Used";...
        "Solar Year 1 Wasted";"Solar Year 2 Wasted";"Solar Year 3 Wasted";"Solar Year 4 Wasted";"Solar Year 5 Wasted";"Solar Year 6 Wasted";"Solar Average Wasted";
        "Solar Year 1 LAEC";"Solar Year 2 LAEC";"Solar Year 3 LAEC";"Solar Year 4 LAEC";"Solar Year 5 LAEC";"Solar Year 6 LAEC";"Solar Average LAEC"; ...
        "Solar Year 1 Xused";"Solar Year 2 Xused";"Solar Year 3 Xused";"Solar Year 4 Xused";"Solar Year 5 Xused";"Solar Year 6 Xused";"Solar Average Xused"; ...
        "Year 1 Wind Energy";"Year 2 Wind Energy";"Year 3 Wind Energy";"Year 4 Wind Energy";"Year 5 Wind Energy";"Year 6 Wind Energy";"Average Wind Energy"; ...
        "Wind Year 1 Used";"Wind Year 2 Used";"Wind Year 3 Used";"Wind Year 4 Used";"Wind Year 5 Used";"Wind Year 6 Used";"Wind Average Used";...
        "Wind Year 1 Wasted";"Wind Year 2 Wasted";"Wind Year 3 Wasted";"Wind Year 4 Wasted";"Wind Year 5 Wasted";"Wind Year 6 Wasted";"Wind Average Wasted";
        "Wind Year 1 LAEC";"Wind Year 2 LAEC";"Wind Year 3 LAEC";"Wind Year 4 LAEC";"Wind Year 5 LAEC";"Wind Year 6 LAEC";"Wind Average LAEC"; ...
        "Wind Year 1 Xused";"Wind Year 2 Xused";"Wind Year 3 Xused";"Wind Year 4 Xused";"Wind Year 5 Xused";"Wind Year 6 Xused";"Wind Average Xused"];

    for y = 1:6 %cycle through each of the years
        battery_uses_year = flip(battery_uses{1,y},1); %extract and flip the optimization result for the battery usage over a year (includes all 6 years at a location)
        %initialize arrays for the usage of energy in each monthly for each of the optimizations between 0-99% curtailment 
        % arrays are for solar and wind separately
        S_jan_used = zeros(100,1);
        S_feb_used = zeros(100,1);
        S_mar_used = zeros(100,1);
        S_apr_used = zeros(100,1);
        S_may_used = zeros(100,1);
        S_jun_used = zeros(100,1);
        S_jul_used = zeros(100,1);
        S_aug_used = zeros(100,1);
        S_sep_used = zeros(100,1);
        S_oct_used = zeros(100,1);
        S_nov_used = zeros(100,1);
        S_dec_used = zeros(100,1);
        W_jan_used = zeros(100,1);
        W_feb_used = zeros(100,1);
        W_mar_used = zeros(100,1);
        W_apr_used = zeros(100,1);
        W_may_used = zeros(100,1);
        W_jun_used = zeros(100,1);
        W_jul_used = zeros(100,1);
        W_aug_used = zeros(100,1);
        W_sep_used = zeros(100,1);
        W_oct_used = zeros(100,1);
        W_nov_used = zeros(100,1);
        W_dec_used = zeros(100,1);
         
        % Start of sub-section to find the optimal amount of solar and wind
        % curtailment respectively at the optimal total curtailment
        k = 100 - (opt_curtailment(y));
        buffCap = H2buffers(k,y) * S(k)/100; %set the buffer capacity based on the optimization result scale by the percent curtailment
        Curtailed = zeros(8760,1); % initialize array of curtailment for each hour over a year in GWh
        demand = single((S(k)/100) * 0.1 * (612.5/640) * ones(8760,1)); %set the hourly demand matrix in GW, scaled according to the percent curtailment
        power = single(panels(end,y) * panelPowers{y} * (1/1E9) + turbines(end,y) * turbinePowers{y} * (1/1E9)); %initialize the hourly power profile
        fract_solar = (panels(end,y) * panelPowers{y})./(panels(end,y)*panelPowers{y} + turbines(end,y)*turbinePowers{y});
        fract_solar(isnan(fract_solar)) = 0;
        Solar_E = [sum(power(1:744).*fract_solar(1:744)),...
            sum(power(745:1416).*fract_solar(745:1416)),...
            sum(power(1417:2160).*fract_solar(1417:2160)),...
            sum(power(2161:2880).*fract_solar(2161:2880)),...
            sum(power(2881:3624).*fract_solar(2881:3624)),...
            sum(power(3625:4344).*fract_solar(3625:4344)),...
            sum(power(4345:5088).*fract_solar(4345:5088)),...
            sum(power(5089:5832).*fract_solar(5089:5832)),...
            sum(power(5833:6552).*fract_solar(5833:6552)),...
            sum(power(6553:7296).*fract_solar(6553:7296)),...
            sum(power(7297:8016).*fract_solar(7297:8016)),...
            sum(power(8017:8760).*fract_solar(8017:8760))];
        Wind_E = [sum(power(1:744).*(1-fract_solar(1:744))),...
            sum(power(745:1416).*(1-fract_solar(745:1416))),...
            sum(power(1417:2160).*(1-fract_solar(1417:2160))),...
            sum(power(2161:2880).*(1-fract_solar(2161:2880))),...
            sum(power(2881:3624).*(1-fract_solar(2881:3624))),...
            sum(power(3625:4344).*(1-fract_solar(3625:4344))),...
            sum(power(4345:5088).*(1-fract_solar(4345:5088))),...
            sum(power(5089:5832).*(1-fract_solar(5089:5832))),...
            sum(power(5833:6552).*(1-fract_solar(5833:6552))),...
            sum(power(6553:7296).*(1-fract_solar(6553:7296))),...
            sum(power(7297:8016).*(1-fract_solar(7297:8016))),...
            sum(power(8017:8760).*(1-fract_solar(8017:8760)))];

        power_orig = power; %set the original power profile to use later
        battery_use = battery_uses_year{k,1}; %extract the particular hourly battery use over a year for a given location and year
        battery_use = (battery_use + (27.5/640)*PReq*(1/1000)) * S(k)/100; %baseline electrical requirement added to battery usage, making it electricity extraction, but keep name as battery
        battery_use(battery_use < 0 ) = 0; %if battery use is negative, set it to zero.
        power = power - battery_use; %power supply after removing the electricity extraction
        elect_curtail = power - H2Sizes(k,y)*S(k)/100; %calculate the amount of electricity curtailed due to it being above the rated capacity of the electrolyser
        elect_curtail(elect_curtail < 0) = 0; %remove entries less than zeros (keep only those where electricity must be curtailed)
        Curtailed = Curtailed + elect_curtail; %add the curtailed electricity to the curtailment array
        power(power > (H2Sizes(k,y)*S(k)/100)) = H2Sizes(k,y)*S(k)/100; %adjust the power supply according to the maximum power set by the electrolyser capacity
        cum_energy = [0;cumsum(power)]; % create a vector of cumulative energy supply for each hour over a year
        cum_demand = [0;cumsum(demand)]; % create a vector of the cumulative energy demand (as hydrogen) for each hour over a year
        Current_def = -1 * (cum_energy(1) - cum_energy(end)) + (cum_demand(1) - cum_demand(end)); % calculate the initial overall energy balance (should be positive)
        
        active_vector = ones(8760,1); %create a vector which describes which hours during the year are "active", meaning more energy during the hour can still be curtailed
        dummy_power = power; %initialize dummy array of power supply
        cum_dummy_power = [0;cumsum(dummy_power)]; % calculate cumulative dummy power

        while Current_def > 0.001 % continue to curtail energy until the overall energy balance over a year is clost to zero
            Avail_power = dummy_power .* active_vector; % calculate an array of the available power
            cum_Avail_power = [0;cumsum(Avail_power)]; % array of cumulative available power
            % create a matrix of the maxaimum percent curtailment which
            % can occur between the starting hour (row entry) and the
            % ending hour (column entry) before the energy deficit over
            % that period is equal to the hydrogen storage capacity
            max_curtail = ((-1 * triu(cum_dummy_power - cum_dummy_power') + triu(cum_demand - cum_demand')) + buffCap*ones(8761,8761))./(-1 * triu(cum_Avail_power - cum_Avail_power'));
                
            max_curtail(max_curtail<=0) = NaN; % remove entries less than zero
            curt = min(min(max_curtail)); % find the overall maximum amount of curtailment before the hydrogen storage capacity is reached
            % if the maximum curtailment is less than the amount of
            % curtailment which will cause the energy balance over the
            % year to close, record the starting and ending indexes.
            % Else, set the percent curtailment the amount which will
            % close the energy balance.    
            if curt < (Current_def/(-1*cum_Avail_power(1) + cum_Avail_power(end)))
                [r,c] = find(max_curtail == curt);
            else
                curt = Current_def/(-1+cum_Avail_power(1) + cum_Avail_power(end));
                r = [0];
                c = [0];
            end
            % add the hourly energy curtailed in this iteration to the
            % array of total curtailed energy
            Curtailed = Curtailed + Avail_power * curt;
            dummy_power = dummy_power - (Avail_power * curt); %update the dummy power array after removing curtailed energy
            cum_dummy_power = [0;cumsum(dummy_power)]; % update cumulative power array
            % if the percent curtailment is not set by closing the
            % yearly energy balance, then remove the active indexes
            % which reached their maximum level of curtailment (the
            % deficit over that period is equal to the hydrogen storage
            % capacity)   
            if r(1,1)~=0
                for n = 1:height(r)
                    active_vector(r(n,1):c(n,1)-1,1) = zeros((c(n,1)-r(n,1)),1);
                end
            end
            %update the overall energy balance over a year
            Current_def = -1 * (cum_dummy_power(1) - cum_dummy_power(end)) + (cum_demand(1) - cum_demand(end));
        end
        S_curt = round(100 * sum(Curtailed .* fract_solar)/sum(power_orig.*fract_solar));
        W_curt = round(100 * sum(Curtailed .* (1-fract_solar))/sum(power_orig.*(1-fract_solar)));
        % End of sub-section to find the optimal amount of solar and wind
        % curtailment respectively at the optimal total curtailment
        
        k = 1;
        parfor k = 1:99 %cycle through each of 1-99% curtailment
            disp(k)
            buffCap = H2buffers(k,y) * S(k)/100; %set the buffer capacity based on the optimization result scale by the percent curtailment
            Curtailed = zeros(8760,1); % initialize array of curtailment for each hour over a year in GWh
            demand = single((S(k)/100) * 0.1 * (612.5/640) * ones(8760,1)); %set the hourly demand matrix in GW, scaled according to the percent curtailment
            power = single(panels(end,y) * panelPowers{y} * (1/1E9) + turbines(end,y) * turbinePowers{y} * (1/1E9)); %initialize the hourly power profile
            fract_solar = (panels(end,y) * panelPowers{y})./(panels(end,y)*panelPowers{y} + turbines(end,y)*turbinePowers{y})
            fract_solar(isnan(fract_solar)) = 0;

            power_orig = power; %set the original power profile to use later
            battery_use = battery_uses_year{k,1}; %extract the particular hourly battery use over a year for a given location and year
            battery_use = (battery_use + (27.5/640)*PReq*(1/1000)) * S(k)/100; %baseline electrical requirement added to battery usage, making is electricity extraction, but keep name as battery
            battery_use(battery_use < 0 ) = 0; %if battery use is negative, set it to zero.
            power = power - battery_use; %power supply after removing the electricity extraction
            elect_curtail = power - H2Sizes(k,y)*S(k)/100; %calculate the amount of electricity curtailed due to it being above the rated capacity of the electrolyser
            elect_curtail(elect_curtail < 0) = 0; %remove entries less than zeros (keep only those where electricity must be curtailed)
            Curtailed = Curtailed + elect_curtail; %add the curtailed electricity to the curtailment array
            power(power > (H2Sizes(k,y)*S(k)/100)) = H2Sizes(k,y)*S(k)/100; %adjust the power supply according to the maximum power set by the electrolyser capacity
            cum_energy = [0;cumsum(power)]; % create a vector of cumulative energy supply for each hour over a year
            cum_demand = [0;cumsum(demand)]; % create a vector of the cumulative energy demand (as hydrogen) for each hour over a year
            Current_def = -1 * (cum_energy(1) - cum_energy(end)) + (cum_demand(1) - cum_demand(end)); % calculate the initial overall energy balance (should be positive)
            
            active_vector = ones(8760,1); %create a vector which describes which hours during the year are "active", meaning more energy during the hour can still be curtailed
            dummy_power = power; %initialize dummy array of power supply
            cum_dummy_power = [0;cumsum(dummy_power)]; % calculate cumulative dummy power

            while Current_def > 0.001 % continue to curtail energy until the overall energy balance over a year is close to zero
                Avail_power = dummy_power .* active_vector; % calculate an array of the available power
                cum_Avail_power = [0;cumsum(Avail_power)]; % array of cumulative available power
                % create a matrix of the maxaimum percent curtailment which
                % can occur between the starting hour (row entry) and the
                % ending hour (column entry) before the energy deficit over
                % that period is equal to the hydrogen storage capacity
                max_curtail = ((-1 * triu(cum_dummy_power - cum_dummy_power') + triu(cum_demand - cum_demand')) + buffCap*ones(8761,8761))./(-1 * triu(cum_Avail_power - cum_Avail_power'));
                
                max_curtail(max_curtail<=0) = NaN; % remove entries less than zero
                curt = min(min(max_curtail)); % find the overall maximum amount of curtailment before the hydrogen storage capacity is reached
                % if the maximum curtailment is less than the amount of
                % curtailment which will cause the energy balance over the
                % year to close, record the starting and ending indexes.
                % Else, set the percent curtailment the amount which will
                % close the energy balance. 
                if curt < (Current_def/(-1*cum_Avail_power(1) + cum_Avail_power(end)))
                    [r,c] = find(max_curtail == curt);
                else
                    curt = Current_def/(-1+cum_Avail_power(1) + cum_Avail_power(end));
                    r = [0];
                    c = [0];
                end
                % add the hourly energy curtailed in this iteration to the
                % array of total curtailed energy
                Curtailed = Curtailed + Avail_power * curt;
                dummy_power = dummy_power - (Avail_power * curt); %update the dummy power array after removing curtailed energy
                cum_dummy_power = [0;cumsum(dummy_power)]; % update cumulative power array
                % if the percent curtailment is not set by closing the
                % yearly energy balance, then remove the active indexes
                % which reached their maximum level of curtailment (the
                % deficit over that period is equal to the hydrogen storage
                % capacity)
                if r(1,1)~=0
                    for n = 1:height(r)
                        active_vector(r(n,1):c(n,1)-1,1) = zeros((c(n,1)-r(n,1)),1);
                    end
                end
                %update the overall energy balance over a year
                Current_def = -1 * (cum_dummy_power(1) - cum_dummy_power(end)) + (cum_demand(1) - cum_demand(end));
            end
            % save the percent of energy used in each month for the current
            % iteration of total percent curtailment
            S_jan_used(k+1,1) = sum((power_orig(1:744) - Curtailed(1:744)).*fract_solar(1:744))/sum(power_orig(1:744).*fract_solar(1:744));
            S_feb_used(k+1,1) = sum((power_orig(745:1416) - Curtailed(745:1416)).*fract_solar(745:1416))/sum(power_orig(745:1416).*fract_solar(745:1416));
            S_mar_used(k+1,1) = sum((power_orig(1417:2160) - Curtailed(1417:2160)).*fract_solar(1417:2160))/sum(power_orig(1417:2160).*fract_solar(1417:2160));
            S_apr_used(k+1,1) = sum((power_orig(2161:2880) - Curtailed(2161:2880)).*fract_solar(2161:2880))/sum(power_orig(2161:2880).*fract_solar(2161:2880));
            S_may_used(k+1,1) = sum((power_orig(2881:3624) - Curtailed(2881:3624)).*fract_solar(2881:3624))/sum(power_orig(2881:3624).*fract_solar(2881:3624));
            S_jun_used(k+1,1) = sum((power_orig(3625:4344) - Curtailed(3625:4344)).*fract_solar(3625:4344))/sum(power_orig(3625:4344).*fract_solar(3625:4344));
            S_jul_used(k+1,1) = sum((power_orig(4345:5088) - Curtailed(4345:5088)).*fract_solar(4345:5088))/sum(power_orig(4345:5088).*fract_solar(4345:5088));
            S_aug_used(k+1,1) = sum((power_orig(5089:5832) - Curtailed(5089:5832)).*fract_solar(5089:5832))/sum(power_orig(5089:5832).*fract_solar(5089:5832));
            S_sep_used(k+1,1) = sum((power_orig(5833:6552) - Curtailed(5833:6552)).*fract_solar(5833:6552))/sum(power_orig(5833:6552).*fract_solar(5833:6552));
            S_oct_used(k+1,1) = sum((power_orig(6553:7296) - Curtailed(6553:7296)).*fract_solar(6553:7296))/sum(power_orig(6553:7296).*fract_solar(6553:7296));
            S_nov_used(k+1,1) = sum((power_orig(7297:8016) - Curtailed(7297:8016)).*fract_solar(7297:8016))/sum(power_orig(7297:8016).*fract_solar(7297:8016));
            S_dec_used(k+1,1) = sum((power_orig(8017:8760) - Curtailed(8017:8760)).*fract_solar(8017:8760))/sum(power_orig(8017:8760).*fract_solar(8017:8760));
            
            W_jan_used(k+1,1) = sum((power_orig(1:744) - Curtailed(1:744)).*(1-fract_solar(1:744)))/sum(power_orig(1:744).*(1-fract_solar(1:744)));
            W_feb_used(k+1,1) = sum((power_orig(745:1416) - Curtailed(745:1416)).*(1-fract_solar(745:1416)))/sum(power_orig(745:1416).*(1-fract_solar(745:1416)));
            W_mar_used(k+1,1) = sum((power_orig(1417:2160) - Curtailed(1417:2160)).*(1-fract_solar(1417:2160)))/sum(power_orig(1417:2160).*(1-fract_solar(1417:2160)));
            W_apr_used(k+1,1) = sum((power_orig(2161:2880) - Curtailed(2161:2880)).*(1-fract_solar(2161:2880)))/sum(power_orig(2161:2880).*(1-fract_solar(2161:2880)));
            W_may_used(k+1,1) = sum((power_orig(2881:3624) - Curtailed(2881:3624)).*(1-fract_solar(2881:3624)))/sum(power_orig(2881:3624).*(1-fract_solar(2881:3624)));
            W_jun_used(k+1,1) = sum((power_orig(3625:4344) - Curtailed(3625:4344)).*(1-fract_solar(3625:4344)))/sum(power_orig(3625:4344).*(1-fract_solar(3625:4344)));
            W_jul_used(k+1,1) = sum((power_orig(4345:5088) - Curtailed(4345:5088)).*(1-fract_solar(4345:5088)))/sum(power_orig(4345:5088).*(1-fract_solar(4345:5088)));
            W_aug_used(k+1,1) = sum((power_orig(5089:5832) - Curtailed(5089:5832)).*(1-fract_solar(5089:5832)))/sum(power_orig(5089:5832).*(1-fract_solar(5089:5832)));
            W_sep_used(k+1,1) = sum((power_orig(5833:6552) - Curtailed(5833:6552)).*(1-fract_solar(5833:6552)))/sum(power_orig(5833:6552).*(1-fract_solar(5833:6552)));
            W_oct_used(k+1,1) = sum((power_orig(6553:7296) - Curtailed(6553:7296)).*(1-fract_solar(6553:7296)))/sum(power_orig(6553:7296).*(1-fract_solar(6553:7296)));
            W_nov_used(k+1,1) = sum((power_orig(7297:8016) - Curtailed(7297:8016)).*(1-fract_solar(7297:8016)))/sum(power_orig(7297:8016).*(1-fract_solar(7297:8016)));
            W_dec_used(k+1,1) = sum((power_orig(8017:8760) - Curtailed(8017:8760)).*(1-fract_solar(8017:8760)))/sum(power_orig(8017:8760).*(1-fract_solar(8017:8760)));
        end
        % create matrix of percent of energy used each month for each of
        % the total percent curtailments
        S_month_used = [S_jan_used,S_feb_used,S_mar_used,S_apr_used,S_may_used,S_jun_used,S_jul_used,S_aug_used,S_sep_used,S_oct_used,S_nov_used,S_dec_used];
        W_month_used = [W_jan_used,W_feb_used,W_mar_used,W_apr_used,W_may_used,W_jun_used,W_jul_used,W_aug_used,W_sep_used,W_oct_used,W_nov_used,W_dec_used];
        %difference between the curtailment in the month at a given curtailment level and the previous level
        % add the final difference in energy used for the case case of 0% curtailment
        S_usedDiff = [S_month_used(2:end,:) - S_month_used(1:end-1,:)]; %difference between the curtailment in the month at a given curtailment level and the previous level
        S_usedDiff = [S_usedDiff ; 1-sum(S_usedDiff)];
        W_usedDiff = [W_month_used(2:end,:) - W_month_used(1:end-1,:)]; %difference between the curtailment in the month at a given curtailment level and the previous level
        W_usedDiff = [W_usedDiff ; 1-sum(W_usedDiff)];
        % remove entries that at NaN due to there being no solar or wind
        % installed
        W_usedDiff(isnan(W_usedDiff)) = 0;
        S_usedDiff(isnan(S_usedDiff)) = 0;

        costStep = CP_costs.topsolar_yearlyLCOU{i}(:,y); %cost of utilizing energy at each percent curtailment level
        LAEC = CP_costs.topsolar_yearlyLAEC{i}(:,y); % total energy cost (except the energy supply cost) at each percent curtailment level
        S_E_cost = CP_costs.topsolar_solaronly_Ecost(i,y); % solar energy supply cost in $/MWh
        W_E_cost = CP_costs.topsolar_windonly_Ecost(i,y); % wind energy supply cost in $/MWh
        % total value of energy in the month is the sum of the difference in curtailment multiplied by the value of the energy wasted
        monthlyValueTotal = sum(S_usedDiff .* (market_value - costStep)) + sum(W_usedDiff .* (market_value - costStep)); %value of energy in the month is the sum of the difference in curtailment multiplied by the value of the energy wasted
        % value of the used energy in the monthly is value of the energy
        % used multiplied by the change in the amount of energy for each
        % percent curtailment, only up to the optimial percent curtailment,
        % normalized by the total energy used in the month
        S_monthlyValueUsed = sum(S_usedDiff(1:(100-opt_curtailment(y)),:) .* (market_value - costStep(1:(100-opt_curtailment(y)),:))) ./ sum(S_usedDiff(1:(100-opt_curtailment(y)),:));
        % similarly for the value of curtailed energy, only analysing the
        % percent curtailment between the optimal energy utilization and
        % 100% energy utilization
        S_monthlyValueWasted = sum(S_usedDiff((101-opt_curtailment(y)):end,:) .* (market_value - costStep((101-opt_curtailment(y)):end,:))) ./ sum(S_usedDiff((101-opt_curtailment(y)):end,:));
        % total cost of the energy used, including the energy supply cost
        % and the cost of utilizing energy
        S_monthlyCostUsed = sum(S_usedDiff(1:(100-opt_curtailment(y)),:) .* ((S_E_cost*(100/(100-S_curt))) + LAEC(1:(100-opt_curtailment(y)),:))) ./ sum(S_usedDiff(1:(100-opt_curtailment(y)),:));
        % total fraction of energy used in each month
        S_monthlyXUsed = sum(S_usedDiff(1:(100-opt_curtailment(y)),:));
        % repeat for the wind energy 
        W_monthlyValueUsed = sum(W_usedDiff(1:(100-opt_curtailment(y)),:) .* (market_value - costStep(1:(100-opt_curtailment(y)),:))) ./ sum(W_usedDiff(1:(100-opt_curtailment(y)),:));
        W_monthlyValueWasted = sum(W_usedDiff((101-opt_curtailment(y)):end,:) .* (market_value - costStep((101-opt_curtailment(y)):end,:))) ./ sum(W_usedDiff((101-opt_curtailment(y)):end,:));
        W_monthlyCostUsed = sum(W_usedDiff(1:(100-opt_curtailment(y)),:) .* ((W_E_cost*(100/(100-W_curt))) + LAEC(1:(100-opt_curtailment(y)),:))) ./ sum(W_usedDiff(1:(100-opt_curtailment(y)),:));
        W_monthlyXUsed = sum(W_usedDiff(1:(100-opt_curtailment(y)),:));
        % remove the NaN entries
        S_monthlyValueUsed(isnan(S_monthlyValueUsed)) = 0;
        S_monthlyValueWasted(isnan(S_monthlyValueWasted)) = 0;
        S_monthlyCostUsed(isnan(S_monthlyCostUsed)) = 0;
        W_monthlyValueUsed(isnan(W_monthlyValueUsed)) = 0;
        W_monthlyValueWasted(isnan(W_monthlyValueWasted)) = 0;
        W_monthlyCostUsed(isnan(W_monthlyCostUsed)) = 0;
        %save the above metrics in the monthlyResults table
        monthlyResults{y,2:end} = monthlyValueTotal; %add results for that year to the record table
        monthlyResults{y+7,2:end} = Solar_E;
        monthlyResults{y+14,2:end} = S_monthlyValueUsed;
        monthlyResults{y+21,2:end} = S_monthlyValueWasted;
        monthlyResults{y+28,2:end} = S_monthlyCostUsed;
        monthlyResults{y+35,2:end} = S_monthlyXUsed;
        monthlyResults{y+42,2:end} = Wind_E;
        monthlyResults{y+49,2:end} = W_monthlyValueUsed;
        monthlyResults{y+56,2:end} = W_monthlyValueWasted;
        monthlyResults{y+63,2:end} = W_monthlyCostUsed;
        monthlyResults{y+70,2:end} = W_monthlyXUsed;
    end
    % find the average of the metrics over the 6 years, removing the zero
    % entries (i.e. only consider the year as counting toward the average
    % if solar or wind were installed for that year)
    for m = 2:13
        monthlyResults{7,m} = mean(nonzeros(monthlyResults{1:6,m}));
        monthlyResults{14,m} = mean(nonzeros(monthlyResults{8:13,m}));
        monthlyResults{21,m} = mean(nonzeros(monthlyResults{15:20,m}));
        monthlyResults{28,m} = mean(nonzeros(monthlyResults{22:27,m}));
        monthlyResults{35,m} = mean(nonzeros(monthlyResults{29:34,m}));
        monthlyResults{42,m} = mean(nonzeros(monthlyResults{36:41,m}));
        monthlyResults{49,m} = mean(nonzeros(monthlyResults{43:48,m}));
        monthlyResults{56,m} = mean(nonzeros(monthlyResults{50:55,m}));
        monthlyResults{63,m} = mean(nonzeros(monthlyResults{57:62,m}));
        monthlyResults{70,m} = mean(nonzeros(monthlyResults{64:69,m}));
        monthlyResults{77,m} = mean(nonzeros(monthlyResults{71:76,m}));
    end
    CP_costs.topsolar_monthlyValues{i} = monthlyResults;
end
save('CP_costs','CP_costs')
clearvars -except CP_costs CP_base_solar CP_base_wind CP_nocurt CP_curt

% 4.iii Calculate the amount of solar and wind utilization in each month for
% combined solar and wind at top wind locations
S = linspace(1,100)'; %list is a scaling array for percent curtailment. 
% It is used to resize all of the process units for each level of curtailment 
% because the optimization assumes constant total production rather than 
% constant solar panel and wind turbine size
market_value = 1000 * (17/1000000)*(3600000/640); %$ per MWh. given a price of green ammonia as $1000/tonne
PReq = 100; %nomial size 100 MW
results = CP_costs.topwind_results;
for i = 1:10 % cycle through all 10 locations
    %load the process design data from the optimization of 0-99%
    %curtailment
    panels = flip(results.panel{i},1);
    turbines = flip(results.turbine{i},1);
    H2Sizes = flip(results.H2Size{i},1); % GW
    BattSizes = flip(results.BattSize{i},1); % GW
    H2buffers = flip(results.H2buffer{i},1) * (0.2045) * 0.000278; % GWh
    batteries = flip(results.batterykWh{i},1) * (1/1E6); % GWh
    battery_uses = (results.batteryUse{i}); % GW
    % extract the panel power and turbine power and the optimal amount of curtailment at the
    % location and year
    loc = CP_costs.topwind_loc(i,1);
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc;
    panelPowers = CP_base_solar.panelPower(ind);
    turbinePowers = CP_base_wind.turbinePower(ind);
    opt_curtailment = round(CP_curt.curtailment(ind)*100);
    
    % initialize a table of monthly results for each of the six years.
    % Results are organized for solar and wind separately. 
    month = ["Value_Type","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
    varTypes = ["string","double","double","double","double","double","double","double","double","double","double","double","double"];
    monthlyResults = table('Size', [77 numel(month)], 'VariableTypes',varTypes,'VariableNames', month);
    monthlyResults.Value_Type = [
        "Year 1 Total";"Year 2 Total";"Year 3 Total";"Year 4 Total";"Year 5 Total";"Year 6 Total";"Average Total"; ...
        "Year 1 Solar Energy";"Year 2 Solar Energy";"Year 3 Solar Energy";"Year 4 Solar Energy";"Year 5 Solar Energy";"Year 6 Solar Energy";"Average Solar Energy"; ...
        "Solar Year 1 Used";"Solar Year 2 Used";"Solar Year 3 Used";"Solar Year 4 Used";"Solar Year 5 Used";"Solar Year 6 Used";"Solar Average Used";...
        "Solar Year 1 Wasted";"Solar Year 2 Wasted";"Solar Year 3 Wasted";"Solar Year 4 Wasted";"Solar Year 5 Wasted";"Solar Year 6 Wasted";"Solar Average Wasted";
        "Solar Year 1 LAEC";"Solar Year 2 LAEC";"Solar Year 3 LAEC";"Solar Year 4 LAEC";"Solar Year 5 LAEC";"Solar Year 6 LAEC";"Solar Average LAEC"; ...
        "Solar Year 1 Xused";"Solar Year 2 Xused";"Solar Year 3 Xused";"Solar Year 4 Xused";"Solar Year 5 Xused";"Solar Year 6 Xused";"Solar Average Xused"; ...
        "Year 1 Wind Energy";"Year 2 Wind Energy";"Year 3 Wind Energy";"Year 4 Wind Energy";"Year 5 Wind Energy";"Year 6 Wind Energy";"Average Wind Energy"; ...
        "Wind Year 1 Used";"Wind Year 2 Used";"Wind Year 3 Used";"Wind Year 4 Used";"Wind Year 5 Used";"Wind Year 6 Used";"Wind Average Used";...
        "Wind Year 1 Wasted";"Wind Year 2 Wasted";"Wind Year 3 Wasted";"Wind Year 4 Wasted";"Wind Year 5 Wasted";"Wind Year 6 Wasted";"Wind Average Wasted";
        "Wind Year 1 LAEC";"Wind Year 2 LAEC";"Wind Year 3 LAEC";"Wind Year 4 LAEC";"Wind Year 5 LAEC";"Wind Year 6 LAEC";"Wind Average LAEC"; ...
        "Wind Year 1 Xused";"Wind Year 2 Xused";"Wind Year 3 Xused";"Wind Year 4 Xused";"Wind Year 5 Xused";"Wind Year 6 Xused";"Wind Average Xused"];

    for y = 1:6 %cycle through each of the years
        battery_uses_year = flip(battery_uses{1,y},1); %extract and flip the optimization result for the battery usage over a year (includes all 6 years at a location)
        %initialize arrays for the usage of energy in each monthly for each of the optimizations between 0-99% curtailment 
        % arrays are for solar and wind separately
        S_jan_used = zeros(100,1);
        S_feb_used = zeros(100,1);
        S_mar_used = zeros(100,1);
        S_apr_used = zeros(100,1);
        S_may_used = zeros(100,1);
        S_jun_used = zeros(100,1);
        S_jul_used = zeros(100,1);
        S_aug_used = zeros(100,1);
        S_sep_used = zeros(100,1);
        S_oct_used = zeros(100,1);
        S_nov_used = zeros(100,1);
        S_dec_used = zeros(100,1);
        W_jan_used = zeros(100,1);
        W_feb_used = zeros(100,1);
        W_mar_used = zeros(100,1);
        W_apr_used = zeros(100,1);
        W_may_used = zeros(100,1);
        W_jun_used = zeros(100,1);
        W_jul_used = zeros(100,1);
        W_aug_used = zeros(100,1);
        W_sep_used = zeros(100,1);
        W_oct_used = zeros(100,1);
        W_nov_used = zeros(100,1);
        W_dec_used = zeros(100,1);
         
        % Start of sub-section to find the optimal amount of solar and wind
        % curtailment respectively at the optimal total curtailment
        k = 100 - (opt_curtailment(y));
        buffCap = H2buffers(k,y) * S(k)/100; %set the buffer capacity based on the optimization result scale by the percent curtailment
        Curtailed = zeros(8760,1); % initialize array of curtailment for each hour over a year in GWh
        demand = single((S(k)/100) * 0.1 * (612.5/640) * ones(8760,1)); %set the hourly demand matrix in GW, scaled according to the percent curtailment
        power = single(panels(end,y) * panelPowers{y} * (1/1E9) + turbines(end,y) * turbinePowers{y} * (1/1E9)); %initialize the hourly power profile
        fract_solar = (panels(end,y) * panelPowers{y})./(panels(end,y)*panelPowers{y} + turbines(end,y)*turbinePowers{y});
        fract_solar(isnan(fract_solar)) = 0;
        Solar_E = [sum(power(1:744).*fract_solar(1:744)),...
            sum(power(745:1416).*fract_solar(745:1416)),...
            sum(power(1417:2160).*fract_solar(1417:2160)),...
            sum(power(2161:2880).*fract_solar(2161:2880)),...
            sum(power(2881:3624).*fract_solar(2881:3624)),...
            sum(power(3625:4344).*fract_solar(3625:4344)),...
            sum(power(4345:5088).*fract_solar(4345:5088)),...
            sum(power(5089:5832).*fract_solar(5089:5832)),...
            sum(power(5833:6552).*fract_solar(5833:6552)),...
            sum(power(6553:7296).*fract_solar(6553:7296)),...
            sum(power(7297:8016).*fract_solar(7297:8016)),...
            sum(power(8017:8760).*fract_solar(8017:8760))];
        Wind_E = [sum(power(1:744).*(1-fract_solar(1:744))),...
            sum(power(745:1416).*(1-fract_solar(745:1416))),...
            sum(power(1417:2160).*(1-fract_solar(1417:2160))),...
            sum(power(2161:2880).*(1-fract_solar(2161:2880))),...
            sum(power(2881:3624).*(1-fract_solar(2881:3624))),...
            sum(power(3625:4344).*(1-fract_solar(3625:4344))),...
            sum(power(4345:5088).*(1-fract_solar(4345:5088))),...
            sum(power(5089:5832).*(1-fract_solar(5089:5832))),...
            sum(power(5833:6552).*(1-fract_solar(5833:6552))),...
            sum(power(6553:7296).*(1-fract_solar(6553:7296))),...
            sum(power(7297:8016).*(1-fract_solar(7297:8016))),...
            sum(power(8017:8760).*(1-fract_solar(8017:8760)))];

        power_orig = power; %set the original power profile to use later
        battery_use = battery_uses_year{k,1}; %extract the particular hourly battery use over a year for a given location and year
        battery_use = (battery_use + (27.5/640)*PReq*(1/1000)) * S(k)/100; %baseline electrical requirement added to battery usage, making it electricity extraction, but keep name as battery
        battery_use(battery_use < 0 ) = 0; %if battery use is negative, set it to zero.
        power = power - battery_use; %power supply after removing the electricity extraction
        elect_curtail = power - H2Sizes(k,y)*S(k)/100; %calculate the amount of electricity curtailed due to it being above the rated capacity of the electrolyser
        elect_curtail(elect_curtail < 0) = 0; %remove entries less than zeros (keep only those where electricity must be curtailed)
        Curtailed = Curtailed + elect_curtail; %add the curtailed electricity to the curtailment array
        power(power > (H2Sizes(k,y)*S(k)/100)) = H2Sizes(k,y)*S(k)/100; %adjust the power supply according to the maximum power set by the electrolyser capacity
        cum_energy = [0;cumsum(power)]; % create a vector of cumulative energy supply for each hour over a year
        cum_demand = [0;cumsum(demand)]; % create a vector of the cumulative energy demand (as hydrogen) for each hour over a year
        Current_def = -1 * (cum_energy(1) - cum_energy(end)) + (cum_demand(1) - cum_demand(end)); % calculate the initial overall energy balance (should be positive)
        
        active_vector = ones(8760,1); %create a vector which describes which hours during the year are "active", meaning more energy during the hour can still be curtailed
        dummy_power = power; %initialize dummy array of power supply
        cum_dummy_power = [0;cumsum(dummy_power)]; % calculate cumulative dummy power

        while Current_def > 0.001 % continue to curtail energy until the overall energy balance over a year is clost to zero
            Avail_power = dummy_power .* active_vector; % calculate an array of the available power
            cum_Avail_power = [0;cumsum(Avail_power)]; % array of cumulative available power
            % create a matrix of the maxaimum percent curtailment which
            % can occur between the starting hour (row entry) and the
            % ending hour (column entry) before the energy deficit over
            % that period is equal to the hydrogen storage capacity
            max_curtail = ((-1 * triu(cum_dummy_power - cum_dummy_power') + triu(cum_demand - cum_demand')) + buffCap*ones(8761,8761))./(-1 * triu(cum_Avail_power - cum_Avail_power'));
                
            max_curtail(max_curtail<=0) = NaN; % remove entries less than zero
            curt = min(min(max_curtail)); % find the overall maximum amount of curtailment before the hydrogen storage capacity is reached
            % if the maximum curtailment is less than the amount of
            % curtailment which will cause the energy balance over the
            % year to close, record the starting and ending indexes.
            % Else, set the percent curtailment the amount which will
            % close the energy balance.    
            if curt < (Current_def/(-1*cum_Avail_power(1) + cum_Avail_power(end)))
                [r,c] = find(max_curtail == curt);
            else
                curt = Current_def/(-1+cum_Avail_power(1) + cum_Avail_power(end));
                r = [0];
                c = [0];
            end
            % add the hourly energy curtailed in this iteration to the
            % array of total curtailed energy
            Curtailed = Curtailed + Avail_power * curt;
            dummy_power = dummy_power - (Avail_power * curt); %update the dummy power array after removing curtailed energy
            cum_dummy_power = [0;cumsum(dummy_power)]; % update cumulative power array
            % if the percent curtailment is not set by closing the
            % yearly energy balance, then remove the active indexes
            % which reached their maximum level of curtailment (the
            % deficit over that period is equal to the hydrogen storage
            % capacity)   
            if r(1,1)~=0
                for n = 1:height(r)
                    active_vector(r(n,1):c(n,1)-1,1) = zeros((c(n,1)-r(n,1)),1);
                end
            end
            %update the overall energy balance over a year
            Current_def = -1 * (cum_dummy_power(1) - cum_dummy_power(end)) + (cum_demand(1) - cum_demand(end));
        end
        S_curt = round(100 * sum(Curtailed .* fract_solar)/sum(power_orig.*fract_solar));
        W_curt = round(100 * sum(Curtailed .* (1-fract_solar))/sum(power_orig.*(1-fract_solar)));
        % End of sub-section to find the optimal amount of solar and wind
        % curtailment respectively at the optimal total curtailment
        
        k = 1;
        parfor k = 1:99 %cycle through each of 1-99% curtailment
            disp(k)
            buffCap = H2buffers(k,y) * S(k)/100; %set the buffer capacity based on the optimization result scale by the percent curtailment
            Curtailed = zeros(8760,1); % initialize array of curtailment for each hour over a year in GWh
            demand = single((S(k)/100) * 0.1 * (612.5/640) * ones(8760,1)); %set the hourly demand matrix in GW, scaled according to the percent curtailment
            power = single(panels(end,y) * panelPowers{y} * (1/1E9) + turbines(end,y) * turbinePowers{y} * (1/1E9)); %initialize the hourly power profile
            fract_solar = (panels(end,y) * panelPowers{y})./(panels(end,y)*panelPowers{y} + turbines(end,y)*turbinePowers{y})
            fract_solar(isnan(fract_solar)) = 0;

            power_orig = power; %set the original power profile to use later
            battery_use = battery_uses_year{k,1}; %extract the particular hourly battery use over a year for a given location and year
            battery_use = (battery_use + (27.5/640)*PReq*(1/1000)) * S(k)/100; %baseline electrical requirement added to battery usage, making is electricity extraction, but keep name as battery
            battery_use(battery_use < 0 ) = 0; %if battery use is negative, set it to zero.
            power = power - battery_use; %power supply after removing the electricity extraction
            elect_curtail = power - H2Sizes(k,y)*S(k)/100; %calculate the amount of electricity curtailed due to it being above the rated capacity of the electrolyser
            elect_curtail(elect_curtail < 0) = 0; %remove entries less than zeros (keep only those where electricity must be curtailed)
            Curtailed = Curtailed + elect_curtail; %add the curtailed electricity to the curtailment array
            power(power > (H2Sizes(k,y)*S(k)/100)) = H2Sizes(k,y)*S(k)/100; %adjust the power supply according to the maximum power set by the electrolyser capacity
            cum_energy = [0;cumsum(power)]; % create a vector of cumulative energy supply for each hour over a year
            cum_demand = [0;cumsum(demand)]; % create a vector of the cumulative energy demand (as hydrogen) for each hour over a year
            Current_def = -1 * (cum_energy(1) - cum_energy(end)) + (cum_demand(1) - cum_demand(end)); % calculate the initial overall energy balance (should be positive)
            
            active_vector = ones(8760,1); %create a vector which describes which hours during the year are "active", meaning more energy during the hour can still be curtailed
            dummy_power = power; %initialize dummy array of power supply
            cum_dummy_power = [0;cumsum(dummy_power)]; % calculate cumulative dummy power

            while Current_def > 0.001 % continue to curtail energy until the overall energy balance over a year is close to zero
                Avail_power = dummy_power .* active_vector; % calculate an array of the available power
                cum_Avail_power = [0;cumsum(Avail_power)]; % array of cumulative available power
                % create a matrix of the maxaimum percent curtailment which
                % can occur between the starting hour (row entry) and the
                % ending hour (column entry) before the energy deficit over
                % that period is equal to the hydrogen storage capacity
                max_curtail = ((-1 * triu(cum_dummy_power - cum_dummy_power') + triu(cum_demand - cum_demand')) + buffCap*ones(8761,8761))./(-1 * triu(cum_Avail_power - cum_Avail_power'));
                
                max_curtail(max_curtail<=0) = NaN; % remove entries less than zero
                curt = min(min(max_curtail)); % find the overall maximum amount of curtailment before the hydrogen storage capacity is reached
                % if the maximum curtailment is less than the amount of
                % curtailment which will cause the energy balance over the
                % year to close, record the starting and ending indexes.
                % Else, set the percent curtailment the amount which will
                % close the energy balance. 
                if curt < (Current_def/(-1*cum_Avail_power(1) + cum_Avail_power(end)))
                    [r,c] = find(max_curtail == curt);
                else
                    curt = Current_def/(-1+cum_Avail_power(1) + cum_Avail_power(end));
                    r = [0];
                    c = [0];
                end
                % add the hourly energy curtailed in this iteration to the
                % array of total curtailed energy
                Curtailed = Curtailed + Avail_power * curt;
                dummy_power = dummy_power - (Avail_power * curt); %update the dummy power array after removing curtailed energy
                cum_dummy_power = [0;cumsum(dummy_power)]; % update cumulative power array
                % if the percent curtailment is not set by closing the
                % yearly energy balance, then remove the active indexes
                % which reached their maximum level of curtailment (the
                % deficit over that period is equal to the hydrogen storage
                % capacity)
                if r(1,1)~=0
                    for n = 1:height(r)
                        active_vector(r(n,1):c(n,1)-1,1) = zeros((c(n,1)-r(n,1)),1);
                    end
                end
                %update the overall energy balance over a year
                Current_def = -1 * (cum_dummy_power(1) - cum_dummy_power(end)) + (cum_demand(1) - cum_demand(end));
            end
            % save the percent of energy used in each month for the current
            % iteration of total percent curtailment
            S_jan_used(k+1,1) = sum((power_orig(1:744) - Curtailed(1:744)).*fract_solar(1:744))/sum(power_orig(1:744).*fract_solar(1:744));
            S_feb_used(k+1,1) = sum((power_orig(745:1416) - Curtailed(745:1416)).*fract_solar(745:1416))/sum(power_orig(745:1416).*fract_solar(745:1416));
            S_mar_used(k+1,1) = sum((power_orig(1417:2160) - Curtailed(1417:2160)).*fract_solar(1417:2160))/sum(power_orig(1417:2160).*fract_solar(1417:2160));
            S_apr_used(k+1,1) = sum((power_orig(2161:2880) - Curtailed(2161:2880)).*fract_solar(2161:2880))/sum(power_orig(2161:2880).*fract_solar(2161:2880));
            S_may_used(k+1,1) = sum((power_orig(2881:3624) - Curtailed(2881:3624)).*fract_solar(2881:3624))/sum(power_orig(2881:3624).*fract_solar(2881:3624));
            S_jun_used(k+1,1) = sum((power_orig(3625:4344) - Curtailed(3625:4344)).*fract_solar(3625:4344))/sum(power_orig(3625:4344).*fract_solar(3625:4344));
            S_jul_used(k+1,1) = sum((power_orig(4345:5088) - Curtailed(4345:5088)).*fract_solar(4345:5088))/sum(power_orig(4345:5088).*fract_solar(4345:5088));
            S_aug_used(k+1,1) = sum((power_orig(5089:5832) - Curtailed(5089:5832)).*fract_solar(5089:5832))/sum(power_orig(5089:5832).*fract_solar(5089:5832));
            S_sep_used(k+1,1) = sum((power_orig(5833:6552) - Curtailed(5833:6552)).*fract_solar(5833:6552))/sum(power_orig(5833:6552).*fract_solar(5833:6552));
            S_oct_used(k+1,1) = sum((power_orig(6553:7296) - Curtailed(6553:7296)).*fract_solar(6553:7296))/sum(power_orig(6553:7296).*fract_solar(6553:7296));
            S_nov_used(k+1,1) = sum((power_orig(7297:8016) - Curtailed(7297:8016)).*fract_solar(7297:8016))/sum(power_orig(7297:8016).*fract_solar(7297:8016));
            S_dec_used(k+1,1) = sum((power_orig(8017:8760) - Curtailed(8017:8760)).*fract_solar(8017:8760))/sum(power_orig(8017:8760).*fract_solar(8017:8760));
            
            W_jan_used(k+1,1) = sum((power_orig(1:744) - Curtailed(1:744)).*(1-fract_solar(1:744)))/sum(power_orig(1:744).*(1-fract_solar(1:744)));
            W_feb_used(k+1,1) = sum((power_orig(745:1416) - Curtailed(745:1416)).*(1-fract_solar(745:1416)))/sum(power_orig(745:1416).*(1-fract_solar(745:1416)));
            W_mar_used(k+1,1) = sum((power_orig(1417:2160) - Curtailed(1417:2160)).*(1-fract_solar(1417:2160)))/sum(power_orig(1417:2160).*(1-fract_solar(1417:2160)));
            W_apr_used(k+1,1) = sum((power_orig(2161:2880) - Curtailed(2161:2880)).*(1-fract_solar(2161:2880)))/sum(power_orig(2161:2880).*(1-fract_solar(2161:2880)));
            W_may_used(k+1,1) = sum((power_orig(2881:3624) - Curtailed(2881:3624)).*(1-fract_solar(2881:3624)))/sum(power_orig(2881:3624).*(1-fract_solar(2881:3624)));
            W_jun_used(k+1,1) = sum((power_orig(3625:4344) - Curtailed(3625:4344)).*(1-fract_solar(3625:4344)))/sum(power_orig(3625:4344).*(1-fract_solar(3625:4344)));
            W_jul_used(k+1,1) = sum((power_orig(4345:5088) - Curtailed(4345:5088)).*(1-fract_solar(4345:5088)))/sum(power_orig(4345:5088).*(1-fract_solar(4345:5088)));
            W_aug_used(k+1,1) = sum((power_orig(5089:5832) - Curtailed(5089:5832)).*(1-fract_solar(5089:5832)))/sum(power_orig(5089:5832).*(1-fract_solar(5089:5832)));
            W_sep_used(k+1,1) = sum((power_orig(5833:6552) - Curtailed(5833:6552)).*(1-fract_solar(5833:6552)))/sum(power_orig(5833:6552).*(1-fract_solar(5833:6552)));
            W_oct_used(k+1,1) = sum((power_orig(6553:7296) - Curtailed(6553:7296)).*(1-fract_solar(6553:7296)))/sum(power_orig(6553:7296).*(1-fract_solar(6553:7296)));
            W_nov_used(k+1,1) = sum((power_orig(7297:8016) - Curtailed(7297:8016)).*(1-fract_solar(7297:8016)))/sum(power_orig(7297:8016).*(1-fract_solar(7297:8016)));
            W_dec_used(k+1,1) = sum((power_orig(8017:8760) - Curtailed(8017:8760)).*(1-fract_solar(8017:8760)))/sum(power_orig(8017:8760).*(1-fract_solar(8017:8760)));
        end
        % create matrix of percent of energy used each month for each of
        % the total percent curtailments
        S_month_used = [S_jan_used,S_feb_used,S_mar_used,S_apr_used,S_may_used,S_jun_used,S_jul_used,S_aug_used,S_sep_used,S_oct_used,S_nov_used,S_dec_used];
        W_month_used = [W_jan_used,W_feb_used,W_mar_used,W_apr_used,W_may_used,W_jun_used,W_jul_used,W_aug_used,W_sep_used,W_oct_used,W_nov_used,W_dec_used];
        %difference between the curtailment in the month at a given curtailment level and the previous level
        % add the final difference in energy used for the case case of 0% curtailment
        S_usedDiff = [S_month_used(2:end,:) - S_month_used(1:end-1,:)]; %difference between the curtailment in the month at a given curtailment level and the previous level
        S_usedDiff = [S_usedDiff ; 1-sum(S_usedDiff)];
        W_usedDiff = [W_month_used(2:end,:) - W_month_used(1:end-1,:)]; %difference between the curtailment in the month at a given curtailment level and the previous level
        W_usedDiff = [W_usedDiff ; 1-sum(W_usedDiff)];
        % remove entries that at NaN due to there being no solar or wind
        % installed
        W_usedDiff(isnan(W_usedDiff)) = 0;
        S_usedDiff(isnan(S_usedDiff)) = 0;

        costStep = CP_costs.topwind_yearlyLCOU{i}(:,y); %cost of utilizing energy at each percent curtailment level
        LAEC = CP_costs.topwind_yearlyLAEC{i}(:,y); % total energy cost (except the energy supply cost) at each percent curtailment level
        S_E_cost = CP_costs.topwind_solaronly_Ecost(i,y); % solar energy supply cost in $/MWh
        W_E_cost = CP_costs.topwind_windonly_Ecost(i,y); % wind energy supply cost in $/MWh
        % total value of energy in the month is the sum of the difference in curtailment multiplied by the value of the energy wasted
        monthlyValueTotal = sum(S_usedDiff .* (market_value - costStep)) + sum(W_usedDiff .* (market_value - costStep)); %value of energy in the month is the sum of the difference in curtailment multiplied by the value of the energy wasted
        % value of the used energy in the monthly is value of the energy
        % used multiplied by the change in the amount of energy for each
        % percent curtailment, only up to the optimial percent curtailment,
        % normalized by the total energy used in the month
        S_monthlyValueUsed = sum(S_usedDiff(1:(100-opt_curtailment(y)),:) .* (market_value - costStep(1:(100-opt_curtailment(y)),:))) ./ sum(S_usedDiff(1:(100-opt_curtailment(y)),:));
        % similarly for the value of curtailed energy, only analysing the
        % percent curtailment between the optimal energy utilization and
        % 100% energy utilization
        S_monthlyValueWasted = sum(S_usedDiff((101-opt_curtailment(y)):end,:) .* (market_value - costStep((101-opt_curtailment(y)):end,:))) ./ sum(S_usedDiff((101-opt_curtailment(y)):end,:));
        % total cost of the energy used, including the energy supply cost
        % and the cost of utilizing energy
        S_monthlyCostUsed = sum(S_usedDiff(1:(100-opt_curtailment(y)),:) .* ((S_E_cost*(100/(100-S_curt))) + LAEC(1:(100-opt_curtailment(y)),:))) ./ sum(S_usedDiff(1:(100-opt_curtailment(y)),:));
        % total fraction of energy used in each month
        S_monthlyXUsed = sum(S_usedDiff(1:(100-opt_curtailment(y)),:));
        % repeat for the wind energy 
        W_monthlyValueUsed = sum(W_usedDiff(1:(100-opt_curtailment(y)),:) .* (market_value - costStep(1:(100-opt_curtailment(y)),:))) ./ sum(W_usedDiff(1:(100-opt_curtailment(y)),:));
        W_monthlyValueWasted = sum(W_usedDiff((101-opt_curtailment(y)):end,:) .* (market_value - costStep((101-opt_curtailment(y)):end,:))) ./ sum(W_usedDiff((101-opt_curtailment(y)):end,:));
        W_monthlyCostUsed = sum(W_usedDiff(1:(100-opt_curtailment(y)),:) .* ((W_E_cost*(100/(100-W_curt))) + LAEC(1:(100-opt_curtailment(y)),:))) ./ sum(W_usedDiff(1:(100-opt_curtailment(y)),:));
        W_monthlyXUsed = sum(W_usedDiff(1:(100-opt_curtailment(y)),:));
        % remove the NaN entries
        S_monthlyValueUsed(isnan(S_monthlyValueUsed)) = 0;
        S_monthlyValueWasted(isnan(S_monthlyValueWasted)) = 0;
        S_monthlyCostUsed(isnan(S_monthlyCostUsed)) = 0;
        W_monthlyValueUsed(isnan(W_monthlyValueUsed)) = 0;
        W_monthlyValueWasted(isnan(W_monthlyValueWasted)) = 0;
        W_monthlyCostUsed(isnan(W_monthlyCostUsed)) = 0;
        %save the above metrics in the monthlyResults table
        monthlyResults{y,2:end} = monthlyValueTotal; %add results for that year to the record table
        monthlyResults{y+7,2:end} = Solar_E;
        monthlyResults{y+14,2:end} = S_monthlyValueUsed;
        monthlyResults{y+21,2:end} = S_monthlyValueWasted;
        monthlyResults{y+28,2:end} = S_monthlyCostUsed;
        monthlyResults{y+35,2:end} = S_monthlyXUsed;
        monthlyResults{y+42,2:end} = Wind_E;
        monthlyResults{y+49,2:end} = W_monthlyValueUsed;
        monthlyResults{y+56,2:end} = W_monthlyValueWasted;
        monthlyResults{y+63,2:end} = W_monthlyCostUsed;
        monthlyResults{y+70,2:end} = W_monthlyXUsed;
    end
    % find the average of the metrics over the 6 years, removing the zero
    % entries (i.e. only consider the year as counting toward the average
    % if solar or wind were installed for that year)
    for m = 2:13
        monthlyResults{7,m} = mean(nonzeros(monthlyResults{1:6,m}));
        monthlyResults{14,m} = mean(nonzeros(monthlyResults{8:13,m}));
        monthlyResults{21,m} = mean(nonzeros(monthlyResults{15:20,m}));
        monthlyResults{28,m} = mean(nonzeros(monthlyResults{22:27,m}));
        monthlyResults{35,m} = mean(nonzeros(monthlyResults{29:34,m}));
        monthlyResults{42,m} = mean(nonzeros(monthlyResults{36:41,m}));
        monthlyResults{49,m} = mean(nonzeros(monthlyResults{43:48,m}));
        monthlyResults{56,m} = mean(nonzeros(monthlyResults{50:55,m}));
        monthlyResults{63,m} = mean(nonzeros(monthlyResults{57:62,m}));
        monthlyResults{70,m} = mean(nonzeros(monthlyResults{64:69,m}));
        monthlyResults{77,m} = mean(nonzeros(monthlyResults{71:76,m}));
    end
    CP_costs.topwind_monthlyValues{i} = monthlyResults;
end
save('CP_costs','CP_costs')
clearvars -except CP_costs CP_base_solar CP_base_wind CP_nocurt CP_curt

% 4.iv Calculate the amount of solar and wind utilization in each month for
% combined solar and wind at top locations
S = linspace(1,100)'; %list is a scaling array for percent curtailment. 
% It is used to resize all of the process units for each level of curtailment 
% because the optimization assumes constant total production rather than 
% constant solar panel and wind turbine size
market_value = 1000 * (17/1000000)*(3600000/640); %$ per MWh. given a price of green ammonia as $1000/tonne
PReq = 100; %nomial size 100 MW
results = CP_costs.top_results;
for i = 1:10 % cycle through all 10 locations
    %load the process design data from the optimization of 0-99%
    %curtailment
    panels = flip(results.panel{i},1);
    turbines = flip(results.turbine{i},1);
    H2Sizes = flip(results.H2Size{i},1); % GW
    BattSizes = flip(results.BattSize{i},1); % GW
    H2buffers = flip(results.H2buffer{i},1) * (0.2045) * 0.000278; % GWh
    batteries = flip(results.batterykWh{i},1) * (1/1E6); % GWh
    battery_uses = (results.batteryUse{i}); % GW
    % extract the panel power and turbine power and the optimal amount of curtailment at the
    % location and year
    loc = CP_costs.top_loc(i,1);
    disp(loc)
    ind = CP_base_solar.location_SYNC == loc;
    panelPowers = CP_base_solar.panelPower(ind);
    turbinePowers = CP_base_wind.turbinePower(ind);
    opt_curtailment = round(CP_curt.curtailment(ind)*100);
    
    % initialize a table of monthly results for each of the six years.
    % Results are organized for solar and wind separately. 
    month = ["Value_Type","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
    varTypes = ["string","double","double","double","double","double","double","double","double","double","double","double","double"];
    monthlyResults = table('Size', [77 numel(month)], 'VariableTypes',varTypes,'VariableNames', month);
    monthlyResults.Value_Type = [
        "Year 1 Total";"Year 2 Total";"Year 3 Total";"Year 4 Total";"Year 5 Total";"Year 6 Total";"Average Total"; ...
        "Year 1 Solar Energy";"Year 2 Solar Energy";"Year 3 Solar Energy";"Year 4 Solar Energy";"Year 5 Solar Energy";"Year 6 Solar Energy";"Average Solar Energy"; ...
        "Solar Year 1 Used";"Solar Year 2 Used";"Solar Year 3 Used";"Solar Year 4 Used";"Solar Year 5 Used";"Solar Year 6 Used";"Solar Average Used";...
        "Solar Year 1 Wasted";"Solar Year 2 Wasted";"Solar Year 3 Wasted";"Solar Year 4 Wasted";"Solar Year 5 Wasted";"Solar Year 6 Wasted";"Solar Average Wasted";
        "Solar Year 1 LAEC";"Solar Year 2 LAEC";"Solar Year 3 LAEC";"Solar Year 4 LAEC";"Solar Year 5 LAEC";"Solar Year 6 LAEC";"Solar Average LAEC"; ...
        "Solar Year 1 Xused";"Solar Year 2 Xused";"Solar Year 3 Xused";"Solar Year 4 Xused";"Solar Year 5 Xused";"Solar Year 6 Xused";"Solar Average Xused"; ...
        "Year 1 Wind Energy";"Year 2 Wind Energy";"Year 3 Wind Energy";"Year 4 Wind Energy";"Year 5 Wind Energy";"Year 6 Wind Energy";"Average Wind Energy"; ...
        "Wind Year 1 Used";"Wind Year 2 Used";"Wind Year 3 Used";"Wind Year 4 Used";"Wind Year 5 Used";"Wind Year 6 Used";"Wind Average Used";...
        "Wind Year 1 Wasted";"Wind Year 2 Wasted";"Wind Year 3 Wasted";"Wind Year 4 Wasted";"Wind Year 5 Wasted";"Wind Year 6 Wasted";"Wind Average Wasted";
        "Wind Year 1 LAEC";"Wind Year 2 LAEC";"Wind Year 3 LAEC";"Wind Year 4 LAEC";"Wind Year 5 LAEC";"Wind Year 6 LAEC";"Wind Average LAEC"; ...
        "Wind Year 1 Xused";"Wind Year 2 Xused";"Wind Year 3 Xused";"Wind Year 4 Xused";"Wind Year 5 Xused";"Wind Year 6 Xused";"Wind Average Xused"];

    for y = 1:6 %cycle through each of the years
        battery_uses_year = flip(battery_uses{1,y},1); %extract and flip the optimization result for the battery usage over a year (includes all 6 years at a location)
        %initialize arrays for the usage of energy in each monthly for each of the optimizations between 0-99% curtailment 
        % arrays are for solar and wind separately
        S_jan_used = zeros(100,1);
        S_feb_used = zeros(100,1);
        S_mar_used = zeros(100,1);
        S_apr_used = zeros(100,1);
        S_may_used = zeros(100,1);
        S_jun_used = zeros(100,1);
        S_jul_used = zeros(100,1);
        S_aug_used = zeros(100,1);
        S_sep_used = zeros(100,1);
        S_oct_used = zeros(100,1);
        S_nov_used = zeros(100,1);
        S_dec_used = zeros(100,1);
        W_jan_used = zeros(100,1);
        W_feb_used = zeros(100,1);
        W_mar_used = zeros(100,1);
        W_apr_used = zeros(100,1);
        W_may_used = zeros(100,1);
        W_jun_used = zeros(100,1);
        W_jul_used = zeros(100,1);
        W_aug_used = zeros(100,1);
        W_sep_used = zeros(100,1);
        W_oct_used = zeros(100,1);
        W_nov_used = zeros(100,1);
        W_dec_used = zeros(100,1);
         
        % Start of sub-section to find the optimal amount of solar and wind
        % curtailment respectively at the optimal total curtailment
        k = 100 - (opt_curtailment(y));
        buffCap = H2buffers(k,y) * S(k)/100; %set the buffer capacity based on the optimization result scale by the percent curtailment
        Curtailed = zeros(8760,1); % initialize array of curtailment for each hour over a year in GWh
        demand = single((S(k)/100) * 0.1 * (612.5/640) * ones(8760,1)); %set the hourly demand matrix in GW, scaled according to the percent curtailment
        power = single(panels(end,y) * panelPowers{y} * (1/1E9) + turbines(end,y) * turbinePowers{y} * (1/1E9)); %initialize the hourly power profile
        fract_solar = (panels(end,y) * panelPowers{y})./(panels(end,y)*panelPowers{y} + turbines(end,y)*turbinePowers{y});
        fract_solar(isnan(fract_solar)) = 0;
        Solar_E = [sum(power(1:744).*fract_solar(1:744)),...
            sum(power(745:1416).*fract_solar(745:1416)),...
            sum(power(1417:2160).*fract_solar(1417:2160)),...
            sum(power(2161:2880).*fract_solar(2161:2880)),...
            sum(power(2881:3624).*fract_solar(2881:3624)),...
            sum(power(3625:4344).*fract_solar(3625:4344)),...
            sum(power(4345:5088).*fract_solar(4345:5088)),...
            sum(power(5089:5832).*fract_solar(5089:5832)),...
            sum(power(5833:6552).*fract_solar(5833:6552)),...
            sum(power(6553:7296).*fract_solar(6553:7296)),...
            sum(power(7297:8016).*fract_solar(7297:8016)),...
            sum(power(8017:8760).*fract_solar(8017:8760))];
        Wind_E = [sum(power(1:744).*(1-fract_solar(1:744))),...
            sum(power(745:1416).*(1-fract_solar(745:1416))),...
            sum(power(1417:2160).*(1-fract_solar(1417:2160))),...
            sum(power(2161:2880).*(1-fract_solar(2161:2880))),...
            sum(power(2881:3624).*(1-fract_solar(2881:3624))),...
            sum(power(3625:4344).*(1-fract_solar(3625:4344))),...
            sum(power(4345:5088).*(1-fract_solar(4345:5088))),...
            sum(power(5089:5832).*(1-fract_solar(5089:5832))),...
            sum(power(5833:6552).*(1-fract_solar(5833:6552))),...
            sum(power(6553:7296).*(1-fract_solar(6553:7296))),...
            sum(power(7297:8016).*(1-fract_solar(7297:8016))),...
            sum(power(8017:8760).*(1-fract_solar(8017:8760)))];

        power_orig = power; %set the original power profile to use later
        battery_use = battery_uses_year{k,1}; %extract the particular hourly battery use over a year for a given location and year
        battery_use = (battery_use + (27.5/640)*PReq*(1/1000)) * S(k)/100; %baseline electrical requirement added to battery usage, making it electricity extraction, but keep name as battery
        battery_use(battery_use < 0 ) = 0; %if battery use is negative, set it to zero.
        power = power - battery_use; %power supply after removing the electricity extraction
        elect_curtail = power - H2Sizes(k,y)*S(k)/100; %calculate the amount of electricity curtailed due to it being above the rated capacity of the electrolyser
        elect_curtail(elect_curtail < 0) = 0; %remove entries less than zeros (keep only those where electricity must be curtailed)
        Curtailed = Curtailed + elect_curtail; %add the curtailed electricity to the curtailment array
        power(power > (H2Sizes(k,y)*S(k)/100)) = H2Sizes(k,y)*S(k)/100; %adjust the power supply according to the maximum power set by the electrolyser capacity
        cum_energy = [0;cumsum(power)]; % create a vector of cumulative energy supply for each hour over a year
        cum_demand = [0;cumsum(demand)]; % create a vector of the cumulative energy demand (as hydrogen) for each hour over a year
        Current_def = -1 * (cum_energy(1) - cum_energy(end)) + (cum_demand(1) - cum_demand(end)); % calculate the initial overall energy balance (should be positive)
        
        active_vector = ones(8760,1); %create a vector which describes which hours during the year are "active", meaning more energy during the hour can still be curtailed
        dummy_power = power; %initialize dummy array of power supply
        cum_dummy_power = [0;cumsum(dummy_power)]; % calculate cumulative dummy power

        while Current_def > 0.001 % continue to curtail energy until the overall energy balance over a year is clost to zero
            Avail_power = dummy_power .* active_vector; % calculate an array of the available power
            cum_Avail_power = [0;cumsum(Avail_power)]; % array of cumulative available power
            % create a matrix of the maxaimum percent curtailment which
            % can occur between the starting hour (row entry) and the
            % ending hour (column entry) before the energy deficit over
            % that period is equal to the hydrogen storage capacity
            max_curtail = ((-1 * triu(cum_dummy_power - cum_dummy_power') + triu(cum_demand - cum_demand')) + buffCap*ones(8761,8761))./(-1 * triu(cum_Avail_power - cum_Avail_power'));
                
            max_curtail(max_curtail<=0) = NaN; % remove entries less than zero
            curt = min(min(max_curtail)); % find the overall maximum amount of curtailment before the hydrogen storage capacity is reached
            % if the maximum curtailment is less than the amount of
            % curtailment which will cause the energy balance over the
            % year to close, record the starting and ending indexes.
            % Else, set the percent curtailment the amount which will
            % close the energy balance.    
            if curt < (Current_def/(-1*cum_Avail_power(1) + cum_Avail_power(end)))
                [r,c] = find(max_curtail == curt);
            else
                curt = Current_def/(-1+cum_Avail_power(1) + cum_Avail_power(end));
                r = [0];
                c = [0];
            end
            % add the hourly energy curtailed in this iteration to the
            % array of total curtailed energy
            Curtailed = Curtailed + Avail_power * curt;
            dummy_power = dummy_power - (Avail_power * curt); %update the dummy power array after removing curtailed energy
            cum_dummy_power = [0;cumsum(dummy_power)]; % update cumulative power array
            % if the percent curtailment is not set by closing the
            % yearly energy balance, then remove the active indexes
            % which reached their maximum level of curtailment (the
            % deficit over that period is equal to the hydrogen storage
            % capacity)   
            if r(1,1)~=0
                for n = 1:height(r)
                    active_vector(r(n,1):c(n,1)-1,1) = zeros((c(n,1)-r(n,1)),1);
                end
            end
            %update the overall energy balance over a year
            Current_def = -1 * (cum_dummy_power(1) - cum_dummy_power(end)) + (cum_demand(1) - cum_demand(end));
        end
        S_curt = round(100 * sum(Curtailed .* fract_solar)/sum(power_orig.*fract_solar));
        W_curt = round(100 * sum(Curtailed .* (1-fract_solar))/sum(power_orig.*(1-fract_solar)));
        % End of sub-section to find the optimal amount of solar and wind
        % curtailment respectively at the optimal total curtailment
        
        k = 1;
        parfor k = 1:99 %cycle through each of 1-99% curtailment
            disp(k)
            buffCap = H2buffers(k,y) * S(k)/100; %set the buffer capacity based on the optimization result scale by the percent curtailment
            Curtailed = zeros(8760,1); % initialize array of curtailment for each hour over a year in GWh
            demand = single((S(k)/100) * 0.1 * (612.5/640) * ones(8760,1)); %set the hourly demand matrix in GW, scaled according to the percent curtailment
            power = single(panels(end,y) * panelPowers{y} * (1/1E9) + turbines(end,y) * turbinePowers{y} * (1/1E9)); %initialize the hourly power profile
            fract_solar = (panels(end,y) * panelPowers{y})./(panels(end,y)*panelPowers{y} + turbines(end,y)*turbinePowers{y})
            fract_solar(isnan(fract_solar)) = 0;

            power_orig = power; %set the original power profile to use later
            battery_use = battery_uses_year{k,1}; %extract the particular hourly battery use over a year for a given location and year
            battery_use = (battery_use + (27.5/640)*PReq*(1/1000)) * S(k)/100; %baseline electrical requirement added to battery usage, making is electricity extraction, but keep name as battery
            battery_use(battery_use < 0 ) = 0; %if battery use is negative, set it to zero.
            power = power - battery_use; %power supply after removing the electricity extraction
            elect_curtail = power - H2Sizes(k,y)*S(k)/100; %calculate the amount of electricity curtailed due to it being above the rated capacity of the electrolyser
            elect_curtail(elect_curtail < 0) = 0; %remove entries less than zeros (keep only those where electricity must be curtailed)
            Curtailed = Curtailed + elect_curtail; %add the curtailed electricity to the curtailment array
            power(power > (H2Sizes(k,y)*S(k)/100)) = H2Sizes(k,y)*S(k)/100; %adjust the power supply according to the maximum power set by the electrolyser capacity
            cum_energy = [0;cumsum(power)]; % create a vector of cumulative energy supply for each hour over a year
            cum_demand = [0;cumsum(demand)]; % create a vector of the cumulative energy demand (as hydrogen) for each hour over a year
            Current_def = -1 * (cum_energy(1) - cum_energy(end)) + (cum_demand(1) - cum_demand(end)); % calculate the initial overall energy balance (should be positive)
            
            active_vector = ones(8760,1); %create a vector which describes which hours during the year are "active", meaning more energy during the hour can still be curtailed
            dummy_power = power; %initialize dummy array of power supply
            cum_dummy_power = [0;cumsum(dummy_power)]; % calculate cumulative dummy power

            while Current_def > 0.001 % continue to curtail energy until the overall energy balance over a year is close to zero
                Avail_power = dummy_power .* active_vector; % calculate an array of the available power
                cum_Avail_power = [0;cumsum(Avail_power)]; % array of cumulative available power
                % create a matrix of the maxaimum percent curtailment which
                % can occur between the starting hour (row entry) and the
                % ending hour (column entry) before the energy deficit over
                % that period is equal to the hydrogen storage capacity
                max_curtail = ((-1 * triu(cum_dummy_power - cum_dummy_power') + triu(cum_demand - cum_demand')) + buffCap*ones(8761,8761))./(-1 * triu(cum_Avail_power - cum_Avail_power'));
                
                max_curtail(max_curtail<=0) = NaN; % remove entries less than zero
                curt = min(min(max_curtail)); % find the overall maximum amount of curtailment before the hydrogen storage capacity is reached
                % if the maximum curtailment is less than the amount of
                % curtailment which will cause the energy balance over the
                % year to close, record the starting and ending indexes.
                % Else, set the percent curtailment the amount which will
                % close the energy balance. 
                if curt < (Current_def/(-1*cum_Avail_power(1) + cum_Avail_power(end)))
                    [r,c] = find(max_curtail == curt);
                else
                    curt = Current_def/(-1+cum_Avail_power(1) + cum_Avail_power(end));
                    r = [0];
                    c = [0];
                end
                % add the hourly energy curtailed in this iteration to the
                % array of total curtailed energy
                Curtailed = Curtailed + Avail_power * curt;
                dummy_power = dummy_power - (Avail_power * curt); %update the dummy power array after removing curtailed energy
                cum_dummy_power = [0;cumsum(dummy_power)]; % update cumulative power array
                % if the percent curtailment is not set by closing the
                % yearly energy balance, then remove the active indexes
                % which reached their maximum level of curtailment (the
                % deficit over that period is equal to the hydrogen storage
                % capacity)
                if r(1,1)~=0
                    for n = 1:height(r)
                        active_vector(r(n,1):c(n,1)-1,1) = zeros((c(n,1)-r(n,1)),1);
                    end
                end
                %update the overall energy balance over a year
                Current_def = -1 * (cum_dummy_power(1) - cum_dummy_power(end)) + (cum_demand(1) - cum_demand(end));
            end
            % save the percent of energy used in each month for the current
            % iteration of total percent curtailment
            S_jan_used(k+1,1) = sum((power_orig(1:744) - Curtailed(1:744)).*fract_solar(1:744))/sum(power_orig(1:744).*fract_solar(1:744));
            S_feb_used(k+1,1) = sum((power_orig(745:1416) - Curtailed(745:1416)).*fract_solar(745:1416))/sum(power_orig(745:1416).*fract_solar(745:1416));
            S_mar_used(k+1,1) = sum((power_orig(1417:2160) - Curtailed(1417:2160)).*fract_solar(1417:2160))/sum(power_orig(1417:2160).*fract_solar(1417:2160));
            S_apr_used(k+1,1) = sum((power_orig(2161:2880) - Curtailed(2161:2880)).*fract_solar(2161:2880))/sum(power_orig(2161:2880).*fract_solar(2161:2880));
            S_may_used(k+1,1) = sum((power_orig(2881:3624) - Curtailed(2881:3624)).*fract_solar(2881:3624))/sum(power_orig(2881:3624).*fract_solar(2881:3624));
            S_jun_used(k+1,1) = sum((power_orig(3625:4344) - Curtailed(3625:4344)).*fract_solar(3625:4344))/sum(power_orig(3625:4344).*fract_solar(3625:4344));
            S_jul_used(k+1,1) = sum((power_orig(4345:5088) - Curtailed(4345:5088)).*fract_solar(4345:5088))/sum(power_orig(4345:5088).*fract_solar(4345:5088));
            S_aug_used(k+1,1) = sum((power_orig(5089:5832) - Curtailed(5089:5832)).*fract_solar(5089:5832))/sum(power_orig(5089:5832).*fract_solar(5089:5832));
            S_sep_used(k+1,1) = sum((power_orig(5833:6552) - Curtailed(5833:6552)).*fract_solar(5833:6552))/sum(power_orig(5833:6552).*fract_solar(5833:6552));
            S_oct_used(k+1,1) = sum((power_orig(6553:7296) - Curtailed(6553:7296)).*fract_solar(6553:7296))/sum(power_orig(6553:7296).*fract_solar(6553:7296));
            S_nov_used(k+1,1) = sum((power_orig(7297:8016) - Curtailed(7297:8016)).*fract_solar(7297:8016))/sum(power_orig(7297:8016).*fract_solar(7297:8016));
            S_dec_used(k+1,1) = sum((power_orig(8017:8760) - Curtailed(8017:8760)).*fract_solar(8017:8760))/sum(power_orig(8017:8760).*fract_solar(8017:8760));
            
            W_jan_used(k+1,1) = sum((power_orig(1:744) - Curtailed(1:744)).*(1-fract_solar(1:744)))/sum(power_orig(1:744).*(1-fract_solar(1:744)));
            W_feb_used(k+1,1) = sum((power_orig(745:1416) - Curtailed(745:1416)).*(1-fract_solar(745:1416)))/sum(power_orig(745:1416).*(1-fract_solar(745:1416)));
            W_mar_used(k+1,1) = sum((power_orig(1417:2160) - Curtailed(1417:2160)).*(1-fract_solar(1417:2160)))/sum(power_orig(1417:2160).*(1-fract_solar(1417:2160)));
            W_apr_used(k+1,1) = sum((power_orig(2161:2880) - Curtailed(2161:2880)).*(1-fract_solar(2161:2880)))/sum(power_orig(2161:2880).*(1-fract_solar(2161:2880)));
            W_may_used(k+1,1) = sum((power_orig(2881:3624) - Curtailed(2881:3624)).*(1-fract_solar(2881:3624)))/sum(power_orig(2881:3624).*(1-fract_solar(2881:3624)));
            W_jun_used(k+1,1) = sum((power_orig(3625:4344) - Curtailed(3625:4344)).*(1-fract_solar(3625:4344)))/sum(power_orig(3625:4344).*(1-fract_solar(3625:4344)));
            W_jul_used(k+1,1) = sum((power_orig(4345:5088) - Curtailed(4345:5088)).*(1-fract_solar(4345:5088)))/sum(power_orig(4345:5088).*(1-fract_solar(4345:5088)));
            W_aug_used(k+1,1) = sum((power_orig(5089:5832) - Curtailed(5089:5832)).*(1-fract_solar(5089:5832)))/sum(power_orig(5089:5832).*(1-fract_solar(5089:5832)));
            W_sep_used(k+1,1) = sum((power_orig(5833:6552) - Curtailed(5833:6552)).*(1-fract_solar(5833:6552)))/sum(power_orig(5833:6552).*(1-fract_solar(5833:6552)));
            W_oct_used(k+1,1) = sum((power_orig(6553:7296) - Curtailed(6553:7296)).*(1-fract_solar(6553:7296)))/sum(power_orig(6553:7296).*(1-fract_solar(6553:7296)));
            W_nov_used(k+1,1) = sum((power_orig(7297:8016) - Curtailed(7297:8016)).*(1-fract_solar(7297:8016)))/sum(power_orig(7297:8016).*(1-fract_solar(7297:8016)));
            W_dec_used(k+1,1) = sum((power_orig(8017:8760) - Curtailed(8017:8760)).*(1-fract_solar(8017:8760)))/sum(power_orig(8017:8760).*(1-fract_solar(8017:8760)));
        end
        % create matrix of percent of energy used each month for each of
        % the total percent curtailments
        S_month_used = [S_jan_used,S_feb_used,S_mar_used,S_apr_used,S_may_used,S_jun_used,S_jul_used,S_aug_used,S_sep_used,S_oct_used,S_nov_used,S_dec_used];
        W_month_used = [W_jan_used,W_feb_used,W_mar_used,W_apr_used,W_may_used,W_jun_used,W_jul_used,W_aug_used,W_sep_used,W_oct_used,W_nov_used,W_dec_used];
        %difference between the curtailment in the month at a given curtailment level and the previous level
        % add the final difference in energy used for the case case of 0% curtailment
        S_usedDiff = [S_month_used(2:end,:) - S_month_used(1:end-1,:)]; %difference between the curtailment in the month at a given curtailment level and the previous level
        S_usedDiff = [S_usedDiff ; 1-sum(S_usedDiff)];
        W_usedDiff = [W_month_used(2:end,:) - W_month_used(1:end-1,:)]; %difference between the curtailment in the month at a given curtailment level and the previous level
        W_usedDiff = [W_usedDiff ; 1-sum(W_usedDiff)];
        % remove entries that at NaN due to there being no solar or wind
        % installed
        W_usedDiff(isnan(W_usedDiff)) = 0;
        S_usedDiff(isnan(S_usedDiff)) = 0;

        costStep = CP_costs.top_yearlyLCOU{i}(:,y); %cost of utilizing energy at each percent curtailment level
        LAEC = CP_costs.top_yearlyLAEC{i}(:,y); % total energy cost (except the energy supply cost) at each percent curtailment level
        S_E_cost = CP_costs.top_solaronly_Ecost(i,y); % solar energy supply cost in $/MWh
        W_E_cost = CP_costs.top_windonly_Ecost(i,y); % wind energy supply cost in $/MWh
        % total value of energy in the month is the sum of the difference in curtailment multiplied by the value of the energy wasted
        monthlyValueTotal = sum(S_usedDiff .* (market_value - costStep)) + sum(W_usedDiff .* (market_value - costStep)); %value of energy in the month is the sum of the difference in curtailment multiplied by the value of the energy wasted
        % value of the used energy in the monthly is value of the energy
        % used multiplied by the change in the amount of energy for each
        % percent curtailment, only up to the optimial percent curtailment,
        % normalized by the total energy used in the month
        S_monthlyValueUsed = sum(S_usedDiff(1:(100-opt_curtailment(y)),:) .* (market_value - costStep(1:(100-opt_curtailment(y)),:))) ./ sum(S_usedDiff(1:(100-opt_curtailment(y)),:));
        % similarly for the value of curtailed energy, only analysing the
        % percent curtailment between the optimal energy utilization and
        % 100% energy utilization
        S_monthlyValueWasted = sum(S_usedDiff((101-opt_curtailment(y)):end,:) .* (market_value - costStep((101-opt_curtailment(y)):end,:))) ./ sum(S_usedDiff((101-opt_curtailment(y)):end,:));
        % total cost of the energy used, including the energy supply cost
        % and the cost of utilizing energy
        S_monthlyCostUsed = sum(S_usedDiff(1:(100-opt_curtailment(y)),:) .* ((S_E_cost*(100/(100-S_curt))) + LAEC(1:(100-opt_curtailment(y)),:))) ./ sum(S_usedDiff(1:(100-opt_curtailment(y)),:));
        % total fraction of energy used in each month
        S_monthlyXUsed = sum(S_usedDiff(1:(100-opt_curtailment(y)),:));
        % repeat for the wind energy 
        W_monthlyValueUsed = sum(W_usedDiff(1:(100-opt_curtailment(y)),:) .* (market_value - costStep(1:(100-opt_curtailment(y)),:))) ./ sum(W_usedDiff(1:(100-opt_curtailment(y)),:));
        W_monthlyValueWasted = sum(W_usedDiff((101-opt_curtailment(y)):end,:) .* (market_value - costStep((101-opt_curtailment(y)):end,:))) ./ sum(W_usedDiff((101-opt_curtailment(y)):end,:));
        W_monthlyCostUsed = sum(W_usedDiff(1:(100-opt_curtailment(y)),:) .* ((W_E_cost*(100/(100-W_curt))) + LAEC(1:(100-opt_curtailment(y)),:))) ./ sum(W_usedDiff(1:(100-opt_curtailment(y)),:));
        W_monthlyXUsed = sum(W_usedDiff(1:(100-opt_curtailment(y)),:));
        % remove the NaN entries
        S_monthlyValueUsed(isnan(S_monthlyValueUsed)) = 0;
        S_monthlyValueWasted(isnan(S_monthlyValueWasted)) = 0;
        S_monthlyCostUsed(isnan(S_monthlyCostUsed)) = 0;
        W_monthlyValueUsed(isnan(W_monthlyValueUsed)) = 0;
        W_monthlyValueWasted(isnan(W_monthlyValueWasted)) = 0;
        W_monthlyCostUsed(isnan(W_monthlyCostUsed)) = 0;
        %save the above metrics in the monthlyResults table
        monthlyResults{y,2:end} = monthlyValueTotal; %add results for that year to the record table
        monthlyResults{y+7,2:end} = Solar_E;
        monthlyResults{y+14,2:end} = S_monthlyValueUsed;
        monthlyResults{y+21,2:end} = S_monthlyValueWasted;
        monthlyResults{y+28,2:end} = S_monthlyCostUsed;
        monthlyResults{y+35,2:end} = S_monthlyXUsed;
        monthlyResults{y+42,2:end} = Wind_E;
        monthlyResults{y+49,2:end} = W_monthlyValueUsed;
        monthlyResults{y+56,2:end} = W_monthlyValueWasted;
        monthlyResults{y+63,2:end} = W_monthlyCostUsed;
        monthlyResults{y+70,2:end} = W_monthlyXUsed;
    end
    % find the average of the metrics over the 6 years, removing the zero
    % entries (i.e. only consider the year as counting toward the average
    % if solar or wind were installed for that year)
    for m = 2:13
        monthlyResults{7,m} = mean(nonzeros(monthlyResults{1:6,m}));
        monthlyResults{14,m} = mean(nonzeros(monthlyResults{8:13,m}));
        monthlyResults{21,m} = mean(nonzeros(monthlyResults{15:20,m}));
        monthlyResults{28,m} = mean(nonzeros(monthlyResults{22:27,m}));
        monthlyResults{35,m} = mean(nonzeros(monthlyResults{29:34,m}));
        monthlyResults{42,m} = mean(nonzeros(monthlyResults{36:41,m}));
        monthlyResults{49,m} = mean(nonzeros(monthlyResults{43:48,m}));
        monthlyResults{56,m} = mean(nonzeros(monthlyResults{50:55,m}));
        monthlyResults{63,m} = mean(nonzeros(monthlyResults{57:62,m}));
        monthlyResults{70,m} = mean(nonzeros(monthlyResults{64:69,m}));
        monthlyResults{77,m} = mean(nonzeros(monthlyResults{71:76,m}));
    end
    CP_costs.top_monthlyValues{i} = monthlyResults;
end
save('CP_costs','CP_costs')
clearvars