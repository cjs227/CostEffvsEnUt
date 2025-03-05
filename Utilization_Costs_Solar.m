% 1. Load the necessary data 
load('SP_base')
load('SP_nocurt')
load('SP_costs') 
% SP_costs structure already needs to have lists of the top10, middle10,
% and bot10 locations. These locations are defined as being in the lowest
% 10% of LCOA, middle 45-55% of LCOA and highest 10% of LCOA for the case
% of no curtailment

PReq=100;
panelCapital = SP_base.panelCapital; % $/kWp
electrolyserCapital = SP_base.electrolyserCapital; % $/GW fed in
ASUCapital = SP_base.ASUCapital; %$/GW fed in
H2compCapital = SP_base.H2compCapital; %$/GW fed in
N2compCapital = SP_base.N2compCapital; %$/GW fed in
HBCapital = SP_base.HBCapital; %dollars for 230 t/day plant
H2storageCapital = SP_base.H2storageCapital; % $/kg
N2storageCapital = SP_base.N2storageCapital; %$/kg
BattstorageCapital = SP_base.BattstorageCapital; %$/kWh
BattpowerCapital = SP_base.BattpowerCapital; %$/GW fed in
panel_OnM = SP_base.panel_OnM; % Dollar per MW per year
elect_OnM = SP_base.elect_OnM;
batt_OnM = SP_base.batt_OnM; %$/MWh
HB_OnM = SP_base.HB_OnM; %dollar per year
discountRate = SP_base.discountRate;
opYear = SP_base.opYear; % years
annual=(SP_base.discountRate/(1-(1+SP_base.discountRate)^(-1*SP_base.opYear)));

% 2a. Solving algorithm to optimize the process design when setting the
% panel size to between 1-99% curtailment. Analysing average locations. 
SP_costs.middle10_results = struct;
for i = 1:10
    loc = SP_costs.middle10_loc(i,1);
    disp(loc)
    ind = SP_base.location == loc;
    panelPowers = SP_base.panelPower(ind);
    years = SP_base.year(ind);
    panelskWp = SP_nocurt.panelkWp(ind);
    clear panel_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6
        panelkWp_nocurt = panelskWp(y);
        panelPower = panelPowers{y};
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
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000)
            else
                panelkWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1} = 0;
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
    SP_costs.middle10_results.panel{i} = panel_save;
    SP_costs.middle10_results.H2Size{i} = H2Size_save;
    SP_costs.middle10_results.N2Size{i} = N2Size_save;
    SP_costs.middle10_results.BattSize{i} = BattSize_save;
    SP_costs.middle10_results.H2buffer{i} = H2buffer_save;
    SP_costs.middle10_results.N2buffer{i} = N2buffer_save;
    SP_costs.middle10_results.batterykWh{i} = battery_save;
    SP_costs.middle10_results.LCOA{i} = LCOA_save;
    SP_costs.middle10_results.batteryUse{i} = batteryUse_save;
end
% 2b. Same as 2a for the top10 locations. 
SP_costs.top10_results = struct;
for i = 1:10
    loc = SP_costs.top10_loc(i,1);
    disp(loc)
    ind = SP_base.location == loc;
    panelPowers = SP_base.panelPower(ind);
    years = SP_base.year(ind);
    panelskWp = SP_nocurt.panelkWp(ind);
    clear panel_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6
        panelkWp_nocurt = panelskWp(y);
        panelPower = panelPowers{y};
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
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000)
            else
                panelkWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1} = 0;
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
    SP_costs.top10_results.panel{i} = panel_save;
    SP_costs.top10_results.H2Size{i} = H2Size_save;
    SP_costs.top10_results.N2Size{i} = N2Size_save;
    SP_costs.top10_results.BattSize{i} = BattSize_save;
    SP_costs.top10_results.H2buffer{i} = H2buffer_save;
    SP_costs.top10_results.N2buffer{i} = N2buffer_save;
    SP_costs.top10_results.batterykWh{i} = battery_save;
    SP_costs.top10_results.LCOA{i} = LCOA_save;
    SP_costs.top10_results.batteryUse{i} = batteryUse_save;
end

% 2c. Same as 2a and 2b for bottom 10 locations.
SP_costs.bot10_results = struct;
for i = 1:10
    loc = SP_costs.bot10_loc(i,1);
    disp(loc)
    ind = SP_base.location == loc;
    panelPowers = SP_base.panelPower(ind);
    years = SP_base.year(ind);
    panelskWp = SP_nocurt.panelkWp(ind);
    clear panel_save H2Size_save N2Size_save BattSize_save H2buffer_save N2buffer_save battery_save LCOA_save batteryUse_save
    for y = 1:6
        panelkWp_nocurt = panelskWp(y);
        panelPower = panelPowers{y};
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
                batteryUse{j,1}=sol.Elec_extract - (27.5/640)*PReq*(1/1000)
            else
                panelkWp(j,1)=0;
                H2SizeGW(j,1)=0;%GW
                N2SizeGW(j,1)=0; %GW
                BattSizeGW(j,1)=0; %GW
                H2bufferkg(j,1)=0; %kg
                N2bufferkg(j,1)=0; %kg
                batterykWh(j,1)=0; %kWh
                LCOA(j,1)=0;
                batteryUse{j,1} = 0;
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
    SP_costs.bot10_results.panel{i} = panel_save;
    SP_costs.bot10_results.H2Size{i} = H2Size_save;
    SP_costs.bot10_results.N2Size{i} = N2Size_save;
    SP_costs.bot10_results.BattSize{i} = BattSize_save;
    SP_costs.bot10_results.H2buffer{i} = H2buffer_save;
    SP_costs.bot10_results.N2buffer{i} = N2buffer_save;
    SP_costs.bot10_results.batterykWh{i} = battery_save;
    SP_costs.bot10_results.LCOA{i} = LCOA_save;
    SP_costs.bot10_results.batteryUse{i} = batteryUse_save;
end

% 3a. Calculate the cost of utilization metrics for top 10 locations using optimization
% results.
clearvars -except SP_costs panelCapital panel_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    panels = flip(SP_costs.top10_results.panel{i},1);
    H2Sizes = flip(SP_costs.top10_results.H2Size{i},1);
    N2Sizes = flip(SP_costs.top10_results.N2Size{i},1);
    BattSizes = flip(SP_costs.top10_results.BattSize{i},1);
    H2buffers = flip(SP_costs.top10_results.H2buffer{i},1);
    N2buffers = flip(SP_costs.top10_results.N2buffer{i},1);
    batteries = flip(SP_costs.top10_results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel size for full energy
        %utilization
        E_cost = SP_costs.top10_results.panel{i}(1,y) * (annual * panelCapital + (1/1000)*panel_OnM) /(876000);
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
        %save the energy cost for the location and year in the SP_value
        %structure
        SP_costs.top10_Ecost(i,y) = E_cost;
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
    %function of energy utilization. add to the SP_costs structure
    costs = table(mean(delta_panelCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["panel_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    SP_costs.top10_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    SP_costs.top10_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    SP_costs.top10_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end

% 3b. Calculate the costs of utilization metrics for the bottom 10 locations using
% optimization results. 
clearvars -except SP_costs panelCapital panel_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    panels = flip(SP_costs.bot10_results.panel{i},1);
    H2Sizes = flip(SP_costs.bot10_results.H2Size{i},1);
    N2Sizes = flip(SP_costs.bot10_results.N2Size{i},1);
    BattSizes = flip(SP_costs.bot10_results.BattSize{i},1);
    H2buffers = flip(SP_costs.bot10_results.H2buffer{i},1);
    N2buffers = flip(SP_costs.bot10_results.N2buffer{i},1);
    batteries = flip(SP_costs.bot10_results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel size for full energy
        %utilization
        E_cost = SP_costs.bot10_results.panel{i}(1,y) * (annual * panelCapital + (1/1000)*panel_OnM) /(876000);
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
        %save the energy cost for the location and year in the SP_value
        %structure
        SP_costs.bot10_Ecost(i,y) = E_cost;
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
    %function of energy utilization. add to the SP_costs structure
    costs = table(mean(delta_panelCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["panel_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    SP_costs.bot10_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    SP_costs.bot10_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    SP_costs.bot10_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end

% 3c. Calculate the utilization cost metrics for the average 10 locations using
% optimization results. 
clearvars -except SP_costs panelCapital panel_OnM electrolyserCapital H2compCapital elect_OnM BattpowerCapital H2storageCapital BattstorageCapital
S = linspace(1,100)';
for i = 1:10 %cycle through each of the 10 locations in the category
    %clear the variables from the previous iteration
    clear delta_panelCost delta_H2SizeCost delta_BattSizeCost delta_H2bufferCost delta_BatteryCost
    % load the results of the optimization for each percent curtailment,
    % flip so that it is in terms of increasing energy utilization
    panels = flip(SP_costs.middle10_results.panel{i},1);
    H2Sizes = flip(SP_costs.middle10_results.H2Size{i},1);
    N2Sizes = flip(SP_costs.middle10_results.N2Size{i},1);
    BattSizes = flip(SP_costs.middle10_results.BattSize{i},1);
    H2buffers = flip(SP_costs.middle10_results.H2buffer{i},1);
    N2buffers = flip(SP_costs.middle10_results.N2buffer{i},1);
    batteries = flip(SP_costs.middle10_results.batterykWh{i},1);
    for y = 1:6 %cycle through each of the years
        %get the energy cost based on the panel size for full energy
        %utilization
        E_cost = SP_costs.middle10_results.panel{i}(1,y) * (annual * panelCapital + (1/1000)*panel_OnM) /(876000);
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
        %save the energy cost for the location and year in the SP_value
        %structure
        SP_costs.middle10_Ecost(i,y) = E_cost;
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
    %function of energy utilization. add to the SP_costs structure
    costs = table(mean(delta_panelCost,2), ...
        mean(delta_H2SizeCost-delta_H2SizeCost_SS,2), ...
        mean(delta_BattSizeCost,2), ...
        mean(delta_H2bufferCost,2), ...
        mean(delta_BatteryCost,2), ...
        mean(delta_H2SizeCost_SS,2), ...
        mean(delta_ASUSizeCost_SS,2), ...
        mean(delta_HBSizeCost_SS,2), ...
        'VariableNames',["panel_cost","H2Size cost extra","BattSize cost","H2buffer cost","Battery cost","H2Size SS","ASU SS","HB SS"]);
    SP_costs.middle10_costs{i} = costs;
    % find the total cost of utilizing energy as a function of percent
    % utilization
    SP_costs.middle10_yearlyLCOU{i} = delta_H2SizeCost ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost ...
        + delta_ASUSizeCost_SS ...
        + delta_HBSizeCost_SS;
    % find the total energy cost as a function of percent energy
    % utilization, not including the cost of energy supply. 
    SP_costs.middle10_yearlyLAEC{i} = delta_H2SizeCost ...
        - delta_H2SizeCost_SS ...
        + delta_BattSizeCost ...
        + delta_H2bufferCost ...
        + delta_BatteryCost;
end
save('SP_costs','SP_costs')
clearvars

% 4. Calculate the monthly costs and values using the previous optimization
% results and utilization costs
load('SP_base')
load('SP_nocurt')
load('SP_curt')
load('SP_costs')
S = linspace(1,100)'; %list is a scaling array for percent curtailment. 
% It is used to resize all of the process units for each level of curtailment 
% because the optimization assumes constant total production rather than 
% constant solar panel size
market_value = 1000 * (17/1000000)*(3600000/640); %$ per MWh. given a price of green ammonia as $1000/tonne
PReq = 100; %nomial size 100 MW
for i = 5:5 % configured to only analyse one location (index 5), could be changed to cycle through several locations
    %load the process design data from the optimization of 0-99%
    %curtailment
    panels = flip(SP_costs.top10_results.panel{i},1);
    H2Sizes = flip(SP_costs.top10_results.H2Size{i},1); % GW
    BattSizes = flip(SP_costs.top10_results.BattSize{i},1); % GW
    H2buffers = flip(SP_costs.top10_results.H2buffer{i},1) * (0.2045) * 0.000278; % GWh
    batteries = flip(SP_costs.top10_results.batterykWh{i},1) * (1/1E6); % GWh
    battery_uses = (SP_costs.top10_results.batteryUse{i}); % GW
    % extract the panel power and the optimal amount of curtailment at the
    % location and year
    loc = SP_costs.top10_loc(i,1);
    disp(loc)
    ind = SP_base.location == loc;
    panelPowers = SP_base.panelPower(ind);
    opt_curtailment = round(SP_curt.curtailment(ind)*100);
    
    % initialize a table of monthly results for each of the six years
    month = ["Value_Type","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
    varTypes = ["string","double","double","double","double","double","double","double","double","double","double","double","double"];
    monthlyResults = table('Size', [35 numel(month)], 'VariableTypes',varTypes,'VariableNames', month);
    monthlyResults.Value_Type = [
        "Year 1 Total";"Year 2 Total";"Year 3 Total";"Year 4 Total";"Year 5 Total";"Year 6 Total";"Average Total"; ...
        "Year 1 Used";"Year 2 Used";"Year 3 Used";"Year 4 Used";"Year 5 Used";"Year 6 Used";"Average Used";...
        "Year 1 Wasted";"Year 2 Wasted";"Year 3 Wasted";"Year 4 Wasted";"Year 5 Wasted";"Year 6 Wasted";"Average Wasted";
        "Year 1 LAEC";"Year 2 LAEC";"Year 3 LAEC";"Year 4 LAEC";"Year 5 LAEC";"Year 6 LAEC";"Average LAEC"; ...
        "Year 1 Xused";"Year 2 Xused";"Year 3 Xused";"Year 4 Xused";"Year 5 Xused";"Year 6 Xused";"Average Xused"];

    for y = 1:6 %cycle through each of the years
        battery_uses_year = flip(battery_uses{1,y},1); %extract and flip the optimization result for the battery usage over a year (includes all 6 years at a location)
        %initialize arrays for the usage of energy in each monthly for each of the optimizations between 0-99% curtailment 
        jan_used = zeros(100,1);
        feb_used = zeros(100,1);
        mar_used = zeros(100,1);
        apr_used = zeros(100,1);
        may_used = zeros(100,1);
        jun_used = zeros(100,1);
        jul_used = zeros(100,1);
        aug_used = zeros(100,1);
        sep_used = zeros(100,1);
        oct_used = zeros(100,1);
        nov_used = zeros(100,1);
        dec_used = zeros(100,1);
        for k = 1:99 %cycle through each of 1-99% curtailment
            disp(k)
            buffCap = H2buffers(k,y) * S(k)/100; %set the buffer capacity based on the optimization result scale by the percent curtailment
            Curtailed = zeros(8760,1); % initialize array of curtailment for each hour over a year in GWh
            demand = single((S(k)/100) * 0.1 * (612.5/640) * ones(8760,1)); %set the hourly demand matrix in GW, scaled according to the percent curtailment
            power = single(panels(end,y) * panelPowers{y} * (1/1E9)); %initialize the hourly power profile provided by solar panels, according to the panel size for no curtailment 
            power_orig = power; %set the original power profile to use later
            battery_use = battery_uses_year{k,1}; %extract the particular hourly battery use over a year for a given location and year
            battery_use = (battery_use + (27.5/640)*PReq*(1/1000)) * S(k)/100; %baseline electrical requirement added to battery usage, making it electricity extraction, but keep name as battery usage
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
            jan_used(k+1,1) = sum(power_orig(1:744) - Curtailed(1:744))/sum(power_orig(1:744));
            feb_used(k+1,1) = sum(power_orig(745:1416) - Curtailed(745:1416))/sum(power_orig(745:1416));
            mar_used(k+1,1) = sum(power_orig(1417:2160) - Curtailed(1417:2160))/sum(power_orig(1417:2160));
            apr_used(k+1,1) = sum(power_orig(2161:2880) - Curtailed(2161:2880))/sum(power_orig(2161:2880));
            may_used(k+1,1) = sum(power_orig(2881:3624) - Curtailed(2881:3624))/sum(power_orig(2881:3624));
            jun_used(k+1,1) = sum(power_orig(3625:4344) - Curtailed(3625:4344))/sum(power_orig(3625:4344));
            jul_used(k+1,1) = sum(power_orig(4345:5088) - Curtailed(4345:5088))/sum(power_orig(4345:5088));
            aug_used(k+1,1) = sum(power_orig(5089:5832) - Curtailed(5089:5832))/sum(power_orig(5089:5832));
            sep_used(k+1,1) = sum(power_orig(5833:6552) - Curtailed(5833:6552))/sum(power_orig(5833:6552));
            oct_used(k+1,1) = sum(power_orig(6553:7296) - Curtailed(6553:7296))/sum(power_orig(6553:7296));
            nov_used(k+1,1) = sum(power_orig(7297:8016) - Curtailed(7297:8016))/sum(power_orig(7297:8016));
            dec_used(k+1,1) = sum(power_orig(8017:8760) - Curtailed(8017:8760))/sum(power_orig(8017:8760));
        end
        % create matrix of percent of energy used each month for each of
        % the total percent curtailments
        month_used = [jan_used,feb_used,mar_used,apr_used,may_used,jun_used,jul_used,aug_used,sep_used,oct_used,nov_used,dec_used];
        usedDiff = [month_used(2:end,:) - month_used(1:end-1,:)]; %difference between the curtailment in the month at a given curtailment level and the previous level
        usedDiff = [usedDiff ; 1-sum(usedDiff)]; % add the final difference in energy used for the case case of 0% curtailment
        
        costStep = SP_costs.top10_yearlyLCOU{i}(:,y); %cost of utilizing energy at each percent curtailment level
        LAEC = SP_costs.top10_yearlyLAEC{i}(:,y); % total energy cost (except the energy supply cost) at each percent curtailment level
        E_cost = SP_costs.top10_Ecost(i,y); % energy supply cost in $/MWh
        
        % total value of energy in the month is the sum of the difference in curtailment multiplied by the value of the energy wasted
        monthlyValueTotal = sum(usedDiff .* (market_value - costStep));
        % value of the used energy in the monthly is value of the energy
        % used multiplied by the change in the amount of energy for each
        % percent curtailment, only up to the optimial percent curtailment,
        % normalized by the total energy used in the month
        monthlyValueUsed = sum(usedDiff(1:(100-opt_curtailment(y)),:) .* (market_value - costStep(1:(100-opt_curtailment(y)),:))) ./ sum(usedDiff(1:(100-opt_curtailment(y)),:));
        % similarly for the value of curtailed energy, only analysing the
        % percent curtailment between the optimal energy utilization and
        % 100% energy utilization
        monthlyValueWasted = sum(usedDiff((101-opt_curtailment(y)):end,:) .* (market_value - costStep((101-opt_curtailment(y)):end,:))) ./ sum(usedDiff((101-opt_curtailment(y)):end,:));
        % total cost of the energy used, including the energy supply cost
        % and the cost of utilizing energy
        monthlyCostUsed = sum(usedDiff(1:(100-opt_curtailment(y)),:) .* ((E_cost*(100/(100-opt_curtailment(y)))) + LAEC(1:(100-opt_curtailment(y)),:))) ./ sum(usedDiff(1:(100-opt_curtailment(y)),:));
        % total fraction of energy used in each month
        monthlyXUsed = sum(usedDiff(1:(100-opt_curtailment(y)),:));
        %save the above metrics in the monthlyResults table
        monthlyResults{y,2:end} = monthlyValueTotal; 
        monthlyResults{y+7,2:end} = monthlyValueUsed;
        monthlyResults{y+14,2:end} = monthlyValueWasted;
        monthlyResults{y+21,2:end} = monthlyCostUsed;
        monthlyResults{y+28,2:end} = monthlyXUsed;
    end
    %calculate the mean of the 6 years of monthlyResults metrics
    monthlyResults{7,2:end} = mean(monthlyResults{1:6,2:end});
    monthlyResults{14,2:end} = mean(monthlyResults{8:13,2:end});
    monthlyResults{21,2:end} = mean(monthlyResults{15:20,2:end});
    monthlyResults{28,2:end} = mean(monthlyResults{22:27,2:end});
    monthlyResults{35,2:end} = mean(monthlyResults{29:34,2:end});
    SP_costs.top10_monthlyValues{i} = monthlyResults;
end
clearvars -except SP_costs SP_base SP_nocurt SP_curt
S = linspace(1,100)'; %list is a scaling array for percent curtailment. 
% It is used to resize all of the process units for each level of curtailment 
% because the optimization assumes constant total production rather than 
% constant solar panel size
market_value = 1000 * (17/1000000)*(3600000/640); %$ per MWh. given a price of green ammonia as $1000/tonne
PReq = 100; %nomial size 100 MW
for i = 9:9 % configured to only analyse one location (index 9), could be changed to cycle through several locations
    %load the process design data from the optimization of 0-99%
    %curtailment
    panels = flip(SP_costs.middle10_results.panel{i},1);
    H2Sizes = flip(SP_costs.middle10_results.H2Size{i},1); % GW
    BattSizes = flip(SP_costs.middle10_results.BattSize{i},1); % GW
    H2buffers = flip(SP_costs.middle10_results.H2buffer{i},1) * (0.2045) * 0.000278; % GWh
    batteries = flip(SP_costs.middle10_results.batterykWh{i},1) * (1/1E6); % GWh
    battery_uses = (SP_costs.middle10_results.batteryUse{i}); % GW
    % extract the panel power and the optimal amount of curtailment at the
    % location and year
    loc = SP_costs.middle10_loc(i,1);
    disp(loc)
    ind = SP_base.location == loc;
    panelPowers = SP_base.panelPower(ind);
    opt_curtailment = round(SP_curt.curtailment(ind)*100);
    
    % initialize a table of monthly results for each of the six years
    month = ["Value_Type","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
    varTypes = ["string","double","double","double","double","double","double","double","double","double","double","double","double"];
    monthlyResults = table('Size', [35 numel(month)], 'VariableTypes',varTypes,'VariableNames', month);
    monthlyResults.Value_Type = [
        "Year 1 Total";"Year 2 Total";"Year 3 Total";"Year 4 Total";"Year 5 Total";"Year 6 Total";"Average Total"; ...
        "Year 1 Used";"Year 2 Used";"Year 3 Used";"Year 4 Used";"Year 5 Used";"Year 6 Used";"Average Used";...
        "Year 1 Wasted";"Year 2 Wasted";"Year 3 Wasted";"Year 4 Wasted";"Year 5 Wasted";"Year 6 Wasted";"Average Wasted";
        "Year 1 LAEC";"Year 2 LAEC";"Year 3 LAEC";"Year 4 LAEC";"Year 5 LAEC";"Year 6 LAEC";"Average LAEC"; ...
        "Year 1 Xused";"Year 2 Xused";"Year 3 Xused";"Year 4 Xused";"Year 5 Xused";"Year 6 Xused";"Average Xused"];

    for y = 1:6 %cycle through each of the years
        battery_uses_year = flip(battery_uses{1,y},1); %extract and flip and optimization result for the battery usage over a year (includes all 6 years at a location)
        %initialize arrays for the usage of energy in each monthly for each of the optimizations between 0-99% curtailment 
        jan_used = zeros(100,1);
        feb_used = zeros(100,1);
        mar_used = zeros(100,1);
        apr_used = zeros(100,1);
        may_used = zeros(100,1);
        jun_used = zeros(100,1);
        jul_used = zeros(100,1);
        aug_used = zeros(100,1);
        sep_used = zeros(100,1);
        oct_used = zeros(100,1);
        nov_used = zeros(100,1);
        dec_used = zeros(100,1);
        for k = 1:99 %cycle through each of 1-99% curtailment
            disp(k)
            buffCap = H2buffers(k,y) * S(k)/100; %set the buffer capacity based on the optimization result scale by the percent curtailment
            Curtailed = zeros(8760,1); % initialize array of curtailment for each hour over a year in GWh
            demand = single((S(k)/100) * 0.1 * (612.5/640) * ones(8760,1)); %set the hourly demand matrix in GW, scaled according to the percent curtailment
            power = single(panels(end,y) * panelPowers{y} * (1/1E9)); %initialize the hourly power profile provided by solar panels, according to the panel size for no curtailment 
            power_orig = power; %set the original power profile to later use
            battery_use = battery_uses_year{k,1}; %extract the particular hourly battery use over a year for a given location and year
            battery_use = (battery_use + (27.5/640)*PReq*(1/1000)) * S(k)/100; %baseline electrical requirement added to battery usage, making is electricity extraction, but keep name as battery usage
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
            % save the percent of energy used in each month for the current
            % iteration of total percent curtailment
            jan_used(k+1,1) = sum(power_orig(1:744) - Curtailed(1:744))/sum(power_orig(1:744));
            feb_used(k+1,1) = sum(power_orig(745:1416) - Curtailed(745:1416))/sum(power_orig(745:1416));
            mar_used(k+1,1) = sum(power_orig(1417:2160) - Curtailed(1417:2160))/sum(power_orig(1417:2160));
            apr_used(k+1,1) = sum(power_orig(2161:2880) - Curtailed(2161:2880))/sum(power_orig(2161:2880));
            may_used(k+1,1) = sum(power_orig(2881:3624) - Curtailed(2881:3624))/sum(power_orig(2881:3624));
            jun_used(k+1,1) = sum(power_orig(3625:4344) - Curtailed(3625:4344))/sum(power_orig(3625:4344));
            jul_used(k+1,1) = sum(power_orig(4345:5088) - Curtailed(4345:5088))/sum(power_orig(4345:5088));
            aug_used(k+1,1) = sum(power_orig(5089:5832) - Curtailed(5089:5832))/sum(power_orig(5089:5832));
            sep_used(k+1,1) = sum(power_orig(5833:6552) - Curtailed(5833:6552))/sum(power_orig(5833:6552));
            oct_used(k+1,1) = sum(power_orig(6553:7296) - Curtailed(6553:7296))/sum(power_orig(6553:7296));
            nov_used(k+1,1) = sum(power_orig(7297:8016) - Curtailed(7297:8016))/sum(power_orig(7297:8016));
            dec_used(k+1,1) = sum(power_orig(8017:8760) - Curtailed(8017:8760))/sum(power_orig(8017:8760));
        end
        % create matrix of percent of energy used each month for each of
        % the total percent curtailments
        month_used = [jan_used,feb_used,mar_used,apr_used,may_used,jun_used,jul_used,aug_used,sep_used,oct_used,nov_used,dec_used];
        usedDiff = [month_used(2:end,:) - month_used(1:end-1,:)]; %difference between the curtailment in the month at a given curtailment level and the previous level
        usedDiff = [usedDiff ; 1-sum(usedDiff)]; % add the final difference in energy used for the case case of 0% curtailment
        
        costStep = SP_costs.middle10_yearlyLCOU{i}(:,y); %cost of utilizing energy at each percent curtailment level
        LAEC = SP_costs.middle10_yearlyLAEC{i}(:,y); % total energy cost (except the energy supply cost) at each percent curtailment level
        E_cost = SP_costs.middle10_Ecost(i,y); % energy supply cost in $/MWh
        
        % total value of energy in the month is the sum of the difference in curtailment multiplied by the value of the energy wasted
        monthlyValueTotal = sum(usedDiff .* (market_value - costStep));
        % value of the used energy in the monthly is value of the energy
        % used multiplied by the change in the amount of energy for each
        % percent curtailment, only up to the optimial percent curtailment,
        % normalized by the total energy used in the month
        monthlyValueUsed = sum(usedDiff(1:(100-opt_curtailment(y)),:) .* (market_value - costStep(1:(100-opt_curtailment(y)),:))) ./ sum(usedDiff(1:(100-opt_curtailment(y)),:));
        % similarly for the value of curtailed energy, only analysing the
        % percent curtailment between the optimal energy utilization and
        % 100% energy utilization
        monthlyValueWasted = sum(usedDiff((101-opt_curtailment(y)):end,:) .* (market_value - costStep((101-opt_curtailment(y)):end,:))) ./ sum(usedDiff((101-opt_curtailment(y)):end,:));
        % total cost of the energy used, including the energy supply cost
        % and the cost of utilizing energy
        monthlyCostUsed = sum(usedDiff(1:(100-opt_curtailment(y)),:) .* ((E_cost*(100/(100-opt_curtailment(y)))) + LAEC(1:(100-opt_curtailment(y)),:))) ./ sum(usedDiff(1:(100-opt_curtailment(y)),:));
        % total fraction of energy used in each month
        monthlyXUsed = sum(usedDiff(1:(100-opt_curtailment(y)),:));
        %save the above metrics in the monthlyResults table
        monthlyResults{y,2:end} = monthlyValueTotal; %add results for that year to the record table
        monthlyResults{y+7,2:end} = monthlyValueUsed;
        monthlyResults{y+14,2:end} = monthlyValueWasted;
        monthlyResults{y+21,2:end} = monthlyCostUsed;
        monthlyResults{y+28,2:end} = monthlyXUsed;
    end
    %calculate the mean of the 6 years of monthlyResults metrics
    monthlyResults{7,2:end} = mean(monthlyResults{1:6,2:end});
    monthlyResults{14,2:end} = mean(monthlyResults{8:13,2:end});
    monthlyResults{21,2:end} = mean(monthlyResults{15:20,2:end});
    monthlyResults{28,2:end} = mean(monthlyResults{22:27,2:end});
    monthlyResults{35,2:end} = mean(monthlyResults{29:34,2:end});
    SP_costs.middle10_monthlyValues{i} = monthlyResults;
end

save('SP_costs','SP_costs')
clearvars