% 1. Load necessary structures, initialize structure WP_curt, and extract cost parameters from the WP_base structure in variables. 
% (this is done so that the linear solver does not call WP_base each iteration)
clearvars
load('WP_base')
load('WP_nocurt')
WP_curt=struct; % initialize WP_curt structure
WP_curt.filename=WP_nocurt.filename; % transfer info from WP_nocurt
WP_curt.year=WP_nocurt.year;
WP_curt.location=WP_nocurt.location;
WP_curt.category=WP_nocurt.category;

PReq=100; % required size of the Hb plant in MW
electrolyserCapital = WP_base.electrolyserCapital; % turbine capital cost $/GW fed in
ASUCapital = WP_base.ASUCapital; % electrolyser capital cost $/GW fed in
H2compCapital = WP_base.H2compCapital; % H2 compressor cost in $/GW fed in
N2compCapital = WP_base.N2compCapital; % N2 compressor cost in $/GW fed in
HBCapital = WP_base.HBCapital; % HB capital cost in dollars for 230 t/day plant
H2storageCapital = WP_base.H2storageCapital; % H2 storage capital cost $/kg
N2storageCapital = WP_base.N2storageCapital; % N2 stoage capital cost in $/kg
BattstorageCapital = WP_base.BattstorageCapital; % battery energy capacity cost in $/kWh
BattpowerCapital = WP_base.BattpowerCapital; % batterey power capacity cost in $/GW fed in
elect_OnM = WP_base.elect_OnM; % electrolyser O&M in percent of total capital cost
HB_OnM = WP_base.HB_OnM; % fixed O&M for HB plant in percent of capital cost per year. The parameter is also used for O&M of compressors and ASU. 
discountRate = WP_base.discountRate; % discount rate of capital
opYear = WP_base.opYear; % operating years of the plant
Annualization=(WP_base.discountRate/(1-(1+WP_base.discountRate)^(-1*WP_base.opYear))); % capital annualization factor 

% create array of turbine capital cost (turbineCapital_store) and operating
% cost in percent of total cost (turbine_OnM_perc_store) according to the
% category of the wind power
for k=1:length(WP_curt.location)
    if WP_curt.category(k) == "onshore"
        turbineCapital_store(k) = WP_base.turbineCapital;
        turbine_OnM_perc_store(k) = WP_base.turbine_OnM_perc;
    elseif WP_curt.category(k) == "offshore" 
        turbineCapital_store(k) = WP_base.turbineCapital_offshore;
        turbine_OnM_perc_store(k) = WP_base.turbine_OnM_offshore_perc;
    elseif WP_curt.category(k) == "floating"
        turbineCapital_store(k) = WP_base.turbineCapital_floating;
        turbine_OnM_perc_store(k) = WP_base.turbine_OnM_floating_perc;
    end
end

% 2. Configure and solve linear problem to find the optimized process design with curtailment for each of the locations/years. The uses parrallelization to speed up execution.
% if needed, the parrallel loop can be exceuted in a normal for loop
% if the solver encounters errors, they are saved in the array "error_index". These iterations would need to be repeated. 
error_index=zeros(length(WP_curt.location),1); %initialize array which will store the indeces of any iterations that encounter an error
parfor i=1:length(WP_curt.location) % iterate through even entry in WP_curt using parrallel for loop
    turbineCapital=turbineCapital_store(i); % extract turbine capital cost from array
    turbine_OnM_perc=turbine_OnM_perc_store(i); % extract turbine O&M from array
    N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
    turbinePower=WP_base.turbinePower{i}; % extract the turbine power in W/kWp
    windProblem=optimproblem('ObjectiveSense','minimize'); % initialize optimization problem
    % Create variables
    E_size_wind=optimvar('E_size_wind',1,'LowerBound',0); % size of turbine relative to HB size (100 MW)
    E_supply=optimvar('E_supply',8760); % energy supply in GW
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
    Energy = E_supply == (PReq*E_size_wind*turbinePower)/1E6; % defining E supply variable in GW
    Energy_min = sum(E_supply) >= PReq*8760/1000; %minimum energy amount to close the energy balance
    powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be less than the energy supply
    
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
    Capex=Annualization*((E_size_wind*PReq*1000*turbineCapital) ...
        +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
        +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
        +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
        +(Batt_size*PReq*BattpowerCapital/1000) ...
        +HBCapital);
    % define operating cost
    Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
        +(turbine_OnM_perc*E_size_wind*PReq*1000*turbineCapital*Annualization) ...
        +HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
        +(HBCapital*HB_OnM);
    windProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
  
    % add constraints to the problem
    windProblem.Constraints.boundedH2 = boundedH2;
    windProblem.Constraints.boundedH2Initial = boundedH2Initial;
    windProblem.Constraints.boundedBatt = boundedBatt;
    windProblem.Constraints.boundedBattInitial = boundedBattInitial;
    windProblem.Constraints.Energy = Energy;
    windProblem.Constraints.powerBalance = powerBalance;
    windProblem.Constraints.Energy_min = Energy_min;
    windProblem.Constraints.H2Balance = H2Balance;
    windProblem.Constraints.BattBalance = BattBalance;
    windProblem.Constraints.H2finalBalance = H2finalBalance;
    windProblem.Constraints.BattfinalBalance = BattfinalBalance;
    windProblem.Constraints.Max_H2 = Max_H2;
    windProblem.Constraints.Max_Batt = Max_Batt;
    [sol,fval,exitflag,output] = solve(windProblem); % execute solver
    if exitflag==1 % if the solver is successful, add the optimized process component sizes to arrays
        turbinekWp(1,i)=sol.E_size_wind*PReq*1000;
        H2SizeGW(1,i)=sol.H2_size*PReq/1000;%GW
        N2SizeGW(1,i)=N2_size*PReq/1000; %GW
        BattSizeGW(1,i)=sol.Batt_size*PReq/1000; %GW
        H2bufferkg(1,i)=sol.H2_buffer*(1/0.2045); %kg
        N2bufferkg(1,i)=0; %kg
        batterykWh(1,i)=sol.Batt_buffer*277.78;
        LCOA(1,i)=fval;
    else % else, add 0 to process component size arrays and add index to error_index array
        turbinekWp(1,i)=0;
        H2SizeGW(1,i)=0;%GW
        N2SizeGW(1,i)=0; %GW
        BattSizeGW(1,i)=0; %GW
        H2bufferkg(1,i)=0; %kg
        N2bufferkg(1,i)=0; %kg
        batterykWh(1,i)=0; %kWh
        LCOA(1,i)=0;
        error_index(i,1)=i;
    end
end
% store the optimized process sizing to the WP_curt structure
WP_curt.turbinekWp=turbinekWp;
WP_curt.H2SizeGW=H2SizeGW;
WP_curt.N2SizeGW=N2SizeGW;
WP_curt.BattSizeGW=BattSizeGW;
WP_curt.H2bufferkg=H2bufferkg;
WP_curt.N2bufferkg=N2bufferkg;
WP_curt.batterykWh=batterykWh;
WP_curt.curtailment = (WP_curt.turbinekWp-WP_nocurt.turbinekWp)./WP_curt.turbinekWp; % find the fraction of curtailment based on the solar panel size in WP_nocurt
save('WP_curt','WP_curt') % save WP_curt structure

% 3. Find the cost of each of the process components with optimizal sizing and sort results according to location
% this code follows the same procedure as in Wind_Data_Processing.m
for k = 1:length(WP_curt.turbinekWp) % create arrays for turbine capital and O&M cost
    if WP_curt.category(k)=="onshore"
        turbineCost(k,1) = Annualization * WP_curt.turbinekWp(k) * WP_base.turbineCapital;
        turbineOnM(k,1) = turbineCost(k,1) * WP_base.turbine_OnM_perc;
    elseif WP_curt.category(k)=="offshore" 
        turbineCost(k,1) = Annualization * WP_curt.turbinekWp(k) * WP_base.turbineCapital_offshore;
        turbineOnM(k,1) =turbineCost(k,1) * WP_base.turbine_OnM_offshore_perc;
    elseif WP_curt.category(k)=="floating"
        turbineCost(k,1) = Annualization * WP_curt.turbinekWp(k) * WP_base.turbineCapital_floating;
        turbineOnM(k,1) = turbineCost(k,1) * WP_base.turbine_OnM_floating_perc;
    end
end
electrolyserCost = Annualization * WP_curt.H2SizeGW' * electrolyserCapital; % electrolyser capital
H2compCost = Annualization * WP_curt.H2SizeGW' * H2compCapital; % H2 compressor capital
ASUCost = Annualization * WP_curt.N2SizeGW' * ASUCapital; % PSA ASU capital 
N2compCost = Annualization * WP_curt.N2SizeGW' * N2compCapital; % N2 compressor capital
H2storeCost = Annualization *  WP_curt.H2bufferkg' * H2storageCapital; % H2 storage capital
N2storeCost = Annualization *  WP_curt.N2bufferkg' * N2storageCapital; % N2 storage capital
BattstoreCost = Annualization * (WP_curt.batterykWh' * BattstorageCapital + WP_curt.BattSizeGW' * BattpowerCapital); % battery capital
HBCost = ones(length(WP_curt.location),1) * Annualization * HBCapital; % HB capital
electOnM = elect_OnM * (WP_curt.H2SizeGW' * electrolyserCapital + WP_curt.H2SizeGW' * H2compCapital); % H2 production O&M
ASUOnM = HB_OnM * (WP_curt.N2SizeGW' * ASUCapital + WP_curt.N2SizeGW' * N2compCapital); % N2 production O&M
HBOnM = ones(length(WP_curt.location),1)* HB_OnM * HBCapital; % HB O&M
WP_curt.totalCapital = (turbineCost + electrolyserCost + H2compCost + ASUCost + N2compCost + N2storeCost + H2storeCost + BattstoreCost + HBCost); % total capital cost in Dollars
WP_curt.LCOA = (WP_curt.totalCapital + electOnM + turbineOnM + HBOnM+ASUOnM)/(230*365); % LCOA in Dollars/tonne
locs=WP_curt.location';
WP_curt.costBreakdown=table(locs,turbineCost,electrolyserCost,H2compCost,ASUCost,N2compCost,H2storeCost,N2storeCost,BattstoreCost,HBCost,turbineOnM,electOnM,ASUOnM,HBOnM); % cost components saved in table

% sort results of LCOA, fraction curtailment, total energy, and turbine size
% by location
LCOA_byLoc=[];
curt_byloc=[];
energy_byLoc=[];
turbinekWp_byLoc=[];
unique_loc=unique(WP_curt.location,'stable');
for i=1:length(unique_loc)
    for ii=1:length(WP_curt.location)
        if unique_loc(i)==WP_curt.location(ii)
            if WP_curt.year{ii}==2013
                LCOA_byLoc(i,1)=WP_curt.LCOA(ii);
                curt_byLoc(i,1)=WP_curt.curtailment(ii);
                energy_byLoc(i,1)=sum(WP_curt.turbinekWp(ii)*WP_base.turbinePower{ii}/1000);
                turbinekWp_byLoc(i,1)=WP_curt.turbinekWp(ii);
            elseif WP_curt.year{ii}==2014
                LCOA_byLoc(i,2)=WP_curt.LCOA(ii);
                curt_byLoc(i,2)=WP_curt.curtailment(ii);
                energy_byLoc(i,2)=sum(WP_curt.turbinekWp(ii)*WP_base.turbinePower{ii}/1000);
                turbinekWp_byLoc(i,2)=WP_curt.turbinekWp(ii);
            elseif WP_curt.year{ii}==2015
                LCOA_byLoc(i,3)=WP_curt.LCOA(ii);
                curt_byLoc(i,3)=WP_curt.curtailment(ii);
                energy_byLoc(i,3)=sum(WP_curt.turbinekWp(ii)*WP_base.turbinePower{ii}/1000);
                turbinekWp_byLoc(i,3)=WP_curt.turbinekWp(ii);
            elseif WP_curt.year{ii}==2016
                LCOA_byLoc(i,4)=WP_curt.LCOA(ii);
                curt_byLoc(i,4)=WP_curt.curtailment(ii);
                energy_byLoc(i,4)=sum(WP_curt.turbinekWp(ii)*WP_base.turbinePower{ii}/1000);
                turbinekWp_byLoc(i,4)=WP_curt.turbinekWp(ii);
            elseif WP_curt.year{ii}==2017
                LCOA_byLoc(i,5)=WP_curt.LCOA(ii);
                curt_byLoc(i,5)=WP_curt.curtailment(ii);
                energy_byLoc(i,5)=sum(WP_curt.turbinekWp(ii)*WP_base.turbinePower{ii}/1000);
                turbinekWp_byLoc(i,5)=WP_curt.turbinekWp(ii);
            elseif WP_curt.year{ii}==2018
                LCOA_byLoc(i,6)=WP_curt.LCOA(ii);
                curt_byLoc(i,6)=WP_curt.curtailment(ii);
                energy_byLoc(i,6)=sum(WP_curt.turbinekWp(ii)*WP_base.turbinePower{ii}/1000);
                turbinekWp_byLoc(i,6)=WP_curt.turbinekWp(ii);
            end
            unique_cat(i)=WP_curt.category(ii);
        end
    end
end
% create array that only include no floating turbines
k=1;
for i=1:length(unique_loc)
    if unique_cat(i) ~="floating"
        unique_loc_noFloat(k,:)=unique_loc(i);
        LCOA_byLoc_noFloat(k,:)=LCOA_byLoc(i,:);
        k=k+1;
    end
end
WP_curt.LCOA_byLoc=LCOA_byLoc;
WP_curt.unique_loc=unique_loc';
WP_curt.curt_byLoc=curt_byLoc';
WP_curt.unique_cat=unique_cat';
WP_curt.unique_loc_noFloat=unique_loc_noFloat;
WP_curt.LCOA_byLoc_noFloat=LCOA_byLoc_noFloat;
WP_curt.turbinekWp_byLoc=turbinekWp_byLoc;
WP_curt.energy_byLoc=energy_byLoc;
save('WP_curt','WP_curt');