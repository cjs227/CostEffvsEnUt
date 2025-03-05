% 1. Load necessary structures, initialize structure WP_ramp, and extract cost parameters from the WP_base structure in variables. 
% (this is done so that the linear solver does not call WP_base each iteration)
clearvars
load('WP_base')
load('WP_nocurt')
load('WP_curt')
WP_ramp=struct; % initialize the WP_ramp structure
WP_ramp.filename=WP_nocurt.filename; % transfer info from the WP_nocurt structure
WP_ramp.year=WP_nocurt.year;
WP_ramp.location=WP_nocurt.location;
WP_ramp.unique_loc=WP_nocurt.unique_loc;
WP_ramp.category=WP_nocurt.category;

PReq=100; % required size of HB plant in MW (in this case it is the equivalent average production with ramping)
electrolyserCapital = WP_base.electrolyserCapital; % electrolyser capital cost in $/GW fed in
ASUCapital = WP_base.ASUCapital; % PSA ASU capital cost in $/GW fed in
H2compCapital = WP_base.H2compCapital; % H2 compressor cost in $/GW fed in
N2compCapital = WP_base.N2compCapital; % N2 compressor cost in $/GW fed in
HBCapital = WP_base.HBCapital; % HB capital cost in dollars for 230 t/day plant
H2storageCapital = WP_base.H2storageCapital; % H2 storage capital cost $/kg
N2storageCapital = WP_base.N2storageCapital; % N2 storage capital cost $/kg
BattstorageCapital = WP_base.BattstorageCapital; % battery storage energy capacity cost $/kWh
BattpowerCapital = WP_base.BattpowerCapital; % battery storage power capacity cost $/GW fed in
elect_OnM = WP_base.elect_OnM; % electrolyser O&M in percent of total capital cost
HB_OnM = WP_base.HB_OnM; % fixed O&M for HB plant in percent of capital cost per year. The parameter is also used for O&M of compressors and ASU. 
discountRate = WP_base.discountRate; % discount rate of capital
opYear = WP_base.opYear; % operating years of the plant
annual=(WP_base.discountRate/(1-(1+WP_base.discountRate)^(-1*WP_base.opYear))); % annualization of capital cost
% create array of turbine capital cost (turbineCapital_store) and operating
% cost in percent of total cost (turbine_OnM_perc_store) according to the
% category of the wind power
for k=1:length(WP_ramp.location)
    if WP_ramp.category(k) == "onshore"
        turbineCapital_store(k) = WP_base.turbineCapital;
        turbine_OnM_perc_store(k) = WP_base.turbine_OnM_perc;
    elseif WP_ramp.category(k) == "offshore" 
        turbineCapital_store(k) = WP_base.turbineCapital_offshore;
        turbine_OnM_perc_store(k) = WP_base.turbine_OnM_offshore_perc;
    elseif WP_ramp.category(k) == "floating"
        turbineCapital_store(k) =  WP_base.turbineCapital_floating;
        turbine_OnM_perc_store(k) = WP_base.turbine_OnM_floating_perc;
    end
end

% 2. Configure and solve linear problem to find the optimized process design with curtailment for each of the locations/years, while allowing ramping of the HB process. This uses parrallelization to speed up execution.
% if needed, the parrallel loop can be exceuted as normal for loop 
% if the solver encounters errors, they are saved in the array "error_index". These iterations would need to be repeated. 
error_index=zeros(length(WP_ramp.location),1); % initialize array to stay index of iterations where solver returns error
minCap = 0.6; % minimum fraction capacity of the HB process (i.e. 60% of the design capacity)
rampRate = 0.2; %maximum change in fraction capacity per hour (i.e. the process can go up by 5% of design capacity per hour)

NH3Prodramp_wind = cell(1,length(WP_ramp.location)); % initialize cell array to store the arrays of HB production over time for each iteration
Extramp_wind = cell(1,length(WP_ramp.location)); % initialize cell array to store the arrays of energy extraction over time for each iteration
% initialize arrays to store the results of the solver
turbinekWp=zeros(1,length(WP_ramp.location));
H2SizeGW=zeros(1,length(WP_ramp.location));%GW
N2SizeGW=zeros(1,length(WP_ramp.location)); %GW
BattSizeGW=zeros(1,length(WP_ramp.location)); %GW
H2bufferkg=zeros(1,length(WP_ramp.location)); %kg
N2bufferkg=zeros(1,length(WP_ramp.location)); %kg
batterykWh=zeros(1,length(WP_ramp.location));
NH3Cap=zeros(1,length(WP_ramp.location)); % store the installed capacity of HB plant
LCOA=zeros(1,length(WP_ramp.location));

% Set an array of indices to covert the 4 hourly timescale of HB production
% levels to the hourly timescale of the power supply
h = 0;
for v = 1:2190
    for g = 1:4
        h = h+1;
        indices(h) = v;
    end
end
parfor i = 1:length(WP_ramp.location)
    turbinePower =WP_base.turbinePower{i}; %extract the turbine power in W/kWp
    N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
    % To speed the solver, upper limits for the turbine capacity, energy
    % extraction capacity, and energy storage capacity is set. 
    Up1 = 1.5 * WP_curt.turbinekWp(1,i)/(PReq*1000); % 1.5X the turbine capacity without ramping
    Up2 = 1.5 * WP_curt.H2SizeGW(1,i)*1000/(PReq); % 1.5X the energy extraction capacity without ramping
    Up3 = 1.5 * WP_curt.H2bufferkg(1,i)/((1/0.2045)); % the energy storagae capacity without ramping

    turbineCapital=turbineCapital_store(i); % extract turbine capital cost for current iteration from array
    turbine_OnM_perc=turbine_OnM_perc_store(i); % extract turbine operating cost for current iteration from array

    % The next five lines solve for the minimum wind turbine capacity (i.e. without curtailment) using
    % the same method as in "Wind_Data_Processing"
    turbinekWpGuess = 1e6; % kWp
    x0 = turbinekWpGuess;
    inventorySum = @(x) sum((turbinePower*x/1e9 - PReq/1000)*3600); % Genergy balance
    x = fzero(@(x) inventorySum (x), x0);
    windMin = x/(1000*PReq); % minimum turbine capacity realtive to PReq

    windProblem=optimproblem('ObjectiveSense','minimize'); % initialize optimization problem
    % Create variables
    NH3cap=optimvar('NH3cap',1,'LowerBound',1,'UpperBound',1/minCap); % size of HB process relative to 230 t/day plant. Maximum design capacity is 1/minCap because the total production is still 230 t/day overall
    E_size_wind=optimvar('E_size_wind',1,'LowerBound',windMin,'UpperBound',Up1); % size of turbines relative to average HB size (100 MW)
    E_supply=optimvar('E_supply',8760); % energy supply in GW
    H2_size=optimvar('H2_size',1,'LowerBound',1,'UpperBound',Up2); % size of H2 generation relative to the HB
    Batt_size = optimvar('Batt_size',1); %size of the battery power intake capacity relative to HB
    H2_buffer = optimvar('H2_buffer',1,'LowerBound',0,'UpperBound',Up3); % H2 buffer storage size in GJ
    Batt_buffer = optimvar('Batt_buffer',1,'LowerBound',0); % Battery buffer storage in GJ
    H2_initial = optimvar('H2_initial',1,'LowerBound',0); % initial level of the H2 buffer in GJ
    Batt_initial = optimvar('Batt_initial',1,'LowerBound',0); % initial level of the battery in GJ
    H2_store = optimvar('H2_store',8760,'LowerBound',0); % level of the H2 buffer at each hour over the year in GJ
    Batt_store = optimvar('Batt_store',8760,'LowerBound',0); % level of the battery at each hour over the year in GJ
    H2_extract = optimvar('H2_extract',8760,'LowerBound',0); % amount of energy fed to H2 production at each hour over the year in GW
    Elec_extract = optimvar('Elec_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)

    NH3prod=optimvar('NH3prod',2190); % production of HB plant over time in terms of MW of energy
    ramp=optimvar('ramp',2189); % ramp of HB plant between time points (one less point than total time points)

    % Create constraints
    Energy = E_supply == (PReq*E_size_wind*turbinePower)/1E6; % defining E supply variable in GW
    Energy_min = sum(E_supply) >= PReq*8760/1000; % minimum energy amount to close the energy balance
    powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be less than the energy supply
    Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
    Max_Batt = Elec_extract - (27.5/640)*NH3prod(indices(1:end))*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity. the battery power capacity could be slighly higher due to ramping HB down, but the HB load is approximated as constant for simplicity. 

    MaxCap = NH3prod <= NH3cap * PReq; % HB plant production must be less than installed capacity
    MinCap = NH3prod >= minCap*NH3cap * PReq; % HB plant production must be greater than minimum capacity
    ProdTotal = sum(NH3prod) * 4 == PReq * 8760; % total ammonia production must be equal to the equivalent in constant HB process (in terms of energy)
    
    H2Balance = optimconstr(8760); %initialize H2 balance of the year
    BattBalance = optimconstr(8760); % initialize battery balance over the year
    %H2 balance at all hours except the first
    H2Balance(2:end) = H2_store(2:end) == H2_store(1:end-1) +...
        (H2_extract(2:end) - NH3prod(indices(2:end))*(612.5/640000))*3600; % in GJ
    %H2 balance at the first hour
    H2Balance(1) = H2_store(1) == H2_initial +...
        (H2_extract(1) - NH3prod(indices(1))*(612.5/640000))*3600; 
    H2finalBalance = H2_store(end,:) == H2_initial; % final energy inventory must equal the initial inventory
    % battery balance at all hours except the first
    BattBalance(2:end) = Batt_store(2:end) == Batt_store(1:end-1) +...
        (Elec_extract(2:end) - NH3prod(indices(2:end))*(27.5/640000))*3600; % in GJ
    %battery balance at the first hour
    BattBalance(1) = Batt_store(1) == Batt_initial +...
        (Elec_extract(1)  - NH3prod(indices(1))*(27.5/640000))*3600; 
    BattfinalBalance = Batt_store(end,:) == Batt_initial; % final energy inventory must equal the initial inventory

    boundedH2 = H2_store <= H2_buffer; % H2 buffer level must be less than the H2 buffer size
    boundedH2Initial = H2_initial <= H2_buffer;
    boundedBatt = Batt_store <= Batt_buffer; % battery buffer level must be less than the battery size
    boundedBattInitial = Batt_initial <= Batt_buffer;
    
    NH3change = NH3prod(2:end) == NH3prod(1:end-1) + ramp(1:end); % HB process production at every time point is the previous time point plus the "ramp" of the process
    RampLow = ramp >= -1*rampRate*NH3cap*PReq; % HB process can only ramp down at up to the maximum ramp rate relative to design capacity
    RampHigh = ramp <= rampRate*NH3cap*PReq; % HB process can only ramp up at up to the maximum ramp rate relative to design capacity

    % define annualized capital cost
    Capex=annual*((E_size_wind*PReq*1000*turbineCapital) ...
            +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
            +(N2_size*PReq*NH3cap*(ASUCapital+N2compCapital)/1000) ...
            +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
            +(Batt_size*PReq*BattpowerCapital/1000) ...
            +HBCapital * NH3cap);
    % define operating cost
    Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
            +(turbine_OnM_perc*E_size_wind*PReq*1000*turbineCapital*annual) ...
            +(HB_OnM*(N2_size*PReq*NH3cap*(ASUCapital+N2compCapital)/1000)) ...
            +(NH3cap * HBCapital*HB_OnM);
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
    windProblem.Constraints.NH3Change = NH3change;
    windProblem.Constraints.RampLow = RampLow;
    windProblem.Constraints.RampHigh = RampHigh;
    windProblem.Constraints.MaxCap = MaxCap;
    windProblem.Constraints.MinCap = MinCap;
    windProblem.Constraints.ProdTotal = ProdTotal;
    [sol,fval,exitflag,output] = solve(windProblem); % execute solver
    if exitflag==1 % if the solver is successful, add the optimized process component sizes to arrays
        turbinekWp(1,i)=sol.E_size_wind*PReq*1000;
        H2SizeGW(1,i)=sol.H2_size*PReq/1000;%GW
        N2SizeGW(1,i)=sol.NH3cap*N2_size*PReq/1000; %GW
        BattSizeGW(1,i)=sol.Batt_size*PReq/1000; %GW
        H2bufferkg(1,i)=sol.H2_buffer*(1/0.2045); %kg
        N2bufferkg(1,i)=0; %kg
        batterykWh(1,i)=sol.Batt_buffer*277.78;
        NH3Cap(1,i)=PReq*sol.NH3cap; %installed capacity of HB plant
        LCOA(1,i)=fval;
        Extramp_wind{1,i} = sol.H2_extract + sol.Elec_extract; % store the energy extraction over time in cell array
        NH3Prodramp_wind{1,i} = sol.NH3prod; % store the HB production over time in cell array
    else % else, add 0 to process component size arrays and add index to error_index array
        turbinekWp(1,i)=0;
        H2SizeGW(1,i)=0;%GW
        N2SizeGW(1,i)=0; %GW
        BattSizeGW(1,i)=0; %GW
        H2bufferkg(1,i)=0; %kg
        N2bufferkg(1,i)=0; %kg
        batterykWh(1,i)=0; %kWh
        NH3Cap(1,i)=0;
        LCOA(1,i)=0;
        Extramp_wind{1,i} = zeros(8760,1);
        NH3Prodramp_wind{1,i} = zeros(8760,1);
        error_index(i,1)=i;
    end
end
% store the optimized process sizing to the WP_ramp structure
WP_ramp.turbinekWp=turbinekWp;
WP_ramp.H2SizeGW=H2SizeGW;
WP_ramp.N2SizeGW=N2SizeGW;
WP_ramp.BattSizeGW=BattSizeGW;
WP_ramp.H2bufferkg=H2bufferkg;
WP_ramp.N2bufferkg=N2bufferkg;
WP_ramp.batterykWh=batterykWh;
WP_ramp.NH3Cap = NH3Cap;
WP_ramp.LCOA=LCOA;

% 3. Find the cost of each of the process components with optimizal sizing and sort results according to location
% this code follows the same procedure as Wind_Data_Processing.mat
WP_ramp.curtailment=(WP_ramp.turbinekWp-WP_nocurt.turbinekWp)./WP_ramp.turbinekWp; % get the level of curtailment for each iteration of the solver
for k = 1:length(WP_ramp.turbinekWp) % create arrays for turbine capital and O&M cost
    if WP_ramp.category(k)=="onshore"
        turbineCost(k,1) = annual * WP_ramp.turbinekWp(k) * WP_base.turbineCapital;
        turbineOnM(k,1) = turbineCost(k,1) * WP_base.turbine_OnM_perc;
    elseif WP_ramp.category(k)=="offshore" 
        turbineCost(k,1) = annual * WP_ramp.turbinekWp(k) * WP_base.turbineCapital_offshore;
        turbineOnM(k,1) =turbineCost(k,1) * WP_base.turbine_OnM_offshore_perc;
    elseif WP_ramp.category(k)=="floating"
        turbineCost(k,1) = annual * WP_ramp.turbinekWp(k) * WP_base.turbineCapital_floating;
        turbineOnM(k,1) = turbineCost(k,1) * WP_base.turbine_OnM_floating_perc;
    end
end
len=length(WP_ramp.location);
electrolyserCost=annual * WP_ramp.H2SizeGW' * electrolyserCapital; % electrolyser capital
H2compCost=annual * WP_ramp.H2SizeGW' * H2compCapital; % H2 compressor capital
ASUCost=annual * WP_ramp.N2SizeGW' * ASUCapital; % PSA ASU capital
N2compCost=annual * WP_ramp.N2SizeGW' * N2compCapital; % N2 compressor capital
H2storeCost=annual * WP_ramp.H2bufferkg' * H2storageCapital; % H2 storage capital
N2storeCost=annual * WP_ramp.N2bufferkg' * N2storageCapital; % N2 storage capital
BattstoreCost=annual * (WP_ramp.batterykWh' * BattstorageCapital + WP_ramp.BattSizeGW' * BattpowerCapital); % battery capital
HBCost=(WP_ramp.NH3Cap'/PReq) * annual * HBCapital; % HB capital
electOnM=elect_OnM * (WP_ramp.H2SizeGW' * electrolyserCapital + WP_ramp.H2SizeGW' * H2compCapital); % H2 production O&M
ASUOnM=HB_OnM * (WP_ramp.N2SizeGW' *ASUCapital + WP_ramp.N2SizeGW' * N2compCapital); % N2 production O&M
HBOnM=(WP_ramp.NH3Cap'/PReq) * HB_OnM * HBCapital; % HB O&M
WP_ramp.totalCapital = (turbineCost + electrolyserCost + H2compCost + ASUCost + N2compCost + N2storeCost + H2storeCost + BattstoreCost + HBCost); % total capital in Dollars
WP_ramp.LCOA = (WP_ramp.totalCapital + electOnM + turbineOnM + HBOnM +ASUOnM)/(230*365); % LCOA in Dollars/tonne
locs=WP_ramp.location';
WP_ramp.costBreakdown=table(locs,turbineCost,electrolyserCost,H2compCost,ASUCost,N2compCost,H2storeCost,N2storeCost,BattstoreCost,HBCost,turbineOnM,electOnM,ASUOnM,HBOnM); % cost components stored in table

% sort results of LCOA, fraction curtailment, total energy, turbine
% size, and HB plant capacity by location
LCOA_byLoc=[];
curtailment_byLoc=[];
energy_byLoc=[];
turbinekWp_byLoc=[];
NH3Cap_byLoc=[];
unique_loc=unique(WP_ramp.location,'stable');
for i=1:length(unique_loc)
    for ii=1:length(WP_ramp.location)
        if unique_loc(i)==WP_ramp.location(ii)
            if WP_ramp.year{ii}==2013
                LCOA_byLoc(i,1)=WP_ramp.LCOA(ii);
                curtailment_byLoc(i,1)=WP_ramp.curtailment(ii);
                energy_byLoc(i,1)=sum(WP_ramp.turbinekWp(ii) * WP_base.turbinePower{ii}/1000);
                turbinekWp_byLoc(i,1)=WP_ramp.turbinekWp(ii);
                NH3Cap_byLoc(i,1)=WP_ramp.NH3Cap(ii);
            elseif WP_ramp.year{ii}==2014
                LCOA_byLoc(i,2)=WP_ramp.LCOA(ii);
                curtailment_byLoc(i,2)=WP_ramp.curtailment(ii);
                energy_byLoc(i,2)=sum(WP_ramp.turbinekWp(ii) * WP_base.turbinePower{ii}/1000);
                turbinekWp_byLoc(i,2)=WP_ramp.turbinekWp(ii);
                NH3Cap_byLoc(i,2)=WP_ramp.NH3Cap(ii);
            elseif WP_ramp.year{ii}==2015
                LCOA_byLoc(i,3)=WP_ramp.LCOA(ii);
                curtailment_byLoc(i,3)=WP_ramp.curtailment(ii);
                energy_byLoc(i,3)=sum(WP_ramp.turbinekWp(ii) * WP_base.turbinePower{ii}/1000);
                turbinekWp_byLoc(i,3)=WP_ramp.turbinekWp(ii);
                NH3Cap_byLoc(i,3)=WP_ramp.NH3Cap(ii);
            elseif WP_ramp.year{ii}==2016
                LCOA_byLoc(i,4)=WP_ramp.LCOA(ii);
                curtailment_byLoc(i,4)=WP_ramp.curtailment(ii);
                energy_byLoc(i,4)=sum(WP_ramp.turbinekWp(ii) * WP_base.turbinePower{ii}/1000);
                turbinekWp_byLoc(i,4)=WP_ramp.turbinekWp(ii);
                NH3Cap_byLoc(i,4)=WP_ramp.NH3Cap(ii);
            elseif WP_ramp.year{ii}==2017
                LCOA_byLoc(i,5)=WP_ramp.LCOA(ii);
                curtailment_byLoc(i,5)=WP_ramp.curtailment(ii);
                energy_byLoc(i,5)=sum(WP_ramp.turbinekWp(ii) * WP_base.turbinePower{ii}/1000);
                turbinekWp_byLoc(i,5)=WP_ramp.turbinekWp(ii);
                NH3Cap_byLoc(i,5)=WP_ramp.NH3Cap(ii);
            elseif WP_ramp.year{ii}==2018
                LCOA_byLoc(i,6)=WP_ramp.LCOA(ii);
                curtailment_byLoc(i,6)=WP_ramp.curtailment(ii);
                energy_byLoc(i,6)=sum(WP_ramp.turbinekWp(ii) * WP_base.turbinePower{ii}/1000);
                turbinekWp_byLoc(i,6)=WP_ramp.turbinekWp(ii);
                NH3Cap_byLoc(i,6)=WP_ramp.NH3Cap(ii);
            end
        end
    end
end
WP_ramp.LCOA_byLoc=LCOA_byLoc;
WP_ramp.unique_loc=unique_loc';
WP_ramp.curtailment_byLoc=curtailment_byLoc;
WP_ramp.energy_byLoc=energy_byLoc;
WP_ramp.turbinekWp_byLoc=turbinekWp_byLoc;
WP_ramp.NH3Cap_byLoc=NH3Cap_byLoc;
clearvars -except WP_ramp Extramp_wind NH3Prodramp_wind
% save WP_ramp Extramp_wind and NH3Prodramp_wind in the file WP_results
save('WP_results')