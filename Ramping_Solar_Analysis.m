% 1. Initialize the structure to save results and load the necessary data
clearvars
load('SP_base')
load('SP_nocurt')
load('SP_curt')
SP_ramp=struct; % initialize the SP_ramp structure
SP_ramp.ID=SP_nocurt.ID; % transfer info from the SP_nocurt structure
SP_ramp.filename=SP_nocurt.filename;
SP_ramp.year=SP_nocurt.year;
SP_ramp.location=SP_nocurt.location;
SP_ramp.unique_loc=SP_nocurt.unique_loc;

PReq=100; % required size of HB plant in MW (in this case it is the equivalent average production with ramping)
panelCapital = SP_base.panelCapital; % panel capital cost in $/kWp
electrolyserCapital = SP_base.electrolyserCapital; % electrolyser capital cost in $/GW fed in
ASUCapital = SP_base.ASUCapital; % PSA ASU capital cost in $/GW fed in
H2compCapital = SP_base.H2compCapital; % H2 compressor cost in $/GW fed in
N2compCapital = SP_base.N2compCapital; % N2 compressor cost in $/GW fed in
HBCapital = SP_base.HBCapital; % HB capital cost in dollars for 230 t/day plant
H2storageCapital = SP_base.H2storageCapital; % H2 storage capital cost $/kg
N2storageCapital = SP_base.N2storageCapital; % N2 storage capital cost $/kg
BattstorageCapital = SP_base.BattstorageCapital; % battery storage energy capacity cost $/kWh
BattpowerCapital = SP_base.BattpowerCapital; % battery storage power capacity cost $/GW fed in
panel_OnM = SP_base.panel_OnM;  % panel O&M in $ per MW per year
elect_OnM = SP_base.elect_OnM; % electrolyser O&M in percent of total capital cost
HB_OnM = SP_base.HB_OnM; % fixed O&M for HB plant in percent of capital cost per year. The parameter is also used for O&M of compressors and ASU. 
discountRate = SP_base.discountRate; % discount rate of capital
opYear = SP_base.opYear; % operating years of the plant
annual=(SP_base.discountRate/(1-(1+SP_base.discountRate)^(-1*SP_base.opYear))); %annaulization factor of capital

error_index=zeros(length(SP_ramp.location),1); % initialize array to stay index of iterations where solver returns error
minCap = 0.6; % minimum fraction capacity of the HB process (i.e. 60% of the design capacity)
rampRate = 0.2; %maximum change in fraction capacity per 4 hours (i.e. the process can go up by 5% of design capacity per hour)
NH3Prodramp_solar = cell(1,length(SP_ramp.location)); % initialize cell array to store the arrays of HB production over time for each iteration
Extramp_solar = cell(1,length(SP_ramp.location)); % initialize cell array to store the arrays of energy extraction over time for each iteration
% initialize arrays to store the results of the solver
panelkWp=zeros(1,length(SP_ramp.location));
H2SizeGW=zeros(1,length(SP_ramp.location));%GW
N2SizeGW=zeros(1,length(SP_ramp.location)); %GW
BattSizeGW=zeros(1,length(SP_ramp.location)); %GW
H2bufferkg=zeros(1,length(SP_ramp.location)); %kg
N2bufferkg=zeros(1,length(SP_ramp.location)); %kg
batterykWh=zeros(1,length(SP_ramp.location));
NH3Cap=zeros(1,length(SP_ramp.location)); % store the installed capacity of HB plant
LCOA=zeros(1,length(SP_ramp.location));

% Set an array of indices to covert the 4 hourly timescale of HB production
% levels to the hourly timescale of the power supply
h = 0;
for v = 1:2190
    for g = 1:4
        h = h+1;
        indices(h) = v;
    end
end

% 2. Configure and solve linear problem to find the optimized process design 
% with curtailment for each of the locations/years, while allowing ramping of the HB process. 
% This uses parrallelization to speed up execution.
% if needed, the parrallel loop can be exceuted as normal for loop 
% if the solver encounters errors, they are saved in the array "error_index". These iterations would need to be repeated. 
parfor i = 1:length(SP_ramp.location)
    N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
    % To speed the solver, upper limits for the panel capacity, energy
    % extraction capacity, and energy storage capacity is set. 
    Up1 = 1.5 * SP_curt.panelkWp(i,1)/(PReq*1000); % 1.5X the panel capacity without ramping
    Up2 = 1.5 * SP_curt.H2SizeGW(i,1)*1000/(PReq); % 1.5X the energy extraction capacity without ramping
    Up3 = 1.5 * SP_curt.H2bufferkg(i,1)/((1/0.2045)); % the energy storagae capacity without ramping
    
    % The next five lines solve for the minimum solar panel capacity (i.e. without curtailment) using
    % the same method as in "Solar_Data_Processing"
    panelkWpGuess = 1e6; % kWp
    x0 = panelkWpGuess;
    inventorySum = @(x) sum((panelPower*x/1e9 - PReq/1000)*3600); % energy balance
    x = fzero(@(x) inventorySum (x), x0);
    solarMin = x/(1000*PReq); % minimum panel capacity realtive to "PReq

    solarProblem=optimproblem('ObjectiveSense','minimize'); % initialize problem
    % Create variables
    NH3cap=optimvar('NH3cap',1,'LowerBound',1,'UpperBound',1/minCap); % size of HB process relative to 230 t/day plant. Maximum design capacity is 1/minCap because the total production is still 230 t/day overall
    E_size_solar=optimvar('E_size_solar',1,'LowerBound',solarMin,'UpperBound',Up1); % size of panels relative to average HB size (100 MW)
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
    Energy = E_supply == (PReq*E_size_solar*panelPower)/1E6; % defining E supply variable in GW
    Energy_min = sum(E_supply) >= PReq*8760/1000; % minimum energy amount to close the energy balance
    powerBalance = Elec_extract + H2_extract <= E_supply; %energy to electricity and H2 production must be less than the energy supply
    Max_H2 = (H2_extract) <= H2_size*PReq/1000; % power for H2 production cannot exceed H2 production capacity
    Max_Batt = Elec_extract - (27.5/640)*NH3prod(indices(1:end))*(1/1000) <= Batt_size*PReq/1000; % power sent to the battery cannot exceent Battery power capacity

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
    Capex=annual*((E_size_solar*PReq*1000*panelCapital) ...
            +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
            +(N2_size*PReq*NH3cap*(ASUCapital+N2compCapital)/1000) ...
            +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
            +(Batt_size*PReq*BattpowerCapital/1000) ...
            +HBCapital * NH3cap);
    % define operating cost 
    Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
            +(panel_OnM*E_size_solar*PReq) ...
            +(HB_OnM*(N2_size*PReq*NH3cap*(ASUCapital+N2compCapital)/1000)) ...
            +(HBCapital*HB_OnM*NH3cap);
    solarProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne
    
    % add constraints to the problem
    solarProblem.Constraints.boundedH2 = boundedH2;
    solarProblem.Constraints.boundedH2Initial = boundedH2Initial;
    solarProblem.Constraints.boundedBatt = boundedBatt;
    solarProblem.Constraints.boundedBattInitial = boundedBattInitial;
    solarProblem.Constraints.Energy = Energy;
    solarProblem.Constraints.powerBalance = powerBalance;
    solarProblem.Constraints.Energy_min = Energy_min;
    solarProblem.Constraints.H2Balance = H2Balance;
    solarProblem.Constraints.BattBalance = BattBalance;
    solarProblem.Constraints.H2finalBalance = H2finalBalance;
    solarProblem.Constraints.BattfinalBalance = BattfinalBalance;
    solarProblem.Constraints.Max_H2 = Max_H2;
    solarProblem.Constraints.Max_Batt = Max_Batt;
    solarProblem.Constraints.NH3Change = NH3change;
    solarProblem.Constraints.RampLow = RampLow;
    solarProblem.Constraints.RampHigh = RampHigh;
    solarProblem.Constraints.MaxCap = MaxCap;
    solarProblem.Constraints.MinCap = MinCap;
    solarProblem.Constraints.ProdTotal = ProdTotal;
    [sol,fval,exitflag,output] = solve(solarProblem); % execute solver
    if exitflag==1 % if the solver is successful, add the optimized process component sizes to arrays
        panelkWp(1,i)=sol.E_size_solar*PReq*1000;
        H2SizeGW(1,i)=sol.H2_size*PReq/1000;%GW
        N2SizeGW(1,i)=N2_size*PReq*sol.NH3cap/1000; %GW
        BattSizeGW(1,i)=sol.Batt_size*PReq/1000; %GW
        H2bufferkg(1,i)=sol.H2_buffer*(1/0.2045); %kg
        N2bufferkg(1,i)=0; %kg
        batterykWh(1,i)=sol.Batt_buffer*277.78;
        NH3Cap(1,i)=PReq*sol.NH3cap; % store the installed capacity of HB plant
        LCOA(1,i)=fval;
        NH3Prodramp_solar{1,i} = sol.NH3prod; % store the HB production over time in cell array
        Extramp_solar{1,i} = sol.H2_extract + sol.Elec_extract;
    else % else, add 0 to process component size arrays and add index to error_index array
        panelkWp(1,i)=0;
        H2SizeGW(1,i)=0;%GW
        N2SizeGW(1,i)=0; %GW
        BattSizeGW(1,i)=0; %GW
        H2bufferkg(1,i)=0; %kg
        N2bufferkg(1,i)=0; %kg
        batterykWh(1,i)=0; %kWh
        NH3Cap(1,i)=0;
        LCOA(1,i)=0;
        NH3Prodramp_solar{1,i} = zeros(8760,1);
        Extramp_solar{1,i} = zeros(8760,1);
        error_index(i,1)=i;
    end
end
% store the optimized process sizing to the SP_ramp structure
SP_ramp.panelkWp=panelkWp;
SP_ramp.H2SizeGW=H2SizeGW;
SP_ramp.N2SizeGW=N2SizeGW;
SP_ramp.BattSizeGW=BattSizeGW;
SP_ramp.H2bufferkg=H2bufferkg;
SP_ramp.N2bufferkg=N2bufferkg;
SP_ramp.batterykWh=batterykWh;
SP_ramp.NH3Cap = NH3Cap;
SP_ramp.LCOA=LCOA;

% 3. Find the cost of each of the process components with optimizal sizing and sort results according to location
% this code follows the same procedure as in Solar_Data_Processing.mat
SP_ramp.curtailment=(SP_ramp.panelkWp-SP_nocurt.panelkWp)./SP_ramp.panelkWp; % get the level of curtailment for each iteration of the solver
len=length(SP_ramp.location);
panelCost=annual * SP_ramp.panelkWp' * panelCapital; % panel capital
electrolyserCost=annual * SP_ramp.H2SizeGW' * electrolyserCapital; % electrolyser capital
H2compCost=annual * SP_ramp.H2SizeGW' * H2compCapital; % H2 compressor capital
ASUCost=annual * SP_ramp.N2SizeGW' * ASUCapital; % PSA ASU capital
N2compCost=annual * SP_ramp.N2SizeGW' * N2compCapital; % N2 compressor capital
H2storeCost=annual * SP_ramp.H2bufferkg' * H2storageCapital; % H2 storeage capital
N2storeCost=annual * SP_ramp.N2bufferkg' * N2storageCapital; % N2 storage capital
BattstoreCost=annual * (SP_ramp.batterykWh' * BattstorageCapital + SP_ramp.BattSizeGW' * BattpowerCapital); % battery capital
HBCost=(SP_ramp.NH3Cap'/PReq) * annual * HBCapital; % HB capital
panelOnM=SP_ramp.panelkWp' * (1/1000) * panel_OnM; % panel O&M
electOnM=elect_OnM * (SP_ramp.H2SizeGW' * electrolyserCapital + SP_ramp.H2SizeGW' * H2compCapital); % H2 production O&M
ASUOnM=HB_OnM * (SP_ramp.N2SizeGW' *ASUCapital + SP_ramp.N2SizeGW' * N2compCapital); % N2 production O&M
HBOnM=(SP_ramp.NH3Cap'/PReq) * HB_OnM * HBCapital; % HB O&M
SP_ramp.totalCapital = (panelCost + electrolyserCost + H2compCost + ASUCost + N2compCost + N2storeCost + H2storeCost + BattstoreCost + HBCost); % total capital cost in Dollars
SP_ramp.LCOA = (SP_ramp.totalCapital + electOnM + panelOnM + HBOnM +ASUOnM)/(230*365); % LCOA in Dollars/tonne
locs=SP_ramp.location';
SP_ramp.costBreakdown=table(locs,panelCost,electrolyserCost,H2compCost,ASUCost,N2compCost,H2storeCost,N2storeCost,BattstoreCost,HBCost,panelOnM,electOnM,ASUOnM,HBOnM); % cost components stored in table

% sort results of LCOA, fraction curtailment, total energy, solar panel
% size, and HB plant capacity by location
LCOA_byLoc=[];
curtailment_byLoc=[];
energy_byLoc=[];
panelkWp_byLoc=[];
NH3Cap_byLoc=[];
unique_loc=unique(SP_ramp.location,'stable');
for i=1:length(unique_loc)
    for ii=1:length(SP_ramp.location)
        if unique_loc(i)==SP_ramp.location(ii)
            if SP_ramp.year{ii}==2013
                LCOA_byLoc(i,1)=SP_ramp.LCOA(ii);
                curtailment_byLoc(i,1)=SP_ramp.curtailment(ii);
                energy_byLoc(i,1)=sum(SP_ramp.panelkWp(ii) * SP_base.panelPower{ii}/1000);
                panelkWp_byLoc(i,1)=SP_ramp.panelkWp(ii);
                NH3Cap_byLoc(i,1)=SP_ramp.NH3Cap(ii);
            elseif SP_ramp.year{ii}==2014
                LCOA_byLoc(i,2)=SP_ramp.LCOA(ii);
                curtailment_byLoc(i,2)=SP_ramp.curtailment(ii);
                energy_byLoc(i,2)=sum(SP_ramp.panelkWp(ii) * SP_base.panelPower{ii}/1000);
                panelkWp_byLoc(i,2)=SP_ramp.panelkWp(ii);
                NH3Cap_byLoc(i,2)=SP_ramp.NH3Cap(ii);
            elseif SP_ramp.year{ii}==2015
                LCOA_byLoc(i,3)=SP_ramp.LCOA(ii);
                curtailment_byLoc(i,3)=SP_ramp.curtailment(ii);
                energy_byLoc(i,3)=sum(SP_ramp.panelkWp(ii) * SP_base.panelPower{ii}/1000);
                panelkWp_byLoc(i,3)=SP_ramp.panelkWp(ii);
                NH3Cap_byLoc(i,3)=SP_ramp.NH3Cap(ii);
            elseif SP_ramp.year{ii}==2016
                LCOA_byLoc(i,4)=SP_ramp.LCOA(ii);
                curtailment_byLoc(i,4)=SP_ramp.curtailment(ii);
                energy_byLoc(i,4)=sum(SP_ramp.panelkWp(ii) * SP_base.panelPower{ii}/1000);
                panelkWp_byLoc(i,4)=SP_ramp.panelkWp(ii);
                NH3Cap_byLoc(i,4)=SP_ramp.NH3Cap(ii);
            elseif SP_ramp.year{ii}==2017
                LCOA_byLoc(i,5)=SP_ramp.LCOA(ii);
                curtailment_byLoc(i,5)=SP_ramp.curtailment(ii);
                energy_byLoc(i,5)=sum(SP_ramp.panelkWp(ii) * SP_base.panelPower{ii}/1000);
                panelkWp_byLoc(i,5)=SP_ramp.panelkWp(ii);
                NH3Cap_byLoc(i,5)=SP_ramp.NH3Cap(ii);
            elseif SP_ramp.year{ii}==2018
                LCOA_byLoc(i,6)=SP_ramp.LCOA(ii);
                curtailment_byLoc(i,6)=SP_ramp.curtailment(ii);
                energy_byLoc(i,6)=sum(SP_ramp.panelkWp(ii) * SP_base.panelPower{ii}/1000);
                panelkWp_byLoc(i,6)=SP_ramp.panelkWp(ii);
                NH3Cap_byLoc(i,6)=SP_ramp.NH3Cap(ii);
            end
        end
    end
end
SP_ramp.LCOA_byLoc=LCOA_byLoc;
SP_ramp.unique_loc=unique_loc';
SP_ramp.curtailment_byLoc=curtailment_byLoc;
SP_ramp.energy_byLoc=energy_byLoc;
SP_ramp.panelkWp_byLoc=panelkWp_byLoc;
SP_ramp.NH3Cap_byLoc=NH3Cap_byLoc;
clearvars -except SP_ramp Extramp_solar NH3Prodramp_solar
% save SP_ramp Extramp_solar and NH3Prodramp_solar in the file SP_ramp_results
save('SP_ramp_results')