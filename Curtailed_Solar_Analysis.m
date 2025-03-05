% 1. Load necessary structures, initialize structure SP_curt, 
% and extract cost parameters from the SP_base structure in variables. 
% (this is done so that the linear solver does not call SP_base each iteration)

clearvars
load('SP_base')
load('SP_nocurt')
SP_curt=struct; % initialize SP_curt structure
SP_curt.ID=SP_nocurt.ID; % transfer some info from SP_nocurt
SP_curt.filename=SP_nocurt.filename;
SP_curt.year=SP_nocurt.year;
SP_curt.location=SP_nocurt.location;
SP_curt.unique_loc=SP_nocurt.unique_loc;

PReq=100; % required size of HB plant in MW
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
panel_OnM = SP_base.panel_OnM; % panel O&M in $ per MW per year
elect_OnM = SP_base.elect_OnM; % electrolyser O&M in percent of total capital cost
HB_OnM = SP_base.HB_OnM; % fixed O&M for HB plant in percent of capital cost per year. The parameter is also used for O&M of compressors and ASU. 
discountRate = SP_base.discountRate; % discount rate of capital
opYear = SP_base.opYear; % operating years of the plant
annual=(SP_base.discountRate/(1-(1+SP_base.discountRate)^(-1*SP_base.opYear))); % annualization of the capital cost

% 2. Configure and solve linear problem to find the optimized process design 
% with curtailment for each of the locations/years. This uses parrallelization 
% to speed up execution.
% if needed, the parrallel loop can be exceuted as normal for loop 
% if the solver encounters errors, they are saved in the array "error_index". 
% These iterations would need to be repeated. 
error_index=zeros(length(SP_base.location),1); %initialize array which will store the indeces of any iterations that encounter an error
parfor i=1:length(SP_curt.filename) % iterate through even entry in SP_curt using parrallel for loop
    panelPower=SP_base.panelPower{i}; % extract the panel power in W/kWp
    N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
    solarProblem=optimproblem('ObjectiveSense','minimize'); % initialize optimization problem
    % Create variables
    E_size_solar=optimvar('E_size_solar',1,'LowerBound',0); % size of panels realtive to HB size (100 MW)
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
    Energy = E_supply == (PReq*E_size_solar*panelPower)/1E6; % defining E supply variable in GW
    Energy_min = sum(E_supply) >= PReq*8760/1000; % minimum energy amount to close the energy balance
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
    Capex=annual*((E_size_solar*PReq*1000*panelCapital) ...
            +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
            +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
            +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
            +(Batt_size*PReq*BattpowerCapital/1000) ...
            +HBCapital);
    % define operating cost 
    Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
            +(panel_OnM*E_size_solar*PReq) ...
            +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
            +(HBCapital*HB_OnM);
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
    [sol,fval,exitflag,output] = solve(solarProblem); % execute solver
    if exitflag==1 % if the solver is successful, add the optimized process component sizes to arrays
        panelkWp(1,i)=sol.E_size_solar*PReq*1000;
        H2SizeGW(1,i)=sol.H2_size*PReq/1000;%GW
        N2SizeGW(1,i)=N2_size*PReq/1000; %GW
        BattSizeGW(1,i)=sol.Batt_size*PReq/1000; %GW
        H2bufferkg(1,i)=sol.H2_buffer*(1/0.2045); %kg
        N2bufferkg(1,i)=0; %kg
        batterykWh(1,i)=sol.Batt_buffer*277.78;
        LCOA(1,i)=fval;
    else % else, add 0 to process component size arrays and add index to error_index array
        panelkWp(1,i)=0;
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
% store the optimized process sizing to the SP_curt structure
SP_curt.panelkWp=panelkWp;
SP_curt.H2SizeGW=H2SizeGW;
SP_curt.N2SizeGW=N2SizeGW;
SP_curt.BattSizeGW=BattSizeGW;
SP_curt.H2bufferkg=H2bufferkg;
SP_curt.N2bufferkg=N2bufferkg;
SP_curt.batterykWh=batterykWh;
SP_curt.LCOA=LCOA;
SP_curt.curtailment=(SP_curt.panelkWp-SP_nocurt.panelkWp)./SP_curt.panelkWp; % find the fraction of curtailment based on the solar panel size in SP_nocurt
save('SP_curt','SP_curt') % save SP_curt structure

% 3. Find the cost of each of the process components with optimizal sizing and sort results according to location
% this code follows the same procedure as in Solar_Data_Processing.mat
len=length(SP_curt.location);
panelCost=annual * SP_curt.panelkWp' * panelCapital; % panel capital cost
electrolyserCost=annual * SP_curt.H2SizeGW' * electrolyserCapital; % electrolyser capial cost
H2compCost=annual * SP_curt.H2SizeGW' * H2compCapital; % H2 compressor cappital
ASUCost=annual * SP_curt.N2SizeGW' * ASUCapital; % ASU capital
N2compCost=annual * SP_curt.N2SizeGW' * N2compCapital; % N2 compressor capital
H2storeCost=annual * SP_curt.H2bufferkg' * H2storageCapital; % H2 storage capital
N2storeCost=annual * SP_curt.N2bufferkg' * N2storageCapital; % N2 storage capital
BattstoreCost=annual * (SP_curt.batterykWh' * BattstorageCapital + SP_curt.BattSizeGW' * BattpowerCapital); % battery capital
HBCost=ones(len,1)* annual * HBCapital; % HB capital
panelOnM=SP_curt.panelkWp' * (1/1000) * panel_OnM; % panel O&M
electOnM=elect_OnM * (SP_curt.H2SizeGW' * electrolyserCapital + SP_curt.H2SizeGW' * H2compCapital); % H2 production O&M
ASUOnM=HB_OnM * (SP_curt.N2SizeGW' *ASUCapital + SP_curt.N2SizeGW' * N2compCapital); % N2 production O&M
HBOnM=ones(len,1)* HB_OnM * HBCapital; % HB O&M
SP_curt.totalCapital = (panelCost + electrolyserCost + H2compCost + ASUCost + N2compCost + N2storeCost + H2storeCost + BattstoreCost + HBCost); % total capital cost in Dollars
SP_curt.LCOA = (SP_curt.totalCapital + electOnM + panelOnM + HBOnM +ASUOnM)/(230*365); % LCOA in Dollars/tonne
locs=SP_curt.location';
SP_curt.costBreakdown=table(locs,panelCost,electrolyserCost,H2compCost,ASUCost,N2compCost,H2storeCost,N2storeCost,BattstoreCost,HBCost,panelOnM,electOnM,ASUOnM,HBOnM); % cost of components store in a table

% sort results of LCOA, fraction curtailment, total energy, and solar panel size
% by location
LCOA_byLoc=[];
curtailment_byLoc=[];
energy_byLoc=[];
panelkWp_byLoc=[];
unique_loc=unique(SP_curt.location,'stable');
for i=1:length(unique_loc)
    for ii=1:length(SP_curt.location)
        if unique_loc(i)==SP_curt.location(ii)
            if SP_curt.year{ii}==2013
                LCOA_byLoc(i,1)=SP_curt.LCOA(ii);
                curtailment_byLoc(i,1)=SP_curt.curtailment(ii);
                energy_byLoc(i,1)=sum(SP_curt.panelkWp(ii) * SP_base.panelPower{ii}/1000);
                panelkWp_byLoc(i,1)=SP_curt.panelkWp(ii);
            elseif SP_curt.year{ii}==2014
                LCOA_byLoc(i,2)=SP_curt.LCOA(ii);
                curtailment_byLoc(i,2)=SP_curt.curtailment(ii);
                energy_byLoc(i,2)=sum(SP_curt.panelkWp(ii) * SP_base.panelPower{ii}/1000);
                panelkWp_byLoc(i,2)=SP_curt.panelkWp(ii);
            elseif SP_curt.year{ii}==2015
                LCOA_byLoc(i,3)=SP_curt.LCOA(ii);
                curtailment_byLoc(i,3)=SP_curt.curtailment(ii);
                energy_byLoc(i,3)=sum(SP_curt.panelkWp(ii) * SP_base.panelPower{ii}/1000);
                panelkWp_byLoc(i,3)=SP_curt.panelkWp(ii);
            elseif SP_curt.year{ii}==2016
                LCOA_byLoc(i,4)=SP_curt.LCOA(ii);
                curtailment_byLoc(i,4)=SP_curt.curtailment(ii);
                energy_byLoc(i,4)=sum(SP_curt.panelkWp(ii) * SP_base.panelPower{ii}/1000);
                panelkWp_byLoc(i,4)=SP_curt.panelkWp(ii);
            elseif SP_curt.year{ii}==2017
                LCOA_byLoc(i,5)=SP_curt.LCOA(ii);
                curtailment_byLoc(i,5)=SP_curt.curtailment(ii);
                energy_byLoc(i,5)=sum(SP_curt.panelkWp(ii) * SP_base.panelPower{ii}/1000);
                panelkWp_byLoc(i,5)=SP_curt.panelkWp(ii);
            elseif SP_curt.year{ii}==2018
                LCOA_byLoc(i,6)=SP_curt.LCOA(ii);
                curtailment_byLoc(i,6)=SP_curt.curtailment(ii);
                energy_byLoc(i,6)=sum(SP_curt.panelkWp(ii) * SP_base.panelPower{ii}/1000);
                panelkWp_byLoc(i,6)=SP_curt.panelkWp(ii);
            end
        end
    end
end
SP_curt.LCOA_byLoc=LCOA_byLoc;
SP_curt.unique_loc=unique_loc';
SP_curt.curtailment_byLoc=curtailment_byLoc;
SP_curt.energy_byLoc=energy_byLoc;
SP_curt.panelkWp_byLoc=panelkWp_byLoc;
save('SP_curt','SP_curt')