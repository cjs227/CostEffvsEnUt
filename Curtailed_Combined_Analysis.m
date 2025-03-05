% 1. Load the necessary data, initialize the combined data structures and
% input the necessary parameters
load('SP_base')
load('WP_base')
CP_base_solar=struct;
CP_base_wind=struct;
%Align the solar and wind locations by looking at the similarity in
%latitude and longitude
iii=0;
for i=1:length(SP_base.location)
    for ii=1:length(WP_base.location)
        if abs(SP_base.Lat(i)-WP_base.LatACT(ii))<0.1 && abs(SP_base.Lon(i)-WP_base.LonACT(ii))<0.1 && SP_base.year{i}==WP_base.year{ii} && WP_base.elevation(ii)>=-10
            iii=iii+1;
            CP_base_solar.location(iii)=SP_base.location(i);
            CP_base_wind.location(iii)=SP_base.location(i);
            CP_base_wind.year(iii)=WP_base.year(ii);
            CP_base_solar.year(iii)=WP_base.year(ii);
            CP_base_wind.Lat(iii)=SP_base.Lat(i);
            CP_base_solar.Lat(iii)=SP_base.Lat(i);
            CP_base_wind.Lon(iii)=SP_base.Lon(i);
            CP_base_solar.Lon(iii)=SP_base.Lon(i);
            CP_base_wind.filename(iii)=WP_base.filename(ii);
            CP_base_solar.filename(iii)=SP_base.filename(i);
            CP_base_wind.turbinePower(iii)=WP_base.turbinePower(ii);
            CP_base_solar.panelPower(iii)=SP_base.panelPower(i);
            CP_base_wind.elevation(iii)=WP_base.elevation(ii);
            break
        end
    end
end
% Tranfer the necessary parameters to the new structures
CP_base_solar.panelCapital=SP_base.panelCapital;
CP_base_solar.electrolyserCapital=SP_base.electrolyserCapital;
CP_base_solar.H2compCapital=SP_base.H2compCapital;
CP_base_solar.ASUCapital=SP_base.ASUCapital;
CP_base_solar.N2compCapital=SP_base.N2compCapital;
CP_base_solar.H2storageCapital=SP_base.H2storageCapital;
CP_base_solar.N2storageCapital=SP_base.N2storageCapital;
CP_base_solar.BattstorageCapital=SP_base.BattstorageCapital;
CP_base_solar.BattpowerCapital=SP_base.BattpowerCapital;
CP_base_solar.HBCapital=SP_base.HBCapital;
CP_base_solar.panel_OnM=SP_base.panel_OnM;
CP_base_solar.elect_OnM=SP_base.elect_OnM;
CP_base_solar.HB_OnM=SP_base.HB_OnM;
CP_base_solar.discountRate=SP_base.discountRate;
CP_base_solar.opYear=SP_base.opYear;
CP_base_wind.turbineCapital=WP_base.turbineCapital;
CP_base_wind.turbineCapital_offshore=WP_base.turbineCapital_offshore;
CP_base_wind.turbineCapital_floating=WP_base.turbineCapital_floating;
CP_base_wind.turbine_OnM=WP_base.turbine_OnM;
CP_base_wind.turbine_OnM_offshore=WP_base.turbine_OnM_offshore;
CP_base_wind.turbine_OnM_floating=WP_base.turbine_OnM_floating;
CP_base_wind.turbine_OnM_perc=WP_base.turbine_OnM_perc;
CP_base_wind.turbine_OnM_offshore_perc=WP_base.turbine_OnM_offshore_perc;
CP_base_wind.turbine_OnM_floating_perc=WP_base.turbine_OnM_floating_perc;
save('CP_base_wind','CP_base_wind')
save('CP_base_solar','CP_base_solar')
clearvars -except CP_base_solar CP_base_wind

% 2. Initialize the run the optimization algorithm for the case of 
% curtailment
CP_curt=struct;
CP_curt.windFilename=CP_base_wind.filename;
CP_curt.solarFilename=CP_base_solar.filename;
CP_curt.year=CP_base_solar.year;
CP_curt.location=CP_base_solar.location;
% open the necessary parameters
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

% Initialize the arrays for saving the optimization results
panelkWp=zeros(1,length(CP_curt.location));
turbinekWp = zeros(1,length(CP_curt.location));
H2SizeGW=zeros(1,length(CP_curt.location));%GW
N2SizeGW=zeros(1,length(CP_curt.location)); %GW
BattSizeGW=zeros(1,length(CP_curt.location)); %GW
H2bufferkg=zeros(1,length(CP_curt.location)); %kg
N2bufferkg=zeros(1,length(CP_curt.location)); %kg
batterykWh=zeros(1,length(CP_curt.location));
LCOA=zeros(1,length(CP_curt.location));
solar_fractInst = zeros(1,length(CP_curt.location));
% Run the optimization in a parrallel for loop
error_index=[];
parfor i=1:length(CP_curt.location)
    panelPower=sub_panelPower{i};
    turbinePower=sub_turbinePower{i};
    N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
    % initialize the optimization variables
    combinedProblem=optimproblem('ObjectiveSense','minimize');
    E_size_solar=optimvar('E_size_solar',1,'LowerBound',0); % size of panels realtive to HB size (100 MW)
    E_size_wind = optimvar('E_size_wind',1,'LowerBound',0);
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
    % Add the constraints
    Energy = E_supply == (PReq*E_size_solar*panelPower+PReq*E_size_wind*turbinePower)/1E6; %GW
    Energy_min = sum(E_supply) >= PReq*8760/1000; %minimum energy amount
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
            +(E_size_wind*PReq*1000*turbineCapital) ...
            +(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
            +(N2_size*PReq*(ASUCapital+N2compCapital)/1000) ...
            +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
            +(Batt_size*PReq*BattpowerCapital/1000) ...
            +HBCapital);
    % define operating cost 
    Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
            +(panel_OnM*E_size_solar*PReq) ...
            +(turbine_OnM*E_size_wind*PReq*1000*turbineCapital*annual) ...
            +(HB_OnM*(N2_size*PReq*(ASUCapital+N2compCapital)/1000)) ...
            +(HBCapital*HB_OnM);
    combinedProblem.Objective =(Capex+Opex)/(230*365); % objective function, LCOA, dollars per tonne

    combinedProblem.Constraints.boundedH2 = boundedH2;
    combinedProblem.Constraints.boundedH2Initial = boundedH2Initial;
    combinedProblem.Constraints.boundedBatt = boundedBatt;
    combinedProblem.Constraints.boundedBattInitial = boundedBattInitial;
    combinedProblem.Constraints.Energy = Energy;
    combinedProblem.Constraints.powerBalance = powerBalance;
    combinedProblem.Constraints.Energy_min = Energy_min;
    combinedProblem.Constraints.H2Balance = H2Balance;
    combinedProblem.Constraints.BattBalance = BattBalance;
    combinedProblem.Constraints.H2finalBalance = H2finalBalance;
    combinedProblem.Constraints.BattfinalBalance = BattfinalBalance;
    combinedProblem.Constraints.Max_H2 = Max_H2;
    combinedProblem.Constraints.Max_Batt = Max_Batt;
    [sol,fval,exitflag,output] = solve(combinedProblem);
    if exitflag==1 
        panelkWp(1,i)=sol.E_size_solar*PReq*1000;
        turbinekWp(1,i)=sol.E_size_wind*PReq*1000;
        H2SizeGW(1,i)=sol.H2_size*PReq/1000;%GW
        N2SizeGW(1,i)=N2_size*PReq/1000; %GW
        BattSizeGW(1,i)=sol.Batt_size*PReq/1000; %GW
        H2bufferkg(1,i)=sol.H2_buffer*(1/0.2045); %kg
        N2bufferkg(1,i)=0; %kg
        batterykWh(1,i)=sol.Batt_buffer*277.78;
        LCOA(1,i)=fval;
        solar_fractInst(1,i)=sol.E_size_solar/(sol.E_size_solar+sol.E_size_wind);
    else
        panelkWp(1,i)=0;
        turbinekWp(1,i)=0;
        H2SizeGW(1,i)=0;%GW
        N2SizeGW(1,i)=0; %GW
        BattSizeGW(1,i)=0; %GW
        H2bufferkg(1,i)=0; %kg
        N2bufferkg(1,i)=0; %kg
        batterykWh(1,i)=0; %kWh
        LCOA(1,i)=0;
        solar_fractInst(1,i)=0;
        error_index(1,i)=i;
    end
end
CP_curt.panelkWp=panelkWp;
CP_curt.turbinekWp=turbinekWp;
CP_curt.H2SizeGW=H2SizeGW;
CP_curt.N2SizeGW=N2SizeGW;
CP_curt.BattSizeGW=BattSizeGW;
CP_curt.H2bufferkg=H2bufferkg;
CP_curt.N2bufferkg=N2bufferkg;
CP_curt.batterykWh=batterykWh;
CP_curt.LCOA=LCOA;
CP_curt.solar_fractInst=solar_fractInst;

% 3. Use the optimization results to calculate the cost of each component
% at each location for each year
% find the total solar and wind energy, and then the total curtailment
for i =1:length(CP_curt.location)
    totalSolar=sum(CP_curt.panelkWp(i)*CP_base_solar.panelPower{i}/1000); %kWh
    totalWind=sum(CP_curt.turbinekWp(i) * CP_base_wind.turbinePower{i}/1000); %kWh
    CP_curt.curtailment(1,i)=(totalSolar+totalWind-(100000*8760))/(totalSolar+totalWind);
    CP_curt.solar_fractEnergy(1,i)=totalSolar/(totalSolar+totalWind);
end
% find the process component costs
panelCost=annual * CP_curt.panelkWp' * panelCapital;
turbineCost=annual * CP_curt.turbinekWp' * turbineCapital;
electrolyserCost=annual * CP_curt.H2SizeGW' * electrolyserCapital;
H2compCost=annual * CP_curt.H2SizeGW' * H2compCapital;
ASUCost=annual * CP_curt.N2SizeGW' * ASUCapital;
N2compCost=annual * CP_curt.N2SizeGW' * N2compCapital;
H2storeCost=annual * CP_curt.H2bufferkg' * H2storageCapital;
N2storeCost=annual * CP_curt.N2bufferkg' * N2storageCapital;
BattstoreCost=annual * (CP_curt.batterykWh' * BattstorageCapital + CP_curt.BattSizeGW' * BattpowerCapital);
HBCost=ones(length(CP_curt.location),1)* annual * HBCapital;
panelOnM=CP_curt.panelkWp' * (1/1000) * panel_OnM;
turbineOnM=turbine_OnM * CP_curt.turbinekWp' * turbineCapital * annual;
electOnM=elect_OnM * (CP_curt.H2SizeGW' * electrolyserCapital + CP_curt.H2SizeGW' * H2compCapital);
ASUOnM=HB_OnM * (CP_curt.N2SizeGW' *ASUCapital + CP_curt.N2SizeGW' * N2compCapital);
HBOnM=ones(length(CP_curt.location),1)* HB_OnM * HBCapital;
CP_curt.totalCapital = (panelCost + turbineCost + electrolyserCost + H2compCost + ASUCost + N2compCost + N2storeCost + H2storeCost + BattstoreCost + HBCost); % Dollars
CP_curt.LCOA = (CP_curt.totalCapital + electOnM + panelOnM + turbineOnM + HBOnM +ASUOnM)/(230*365); % Dollars/tonne
locs=CP_curt.location';
CP_curt.costBreakdown=table(locs,panelCost,turbineCost,electrolyserCost,H2compCost,ASUCost,N2compCost,H2storeCost,N2storeCost,BattstoreCost,HBCost,panelOnM,turbineOnM,electOnM,ASUOnM,HBOnM);

% organize the cost and curtailment metrics by locations
LCOA_byLoc=[];
curtailment_byLoc=[];
fractSolar_byLoc=[];
unique_loc=unique(CP_curt.location,'stable');
for i=1:length(unique_loc)
    for ii=1:length(CP_curt.location)
        if unique_loc(i)==CP_curt.location(ii)
            if CP_curt.year{ii}==2013
                LCOA_byLoc(i,1)=CP_curt.LCOA(ii);
                curtailment_byLoc(i,1)=CP_curt.curtailment(ii);
                fractSolar_byLoc(i,1)=CP_curt.solar_fractEnergy(ii);
            elseif CP_curt.year{ii}==2014
                LCOA_byLoc(i,2)=CP_curt.LCOA(ii);
                curtailment_byLoc(i,2)=CP_curt.curtailment(ii);
                fractSolar_byLoc(i,2)=CP_curt.solar_fractEnergy(ii);
            elseif CP_curt.year{ii}==2015
                LCOA_byLoc(i,3)=CP_curt.LCOA(ii);
                curtailment_byLoc(i,3)=CP_curt.curtailment(ii);
                fractSolar_byLoc(i,3)=CP_curt.solar_fractEnergy(ii);
            elseif CP_curt.year{ii}==2016
                LCOA_byLoc(i,4)=CP_curt.LCOA(ii);
                curtailment_byLoc(i,4)=CP_curt.curtailment(ii);
                fractSolar_byLoc(i,4)=CP_curt.solar_fractEnergy(ii);
            elseif CP_curt.year{ii}==2017
                LCOA_byLoc(i,5)=CP_curt.LCOA(ii);
                curtailment_byLoc(i,5)=CP_curt.curtailment(ii);
                fractSolar_byLoc(i,5)=CP_curt.solar_fractEnergy(ii);
            elseif CP_curt.year{ii}==2018
                LCOA_byLoc(i,6)=CP_curt.LCOA(ii);
                curtailment_byLoc(i,6)=CP_curt.curtailment(ii);
                fractSolar_byLoc(i,6)=CP_curt.solar_fractEnergy(ii);
            end
        end
    end
end
CP_curt.LCOA_byLoc=LCOA_byLoc;
CP_curt.unique_loc=unique_loc';
CP_curt.curtailment_byLoc=curtailment_byLoc;
CP_curt.fractSolar_byLoc=fractSolar_byLoc;
save('CP_curt','CP_curt')
clearvars