% 1. Read all the csv files in "Raw_Data_Solar" and store properties in the 
% structure SP_base. Each of the six years from 2013-2018 have a separate entry.

SP_base = struct; % baseline solar power structure
fullpath = fullfile("Raw_Solar_Data","*2013_2018.csv"); 
theFiles = dir(fullpath);
load('gebco_variables')
k=0;
for i = 1:length(theFiles)
    fullFileName = fullfile(theFiles(i).folder,theFiles(i).name);
    panelPowerTemp = readtable(fullFileName);
    split_file=strsplit(string(theFiles(i).name),"_");
    Lat=str2double(split_file(1,2)); % get latitude from the filename
    Lon=str2double(split_file(1,3)); %get longitude from the filename
    [lat_diff,lat_min]=min(abs(lat_gebco-Lat)); % get nearest gebco latitude
    [lon_diff,lon_min]=min(abs(lon_gebco-Lon)); % get nearest gebco longitude
    elevation=elevation_gebco(lon_min,lat_min); % get elevation based on hearest gebco latitude and longitude
    if elevation>-10 % only process locations that are onshore. an elevation threshold of -10 is used to include locations along the coast/lowland areas
        k=k+1;
        %extract data for the year 2013
        SP_base.ID(6*k-5)=k; %ID of the entry
        SP_base.year{6*k-5} = 2013; %year of the entry
        SP_base.panelPower{6*k-5} = table2array(panelPowerTemp(1:8760,2)); %power in W/kWp for the year
        SP_base.location(6*k-5) = string(theFiles(i).name(12:24)); %location as a string
        SP_base.Lat(6*k-5)=Lat; %Latitude of the location
        SP_base.Lon(6*k-5)=Lon; %Longitude of the location
        %extract data for the year 2014
        SP_base.ID(6*k-4)=k;
        SP_base.year{6*k-4} = 2014;
        SP_base.panelPower{6*k-4} = table2array(panelPowerTemp(8761:17520,2));
        SP_base.location(6*k-4) = string(theFiles(i).name(12:24));
        SP_base.Lat(6*k-4)=Lat;
        SP_base.Lon(6*k-4)=Lon;
        %extract data for the year 2015
        SP_base.ID(6*k-3)=k;
        SP_base.year{6*k-3} = 2015;
        SP_base.panelPower{6*k-3} = table2array(panelPowerTemp(17545:26304,2)); %this year is shifted by 1 day due to misapplication of leap-years
        SP_base.location(6*k-3) = string(theFiles(i).name(12:24));
        SP_base.Lat(6*k-3)=Lat;
        SP_base.Lon(6*k-3)=Lon;
        %extract data for the year 2016
        SP_base.ID(6*k-2)=k;
        SP_base.year{6*k-2} = 2016;
        SP_base.panelPower{6*k-2} = table2array(panelPowerTemp(26305:35064,2)); %this year is shifted by 1 day due to misapplication of leap-years
        SP_base.location(6*k-2) = string(theFiles(i).name(12:24));
        SP_base.Lat(6*k-2)=Lat;
        SP_base.Lon(6*k-2)=Lon;
        %extract data for the year 2017
        SP_base.ID(6*k-1)=k;
        SP_base.year{6*k-1} = 2017;
        SP_base.panelPower{6*k-1} = table2array(panelPowerTemp(35065:43824,2));
        SP_base.location(6*k-1) = string(theFiles(i).name(12:24));
        SP_base.Lat(6*k-1)=Lat;
        SP_base.Lon(6*k-1)=Lon;
        %extract data for the year 2018
        SP_base.ID(6*k)=k;
        SP_base.year{6*k} = 2018;
        SP_base.panelPower{6*k} = table2array(panelPowerTemp(43825:52584,2));
        SP_base.location(6*k) = string(theFiles(i).name(12:24));
        SP_base.Lat(6*k)=Lat;
        SP_base.Lon(6*k)=Lon;
    end
end
clearvars -except SP_base

% 2. Due to outlier points in the power profiles which are unreasonably 
% high, each of the power profiles is searched for points with >1000 and 
% these entries are replaced with the next point. 
for i=1:length(SP_base.panelPower)
    if any(SP_base.panelPower{i}>1000)
        power=SP_base.panelPower{i};
        for j=1:height(power)
            if power(j)>1000
                power(j)=power(j+1);
            end
        end
        SP_base.panelPower{i}=power;
    end
end

% 3. Input cost parameters into the structure SP_base
% (the sources and explanation of the numbers are in the manuscript)
panelCapital = 1.1*0.47E6/1000; % capital cost of solar panels in $/kWp
electrolyserCapital = (400/409)*700/(1E-6); % capital cost of electrolysers in $/GW fed in
ASUCapital = 147000*(1/45)/(1E-6); % capital cost of PSA ASU in $/GW fed in
H2compCapital = (9/409)*6670/(1E-6); % capital cost of hydrogen compressors for HB in $/GW fed in
N2compCapital = 176400*(1/45)*(1/1E-6); %capital cost of nitrogen compressors in HB in $/GW fed in
HBCapital = 125000*230; % capital cost of HB plant in $ for 230 t/day plant
H2storageCapital = 900; % captial cost of H2 storage in $/kg
N2storageCapital = 32; % capital cost of N2 storage in $/kg
BattstorageCapital = 150*1.1; % capital cost of battery energy capacity in $/kWh
BattpowerCapital = 270*1.1/(1E-6); % capital cost of battery power capacity in $/GW fed in
panel_OnM = 1.1*11300; % fixed O&M for solar panels in $ per MW per year
elect_OnM = 0.02; % fixed O&M for electrolyser in percent of capital cost per year
HB_OnM=0.02; % fixed O&M for HB plant in percent of capital cost per year. The parameter is also used for O&M of compressors and ASU. 
discountRate = 0.07; % discount rate of capital
opYear = 25; % lifetime of the plant
SP_base.panelCapital=panelCapital;
SP_base.electrolyserCapital=electrolyserCapital;
SP_base.ASUCapital=ASUCapital;
SP_base.H2compCapital = H2compCapital;
SP_base.N2compCapital = N2compCapital;
SP_base.HBCapital=HBCapital;
SP_base.H2storageCapital=H2storageCapital;
SP_base.N2storageCapital=N2storageCapital;
SP_base.BattstorageCapital=BattstorageCapital;
SP_base.BattpowerCapital=BattpowerCapital;
SP_base.panel_OnM=panel_OnM;
SP_base.elect_OnM=elect_OnM;
SP_base.HB_OnM=HB_OnM;
SP_base.discountRate=discountRate;
SP_base.opYear=opYear;
annual=(SP_base.discountRate/(1-(1+SP_base.discountRate)^(-1*SP_base.opYear))); % annualization of the capital cost

% 4. Run an algorithm to solve for the panel size (kWp) H2 production size (GW), 
% battery power size (GW), H2 buffer size (kg) and battery size (kWh) for the case of no curtailment. 
% (i.e. close the energy balance). Store these results in stucture SP_nocurt. 
% Plant size (constant energy requirement) is 100MW. 

SP_nocurt=struct; %initialize stucture 
SP_nocurt.ID=SP_base.ID; %import some values from SP_base
SP_nocurt.filename=SP_base.filename;
SP_nocurt.year=SP_base.year;
SP_nocurt.location=SP_base.location;
PReq = 100; % MW of constant supply required for plant

parfor i = 1:length(SP_nocurt.filename) %cycle through each of the power profile in SP_base
    panelPower = SP_base.panelPower{i};% extract the panel power in W/kWp
    N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
    %find the require panel size to generate the minimum amount of energy
    panelkWpGuess = 1e6; % kWp
    x0 = panelkWpGuess;
    inventorySum = @(x) sum((panelPower*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
    x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
    panelkWp(i) = x; %store panel kWp 

    E_supply = (x*panelPower)/1E9;%energy supply in GW
    % Using the required panel size, solve for the relative size of the H2
    % production and battery power, allowing for the relative power
    % commitment to each to vary over time
    solarProblem=optimproblem('ObjectiveSense','minimize'); % initialize optimization problem
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
    Elec_extract = optimvar('Batt_extract',8760,'LowerBound',0); % amount of energy committed as electricity (both direct fed and to the battery)
    % Create constraints
    powerBalance = Elec_extract + H2_extract == E_supply; %energy to electricity and H2 production must be equal to the energy supply
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
    Capex=annual*((H2_size*PReq*(electrolyserCapital+H2compCapital)/1000) ...
            +((H2_buffer*(1/0.2045)*H2storageCapital+Batt_buffer*277.78*BattstorageCapital)) ...
            +(Batt_size*PReq*BattpowerCapital/1000));
    % define operating cost 
    Opex = elect_OnM*(H2_size*PReq*(electrolyserCapital+H2compCapital)/1000);
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
    if exitflag==1 % if the solver is successful, add the optimized process component sizes to arrays
        H2SizeGW(1,i)=sol.H2_size*PReq/1000;%GW
        N2SizeGW(1,i)=N2_size*PReq/1000; %GW
        BattSizeGW(1,i)=sol.Batt_size*PReq/1000; %GW
        H2bufferkg(1,i)=sol.H2_buffer*(1/0.2045); %kg
        N2bufferkg(1,i)=0; %kg
        batterykWh(1,i)=sol.Batt_buffer*277.78;
    else % else, add 0 to process component size arrays and add index to error_index array
        H2SizeGW(1,i)=0;%GW
        N2SizeGW(1,i)=0; %GW
        BattSizeGW(1,i)=0; %GW
        H2bufferkg(1,i)=0; %kg
        N2bufferkg(1,i)=0; %kg
        batterykWh(1,i)=0; %kWh
        error_index(i,1)=i;
    end
end
SP_nocurt.panelkWp = panelkWp;
SP_nocurt.H2SizeGW = H2SizeGW;
SP_nocurt.N2SizeGW = N2SizeGW;
SP_nocurt.BattSizeGW = BattSizeGW;
SP_nocurt.H2bufferkg = H2bufferkg;
SP_nocurt.N2bufferkg = N2bufferkg;
SP_nocurt.batterykWh = batterykWh;

% 5. Calculate the cost of each process component for the case of No Curtailed energy, 
% and store in SP_nocurt
annual=discountRate/(1-((1+discountRate)^(-opYear))); % annualization factor for capital costs
len=length(SP_nocurt.ID);
panelCost=annual * SP_nocurt.panelkWp' * panelCapital; % panel capital
electrolyserCost=annual * SP_nocurt.H2SizeGW' * electrolyserCapital; % electrolyser capital
H2compCost=annual * SP_nocurt.H2SizeGW' * H2compCapital; % H2 compressor capital
ASUCost=annual * SP_nocurt.N2SizeGW' * ASUCapital; % ASU capital
N2compCost=annual * SP_nocurt.N2SizeGW' * N2compCapital; % N2 compressor capital
H2storeCost=annual * SP_nocurt.H2bufferkg' * H2storageCapital; % H2 storage capital
N2storeCost=annual * SP_nocurt.N2bufferkg' * N2storageCapital; % N2 storage capital
BattstoreCost=annual * (SP_nocurt.batterykWh' * BattstorageCapital + SP_nocurt.BattSizeGW' * BattpowerCapital); % battery capital
HBCost=ones(len,1)* annual * HBCapital; % HB process capital
panelOnM=SP_nocurt.panelkWp' * (1/1000) * panel_OnM; % panel O&M
electOnM=elect_OnM * (SP_nocurt.H2SizeGW' * electrolyserCapital + SP_nocurt.H2SizeGW' * H2compCapital); % H2 production O&M
ASUOnM=HB_OnM * (SP_nocurt.N2SizeGW' *ASUCapital + SP_nocurt.N2SizeGW' * N2compCapital); % N2 production O&M
HBOnM=ones(len,1)* HB_OnM * HBCapital; % HB O&M
SP_nocurt.totalCapital = (panelCost + electrolyserCost + H2compCost + ASUCost + N2compCost + N2storeCost + H2storeCost + BattstoreCost + HBCost); % total annualized capital cost
SP_nocurt.LCOA = (SP_nocurt.totalCapital + electOnM + panelOnM + HBOnM)/(230*365); % levelized cost of ammonia in Dollars/tonne
locs=SP_nocurt.location';
SP_nocurt.costBreakdown=table(locs,panelCost,electrolyserCost,H2compCost,ASUCost,N2compCost,H2storeCost,N2storeCost,BattstoreCost,HBCost,panelOnM,electOnM,ASUOnM,HBOnM); %all of the cost results stored in a table
LCOA_byLoc=[]; % initialized array of LCOA sorted by location
unique_loc=unique(SP_nocurt.location,'stable'); % create array of unique locations
for i=1:length(unique_loc) % cycle through every unique location and original location array
    for ii=1:length(SP_nocurt.location)
        if unique_loc(i)==SP_nocurt.location(ii) % check if the locations are the same
            % add LCOA to array LCOA_byLoc by column according to the year
            if SP_nocurt.year{ii}==2013
                LCOA_byLoc(i,1)=SP_nocurt.LCOA(ii);
            elseif SP_nocurt.year{ii}==2014
                LCOA_byLoc(i,2)=SP_nocurt.LCOA(ii);
            elseif SP_nocurt.year{ii}==2015
                LCOA_byLoc(i,3)=SP_nocurt.LCOA(ii);
            elseif SP_nocurt.year{ii}==2016
                LCOA_byLoc(i,4)=SP_nocurt.LCOA(ii);
            elseif SP_nocurt.year{ii}==2017
                LCOA_byLoc(i,5)=SP_nocurt.LCOA(ii);
            elseif SP_nocurt.year{ii}==2018
                LCOA_byLoc(i,6)=SP_nocurt.LCOA(ii);
            end
        end
    end
end
SP_nocurt.LCOA_byLoc=LCOA_byLoc;
SP_nocurt.unique_loc=unique_loc';

save('SP_nocurt','SP_nocurt'); % save the SP_nocurt structure 
save('SP_base','SP_base'); % save the SP_base structure