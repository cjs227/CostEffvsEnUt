% 1. Read all the csv files in "Raw_Data_Wind" and store properties in the 
% structure WP_base. Each of the six years from 2013-2018 have a separate entry.
WP_base = struct; % baseline wind pwoer structure
fullpath = fullfile("Raw_Data_Wind","*.nc"); 
theFiles = dir(fullpath);
load('gebco_variables')
i=0;
for j = 1:length(theFiles) % cycle through the NC file for each location in the folder
    fullFileName = fullfile(theFiles(j).folder,theFiles(j).name);
    LatACT=ncread(fullFileName,'XLAT'); % get actual latitude from NC file
    LonACT=ncread(fullFileName,'XLON'); % get actual longitude from the NC file
    if anynan(LatACT) || anynan(LonACT) % check there is not an error
        continue
    else
        i=i+1;
    end
    LatEXP=str2double(extractBetween(fullFileName,"Lat_","_Lon")); % get expected latitude based on file name
    LonEXP=str2double(extractBetween(fullFileName,"Lon_",".nc")); % get expected longitude based on the file name
    if abs(LatEXP-LatACT)>0.5 || abs(LonEXP-LonACT)>0.5 % if the actual deviates from the expected, display the longitudes and latitudes for checking
        disp('Location Deviation') 
        disp(LatACT)
        disp(LatEXP)
        disp(LonACT)
        disp(LonEXP)
    end
    [lat_diff,lat_min]=min(abs(lat_gebco-LatACT)); % get nearest gebco latitude
    [lon_diff,lon_min]=min(abs(lon_gebco-LonACT)); % get nearest gebco longitude
    elevation=elevation_gebco(lon_min,lat_min); % get elevation based on hearest gebco latitude and longitude
    
    % determine the catagory of the wind power based on the elevation
    % If the elevation is less than -1000, the location is not included
    if elevation>-1000 && elevation<-50
        category="floating";
    elseif elevation>=-50 && elevation<0
        category="offshore";
    elseif elevation>=0
        category="onshore";
    else
        i=i-1;
        continue
    end
    wind_speed=ncread(fullFileName,'WS10'); % read the wind speed profile at a height of 10 m
    for k=1:8760 % get the indeces of the profiles for the six years on 1 hour basis. Profiles in the nc files are in 30 min intervals. 
        i1=(k-1)*2+1;
        i2=(k-1)*2+17521;
        i3=(k-1)*2+35041;
        i4=(k-1)*2+52561;
        i5=(k-1)*2+70129;
        i6=(k-1)*2+87649;
        
        % extract wind power profiles depending on whether the category is
        % onshore or offshore/floating
        if category=="onshore"
            % wind speed profiles for the 6 years at a height of 100 m (10X
            % higher)
            W1=((wind_speed(i1)+wind_speed(i1+1))/2)*((10)^0.143);
            W2=((wind_speed(i2)+wind_speed(i2+1))/2)*((10)^0.143);
            W3=((wind_speed(i3)+wind_speed(i3+1))/2)*((10)^0.143);
            W4=((wind_speed(i4)+wind_speed(i4+1))/2)*((10)^0.143);
            W5=((wind_speed(i5)+wind_speed(i5+1))/2)*((10)^0.143);
            W6=((wind_speed(i6)+wind_speed(i6+1))/2)*((10)^0.143);
            % for each of the 6 years, convert wind speed to power in W/kWp
            % according to a cut in speed of 3m/s, a rated speed of 13 m/s
            % and a cutout speed of 25 m/s.
            if W1>3 && W1<=13
                windPower_2013(k,1)=((W1^3-3^3)/(13^3-3^3))*1000;
            elseif W1>13 && W1<=25
                windPower_2013(k,1)=1000;
            else
                windPower_2013(k,1)=0;
            end
            if W2>3 && W2<=13
                windPower_2014(k,1)=((W2^3-3^3)/(13^3-3^3))*1000;
            elseif W2>13 && W2<=25
                windPower_2014(k,1)=1000;
            else
                windPower_2014(k,1)=0;
            end
            if W3>3 && W3<=13
                windPower_2015(k,1)=((W3^3-3^3)/(13^3-3^3))*1000;
            elseif W3>13 && W3<=25
                windPower_2015(k,1)=1000;
            else
                windPower_2015(k,1)=0;
            end
            if W4>3 && W4<=13
                windPower_2016(k,1)=((W4^3-3^3)/(13^3-3^3))*1000;
            elseif W4>13 && W4<=25
                windPower_2016(k,1)=1000;
            else
                windPower_2016(k,1)=0;
            end
            if W5>3 && W5<=13
                windPower_2017(k,1)=((W5^3-3^3)/(13^3-3^3))*1000;
            elseif W5>13 && W5<=25
                windPower_2017(k,1)=1000;
            else
                windPower_2017(k,1)=0;
            end
            if W6>3 && W6<=13
                windPower_2018(k,1)=((W6^3-3^3)/(13^3-3^3))*1000;
            elseif W6>13 && W6<=25
                windPower_2018(k,1)=1000;
            else
                windPower_2018(k,1)=0;
            end
        elseif category=="floating" || category=="offshore"
            % wind speed profiles for the 6 years at a height of 130 m (13X
            % higher)
            W1=((wind_speed(i1)+wind_speed(i1+1))/2)*((13)^0.143);
            W2=((wind_speed(i2)+wind_speed(i2+1))/2)*((13)^0.143);
            W3=((wind_speed(i3)+wind_speed(i3+1))/2)*((13)^0.143);
            W4=((wind_speed(i4)+wind_speed(i4+1))/2)*((13)^0.143);
            W5=((wind_speed(i5)+wind_speed(i5+1))/2)*((13)^0.143);
            W6=((wind_speed(i6)+wind_speed(i6+1))/2)*((13)^0.143);
            % for each of the 6 years, convert wind speed to power in W/kWp
            % according to a cut in speed of 3m/s, a rated speed of 13 m/s
            % and a cutout speed of 31 m/s.
            if W1>3 && W1<=13
                windPower_2013(k,1)=((W1^3-3^3)/(13^3-3^3))*1000;
            elseif W1>13 && W1<=31
                windPower_2013(k,1)=1000;
            else
                windPower_2013(k,1)=0;
            end
            if W2>3 && W2<=13
                windPower_2014(k,1)=((W2^3-3^3)/(13^3-3^3))*1000;
            elseif W2>13 && W2<=31
                windPower_2014(k,1)=1000;
            else
                windPower_2014(k,1)=0;
            end
            if W3>3 && W3<=13
                windPower_2015(k,1)=((W3^3-3^3)/(13^3-3^3))*1000;
            elseif W3>13 && W3<=31
                windPower_2015(k,1)=1000;
            else
                windPower_2015(k,1)=0;
            end
            if W4>3 && W4<=13
                windPower_2016(k,1)=((W4^3-3^3)/(13^3-3^3))*1000;
            elseif W4>13 && W4<=31
                windPower_2016(k,1)=1000;
            else
                windPower_2016(k,1)=0;
            end
            if W5>3 && W5<=13
                windPower_2017(k,1)=((W5^3-3^3)/(13^3-3^3))*1000;
            elseif W5>13 && W5<=31
                windPower_2017(k,1)=1000;
            else
                windPower_2017(k,1)=0;
            end
            if W6>3 && W6<=13
                windPower_2018(k,1)=((W6^3-3^3)/(13^3-3^3))*1000;
            elseif W6>13 && W6<=31
                windPower_2018(k,1)=1000;
            else
                windPower_2018(k,1)=0;
            end
        end

    end
    
    % Add the wind power profiles in WP_base with a new entry for each year
    % and location
    WP_base.year{6*i-5} = 2013; % year
    WP_base.turbinePower{6*i-5} = windPower_2013; % power profile in W/kWP
    WP_base.LatACT(6*i-5)=LatACT; % actual latitude
    WP_base.LonACT(6*i-5)=LonACT; % actual longitude
    WP_base.LatEXP(6*i-5)=LatEXP; % expected latitude
    WP_base.LonEXP(6*i-5)=LonEXP; % expected longitude
    WP_base.location(6*i-5) = LatACT+"_"+LonACT; % location as a string
    WP_base.elevation(6*i-5)=elevation; % elevation
    WP_base.category(6*i-5)=category; % category of onshore, offshore, or floating
    
    WP_base.year{6*i-4} = 2014;
    WP_base.turbinePower{6*i-4} = windPower_2014;
    WP_base.LatACT(6*i-4)=LatACT;
    WP_base.LonACT(6*i-4)=LonACT;
    WP_base.LatEXP(6*i-4)=LatEXP;
    WP_base.LonEXP(6*i-4)=LonEXP;
    WP_base.location(6*i-4) = LatACT+"_"+LonACT;
    WP_base.elevation(6*i-4)=elevation;
    WP_base.category(6*i-4)=category;

    WP_base.year{6*i-3} = 2015;
    WP_base.turbinePower{6*i-3} = windPower_2015;
    WP_base.LatACT(6*i-3)=LatACT;
    WP_base.LonACT(6*i-3)=LonACT;
    WP_base.LatEXP(6*i-3)=LatEXP;
    WP_base.LonEXP(6*i-3)=LonEXP;
    WP_base.location(6*i-3) = LatACT+"_"+LonACT;
    WP_base.elevation(6*i-3)=elevation;
    WP_base.category(6*i-3)=category;

    WP_base.year{6*i-2} = 2016;
    WP_base.turbinePower{6*i-2} = windPower_2016;
    WP_base.LatACT(6*i-2)=LatACT;
    WP_base.LonACT(6*i-2)=LonACT;
    WP_base.LatEXP(6*i-2)=LatEXP;
    WP_base.LonEXP(6*i-2)=LonEXP;
    WP_base.location(6*i-2) = LatACT+"_"+LonACT;
    WP_base.elevation(6*i-2)=elevation;
    WP_base.category(6*i-2)=category;

    WP_base.year{6*i-1} = 2017;
    WP_base.turbinePower{6*i-1} = windPower_2017;
    WP_base.LatACT(6*i-1)=LatACT;
    WP_base.LonACT(6*i-1)=LonACT;
    WP_base.LatEXP(6*i-1)=LatEXP;
    WP_base.LonEXP(6*i-1)=LonEXP;
    WP_base.location(6*i-1) = LatACT+"_"+LonACT;
    WP_base.elevation(6*i-1)=elevation;
    WP_base.category(6*i-1)=category;

    WP_base.year{6*i} = 2018;
    WP_base.turbinePower{6*i} = windPower_2018;
    WP_base.LatACT(6*i)=LatACT;
    WP_base.LonACT(6*i)=LonACT;
    WP_base.LatEXP(6*i)=LatEXP;
    WP_base.LonEXP(6*i)=LonEXP;
    WP_base.location(6*i) = LatACT+"_"+LonACT;
    WP_base.elevation(6*i)=elevation;
    WP_base.category(6*i)=category;
end
clearvars -except WP_base

% 2. Input cost parameters into the structure SP_base
% (the sources and explanation of the numbers are in the manuscript)
turbineCapital = 1.09 * 1.1 * 1.12E6/1000; % capital cost of onshore turbine in Dollar/kWp
turbineCapital_offshore = 1.1 * 1.87E6/1000; % capital cost of offshore turbine in Dollar/kWp
turbineCapital_floating = 1.1 * 3E6/1000; % capital cost of floating turbine in Dollar/kWp
electrolyserCapital = (400/409)*700/(1E-6); % capital cost of electrolysers in $/GW fed in
ASUCapital = 147000*(1/45)/(1E-6); % capital cost of PSA ASU in $/GW fed in
H2compCapital = (9/409)*6670/(1E-6); % capital cost of hydrogen compressors for HB in $/GW fed in
N2compCapital = 176400*(1/45)*(1/1E-6); %capital cost of nitrogen compressors in HB in $/GW fed in
HBCapital = 125000*230; % capital cost of HB plant in $ for 230 t/day plant
H2storageCapital = 900; % captial cost of H2 storage in $/kg
N2storageCapital = 32; % capital cost of N2 storage in $/kgg
BattstorageCapital = 150*1.1; % capital cost of battery energy capacity in $/kWh
BattpowerCapital = 270*1.1/(1E-6); % capital cost of battery power capacity in $/GW fed in
turbine_OnM = 1.1 * 14000; % fixed O&M of onshore turbine in Dollar per MW per year
turbine_OnM_offshore = 1.1 * 50000; % fixed O&M of offshore turbine
turbine_OnM_floating = 0.3; % fixed O&M of floating turbine in percent of ANNUALIZED capital cost such that O&M is 20% of total cost
elect_OnM = 0.02; % fixed O&M for electrolyser in percent of capital cost per year
HB_OnM=0.02; % fixed O&M for HB plant in percent of capital cost per year. The parameter is also used for O&M of compressors and ASU. 
discountRate = 0.07; % discount rate of capital
opYear = 25; % lifetime of the plant

WP_base.turbineCapital=turbineCapital;
WP_base.turbineCapital_offshore=turbineCapital_offshore;
WP_base.turbineCapital_floating=turbineCapital_floating;
WP_base.electrolyserCapital=electrolyserCapital;
WP_base.ASUCapital=ASUCapital;
WP_base.H2compCapital=H2compCapital;
WP_base.N2compCapital=N2compCapital;
WP_base.HBCapital=HBCapital;
WP_base.H2storageCapital=H2storageCapital;
WP_base.N2storageCapital=N2storageCapital;
WP_base.BattstorageCapital=BattstorageCapital;
WP_base.BattpowerCapital=BattpowerCapital;
WP_base.turbine_OnM=turbine_OnM;
WP_base.turbine_OnM_offshore=turbine_OnM_offshore;
WP_base.turbine_OnM_floating=turbine_OnM_floating;
WP_base.elect_OnM=elect_OnM;
WP_base.HB_OnM=HB_OnM;
WP_base.discountRate=discountRate;
WP_base.opYear=opYear;
% convert the O&M of the turbines in percent of total cost to make it in
% line with the case of floating turbine
WP_base.turbine_OnM_perc = WP_base.turbine_OnM/(1000*annual*WP_base.turbineCapital);
WP_base.turbine_OnM_offshore_perc = WP_base.turbine_OnM_offshore/(1000*annual*WP_base.turbineCapital_offshore);
WP_base.turbine_OnM_floating_perc = WP_base.turbine_OnM_floating;

% 3. Run an algorithm to solve for the turbine size, H2 production size, 
% battery power size, H2 buffer size, and battery buffer size for the case of 
% no curtailment. (i.e. close the energy balance). Store these results in stucture WP_nocurt. 
% Plant size (constant energy requirement) is 100MW. 

WP_nocurt=struct; % intialize strucutre for no curtailment
WP_nocurt.filename=WP_base.filename; % transfer some variables from WP_base
WP_nocurt.year=WP_base.year;
WP_nocurt.location=WP_base.location;
WP_nocurt.category=WP_base.category;

PReq = 100; % MW of constant supply required for plant
parfor i = 1:length(WP_base.turbinePower) % cycle through each of the power profiles in WP_base
    turbinePower = WP_base.turbinePower{i};
    N2_size = (22.5/640); %set the size of the N2 generation to be stochiometric to the ammonia production (relative size, just as for H2 size and Batt size)
    %find the require turbine size to generate the minimum amount of energy
    turbinekWpGuess = 1e6; % kWp
    x0 = turbinekWpGuess;
    inventorySum = @(x) sum((turbinePower*x/1e9 - PReq/1000)*3600); % function defining energy balance in GJ
    x = fzero(@(x) inventorySum (x), x0); % solve for panel kWp to close energy balance
    turbinekWp(i) = x; % store panel kWp in SP_nocurt
    
    E_supply = (x*turbinePower)/1E9;%energy supply in GW
    % Using the required turbine size, solve for the relative size of the H2
    % production and battery power, allowing for the relative power
    % commitment to each to vary over time
    windProblem=optimproblem('ObjectiveSense','minimize'); % initialize optimization problem
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
WP_nocurt.turbinekWp = turbinekWp;
WP_nocurt.H2SizeGW = H2SizeGW;
WP_nocurt.N2SizeGW = N2SizeGW;
WP_nocurt.BattSizeGW = BattSizeGW;
WP_nocurt.H2bufferkg = H2bufferkg;
WP_nocurt.N2bufferkg = N2bufferkg;
WP_nocurt.batterykWh = batterykWh;

% 4. Calculate the cost of each process component for the case of 
% No Curtailed energy, and store in WP_nocurt
len=length(WP_nocurt.filename);
annual=discountRate/(1-((1+discountRate)^(-opYear))); % annualization factor for capital costs
% create array of turbine capital and O&M for each entry in WP_nocurt
% according to the category of the entry
turbineCost=zeros(len,1);
turbineOnM=zeros(len,1);
for k = 1:length(WP_nocurt.turbinekWp)
    if WP_nocurt.category(k)=="onshore"
        turbineCost(k) = annual * WP_nocurt.turbinekWp(k) * turbineCapital;
        turbineOnM(k) = WP_nocurt.turbinekWp(k) * (1/1000) * turbine_OnM;
    elseif WP_nocurt.category(k)=="offshore" 
        turbineCost(k) = annual * WP_nocurt.turbinekWp(k) * turbineCapital_offshore;
        turbineOnM(k) = WP_nocurt.turbinekWp(k) * (1/1000) * turbine_OnM_offshore;
    elseif WP_nocurt.category(k)=="floating"
        turbineCost(k) = annual * WP_nocurt.turbinekWp(k) * turbineCapital_floating;
        turbineOnM(k) = annual * WP_nocurt.turbinekWp(k) * turbineCapital_floating * turbine_OnM_floating;
    end
end
electrolyserCost = annual * WP_nocurt.H2SizeGW' * electrolyserCapital; % electrolyser capital
H2compCost = annual * WP_nocurt.H2SizeGW' * H2compCapital; % H2 compressor capital
ASUCost = annual * WP_nocurt.N2SizeGW' * ASUCapital; % PSA ASU capital
N2compCost = annual * WP_nocurt.N2SizeGW' * N2compCapital; % N2 compressor capital
H2storeCost = annual * WP_nocurt.H2bufferkg' * H2storageCapital; % H2 storage capital
N2storeCost = annual * WP_nocurt.N2bufferkg' * N2storageCapital; % N2 storage capital
BattstoreCost = annual * (WP_nocurt.batterykWh' * BattstorageCapital + WP_nocurt.BattSizeGW' * BattpowerCapital); % battery capital
HBCost = ones(len,1) * annual * HBCapital; % HB capital
electOnM = elect_OnM * (WP_nocurt.H2SizeGW' * electrolyserCapital + WP_nocurt.H2SizeGW' * H2compCapital); % H2 production O&M
ASUOnM = HB_OnM * (WP_nocurt.N2SizeGW' * ASUCapital + WP_nocurt.N2SizeGW' * N2compCapital); % N2 production O&M
HBOnM = ones(len,1)* HB_OnM * HBCapital; % HB O&M
WP_nocurt.totalCapital = (turbineCost + electrolyserCost + H2compCost + ASUCost + N2compCost + N2storeCost + H2storeCost + BattstoreCost + HBCost); % total annualized capital in Dollars
WP_nocurt.LCOA = (WP_nocurt.totalCapital + electOnM + turbineOnM + HBOnM+ASUOnM)/(230*365); % LCOA in Dollars/tonne
locs=WP_nocurt.location';
WP_nocurt.costBreakdown=table(locs,turbineCost,electrolyserCost,H2compCost,ASUCost,N2compCost,H2storeCost,N2storeCost,BattstoreCost,HBCost,turbineOnM,electOnM,ASUOnM,HBOnM); % cost parameters stored in table
LCOA_byLoc=[];% initialized array of LCOA sorted by location
unique_loc=unique(WP_nocurt.location,'stable'); % create array of unique locations
for i=1:length(unique_loc) % cycle through every unique location and original location array
    for ii=1:length(WP_nocurt.location)
        if unique_loc(i)==WP_nocurt.location(ii) % check if the locations are the same
            % add LCOA to array LCOA_byLoc by column according to the year
            if WP_nocurt.year{ii}==2013
                LCOA_byLoc(i,1)=WP_nocurt.LCOA(ii);
            elseif WP_nocurt.year{ii}==2014
                LCOA_byLoc(i,2)=WP_nocurt.LCOA(ii);
            elseif WP_nocurt.year{ii}==2015
                LCOA_byLoc(i,3)=WP_nocurt.LCOA(ii);
            elseif WP_nocurt.year{ii}==2016
                LCOA_byLoc(i,4)=WP_nocurt.LCOA(ii);
            elseif WP_nocurt.year{ii}==2017
                LCOA_byLoc(i,5)=WP_nocurt.LCOA(ii);
            elseif WP_nocurt.year{ii}==2018
                LCOA_byLoc(i,6)=WP_nocurt.LCOA(ii);
            end
            unique_cat(i)=WP_nocurt.category(ii);
        end
    end
end
k=1;
% create additional array of only the locations and LCOA that are not
% floating
for i=1:length(unique_loc)
    if unique_cat(i) ~="floating"
        unique_loc_noFloat(k,:)=unique_loc(i);
        LCOA_byLoc_noFloat(k,:)=LCOA_byLoc(i,:);
        k=k+1;
    end
end
WP_nocurt.LCOA_byLoc=LCOA_byLoc;
WP_nocurt.unique_loc=unique_loc';
WP_nocurt.unique_cat=unique_cat';
WP_nocurt.unique_loc_noFloat=unique_loc_noFloat;
WP_nocurt.LCOA_byLoc_noFloat=LCOA_byLoc_noFloat;

save('WP_base','WP_base') % save WP_base structure
save('WP_nocurt','WP_nocurt'); % save WP_nocurt structure
clearvars