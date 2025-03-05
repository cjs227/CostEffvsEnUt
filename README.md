#Guide to code used in the manuscript "Cost efficiency vs energy utilization in green ammonia production from intermittent renewables"

Written and executed in MATLAB 2023.

##Raw solar data processed using **Solar_Data_Processing.m**
This script requires solar profile CSV data downloaded from the PVGIS platform (https://re.jrc.ec.europa.eu/pvg_tools/en/) to be in the folder **Raw_Data_Solar**. It also requires the MATLAB mat file **gebco_variables**,  
This script executes 5 sections and save **SP_base** and **SP_nocurt**. 
1. Extract all the data from the CSV files and sort in the stucture SP_base
2. Remove outliers from the power profiles in SP_base structure
3. Input cost parameters in the structure SP_base
4. Solve for the panel, H2 production, N2 production, H2 storage, and battery size for the case of No curtailment. Store in the structure SP_nocurt.
5. Calculate the cost of process components for each case and LCOA. Sort LCOA according to location. 

##Raw wind data processed  **Wind_Data_Processing.m**. 
This script requires wind profile NC data downloaded from the New European Wind Atlas platform (https://map.neweuropeanwindatlas.eu/) to be in the folder **Raw_Data_Wind**. Is also requires the MATLAB mat file **gebco_varibales**. 
This script executes 5 sections and save **WP_base** and **WP_nocurt**. 
1. Extract all the data from the NC files and sort in the stucture WP_base 
2. Input cost parameters in the structure WP_base
3. Solve for the turbine, H2 production, N2 production, H2 storage, and battery size for the case of No curtailment. Store in the stucture WP_nocurt.
4. Calculate the cost of process components for each case and LCOA. Sort LCOA according to location. 

##Capacity planning with curtailment of solar energy in **Curtailed_Solar_Analysis.m**
It requires the mat files **SP_base** and **SP_nocurt** in the path.
This script executes 3 sections and saves SP_curt.
1. Load the necessary mat variables and extract cost parameters. 
2. Execute the linear solver for optimized process design for each entry in SP_base using a parallelized loop.
3. Calculate the cost components for each optimization solution and store in SP_curt. 

##Capacity planning with curtailment of wind energy in **Curtailed_Wind_Analysis.m**
It requires the mat files **WP_base** and **WP_nocurt** in the path.
This script executes 3 sections and saves **WP_curt**.
1. Load the necessary mat variables and extracting cost parameters. 
2. Execute the linear solver for optimized process design for each entry in WP_base using a parallelized loop.
3. Calculate the cost components for each optimization solution and store in WP_curt.

##Capacity planning with curtailment and combined solar/wind in **Curtailed_Combined_Analysis.m**
It requires the mat files **SP_base** and **WP_base** in the path.
This script executes 7 sections and saves **CP_base_solar**, **CP_base_wind**, and **CP_curt**.
1. Load data for solar and wind and save in combined structures CP_base_wind and CP_base_solar
2. Initialize structure for "curtailment", CP_curt. Configure and solve linear optimization problem for optimal process design with combined solar and wind with curtailment using a parallelized loop. 
3. Find the cost components for the results of each optimization solution. 

##Capacity planning with solar energy, curtailment, and ramping of the HB process in **Ramping_Solar_Analysis.m**
It requires the mat files **SP_base**, **SP_nocurt** and **SP_curt** in the path. 
This script has 3 sections and saves **SP_ramp** (structure), **NH3Prodramp_solar** (cell array of HB production over time in optimal solution) and **Extramp_solar** (cell array of energy extraction over time in optimal solution) in the file **SP_ramp_results.mat**
1. Load necessary data and initialize structure SP_ramp.
2. Configure and solve linear optimization problem for optimal process design with curtailment while ramping the HB process between 60-100% of design capacity at a maximum rate of 5% design capaicty per hour using a parallelized loop. 
3. Find the cost components for the results of each optimization solution. 

##Capacity planning with wind energy, curtailment, and ramping of the HB process in **Ramping_Wind_Analysis.m**
It requires the mat files **WP_base**, **WP_nocurt** and **WP_curt** in the path.
This script has 3 sections and saves **WP_ramp** (structure), **NH3Prodramp_wind** (cell array of HB production over time in optimal solution) and **Extramp_wind** (cell array of energy extraction over time in optimal solution) in he file **WP_ramp_results.mat**
1. Load necessary date and intialize strucutre WP_ramp.
2. Configure and solve linear optimization problem for optimal process design with curtailment while ramping the HB process between 60-100% of design capacity at a maximum rate of 5% design capaicty per hour using a parallelized loop.
3. Find the cost components for the results of each optimization solution. 

##Cost and value analysis of solar energy in **Utilization_Costs_Solar.m**
It requires the mat files **SP_base**, **SP_nocurt**, **SP_curt** and **SP_costs** in the path. The **SP_costs** structure must already contain vectors which contain the **lat_lon** codes for the location in the top10, middle10, and bottom10. (e.g. SP_costs.top10_loc)
This script has 8 sections and saves **SP_costs** (structure). 
1. Load the necessary data.
2. 
a. Solve for the optimized process design for each year at each location categorized as middle10, for curtailment level ranging from 0-99% in steps of 1% curtailment.
b. Same as 2a for locations categorized as top10. 
c. Same as 2a for locations categoried as bottom10. 
3. 
a. Using the optimizations in 2a, Calculate the cost of utilizing energy for each process component as well as the total cost of utilization and total energy cost. 
b. Using the optimizations in 2b, Calculate the cost of utilizing energy for each process component as well as the total cost of utilization and total energy cost. 
c. Using the optimizations in 2c, Calculate the cost of utilizing energy for each process component as well as the total cost of utilization and total energy cost.
4. Using the optimization results and utilization costs in section 3 to find the monthly costs and values for two select locations.  

##Cost and value analysis of wind energy in **Utilization_Costs_Wind.m**
It requires the mat files **WP_base**, **WP_nocurt**, **WP_curt** and **WP_costs** in the path. The **WP_costs** structure must already contain vectors which contain the **lat_lon** codes for the location in the top10, middle10, and bottom10. (e.g. WP_costs.top10_loc)
This script has 8 sections and saves **WP_costs** (structure). 
1. Load the necessary data.
2. 
a. Solve for the optimized process design for each year at each location categorized as middle10, for curtailment level ranging from 0-99% in steps of 1% curtailment.
b. Same as 2a for locations categorized as top10. 
c. Same as 2a for locations categoried as bottom10. 
3. 
a. Using the optimizations in 2a, Calculate the cost of utilizing energy for each process component as well as the total cost of utilization and total energy cost. 
b. Using the optimizations in 2b, Calculate the cost of utilizing energy for each process component as well as the total cost of utilization and total energy cost. 
c. Using the optimizations in 2c, Calculate the cost of utilizing energy for each process component as well as the total cost of utilization and total energy cost.
4. Using the optimization results and utilization costs in section 3 to find the monthly costs and values for two select locations.

##Cost and value analysis of solar energy with ramping to 40% min. capacity in **CostswRamp_Solar.m**
It requires the mat files **SP_base**, **SP_nocurt**, **SP_curt** and **SP_costs** in the path. 
This script has 9 sections and saves **SP_rampcosts** (structure). 
1. Load the necssary data. 
2. 
a. Solve for the optimized process design with ramping to 40% mininimum capacity for each year at each location categorized as top10, for curtailment level ranging from 0-99% in steps of 1% curtailment. 
b. Same as 2a for locations categorized as middle10. 
c. Same as 2a for locations categorized as bottom10. 
3. 
a. Using the optimizations in 2a, calculate the cost of utilizing energy for each process component as well as the total cost of utilization and total energy cost.
b. Using the optimizations in 2b, calculate the cost of utilizing energy for each process component as well as the total cost of utilization and total energy cost.
c. Using the optimizations in 2c, calculate the cost of utilizing energy for each process component as well as the total cost of utilization and total energy cost.
4. 
a. Using the optimization results and utilization costs in 2a to find the monthly costs and values for a select top10 location.
b. Using the optimization results and utilization costs in 2b to find the monthly costs and values for a select middle10 location.

##Cost and value analysis of wind energy with ramping to 40% min. capacity in **CostswRamp_Wind.m**
It requires the mat files **WP_base**, **WP_nocurt**, **WP_curt** and **WP_costs** in the path. 
This script has 9 sections and saves **WP_rampcosts** (structure). 
1. Load the necessary data. 
2. 
a. Solve for the optimized process design with ramping to 40% mininimum capacity for each year at each location categorized as top10, for curtailment level ranging from 0-99% in steps of 1% curtailment. 
b. Same as 2a for locations categorized as middle10. 
c. Same as 2a for locations categorized as bottom10. 
3. 
a. Using the optimizations in 2a, calculate the cost of utilizing energy for each process component as well as the total cost of utilization and total energy cost.
b. Using the optimizations in 2b, calculate the cost of utilizing energy for each process component as well as the total cost of utilization and total energy cost.
c. Using the optimizations in 2c, calculate the cost of utilizing energy for each process component as well as the total cost of utilization and total energy cost.
4. 
a. Using the optimization results and utilization costs in 2a to find the monthly costs and values for a select top10 location.
b. Using the optimization results and utilization costs in 2b to find the monthly costs and values for a select middle10 location.

##Cost and value analysis with combined solar/wind energy in **Utilization_Costs_Combined.m**
It requires **CP_base_solar.mat**, **CP_base_wind.mat** and **CP_costs**. The CP_costs structure must already contain vectors which contain the **lat_lon** codes for the locations categorized as top, topsolar, topwind and average. These are defined as the lowest 15% of LCOA for solar and wind individually, lowest 10% of solar LCOA and 45-55% percentile of wind LCOA, lowest 10% of wind LCOA and 45-55% percentile of solar LCOA, and 45-55% percentile of solar and wind LCOA individually, respectively. 
This script has 29 sections and saves the structure **CP_costs**. 
1. Load the necessary data. 
2. 
a-1. Solve for the optimized process design with combined solar and wind at average locations. 
a-2. Same locations and procedure as 2a-1, but only using solar energy at the location. 
a-3. Same locations and procedure as 2a-1, but only using wind energy at the location.
b-1 - b-3. Same as a-1 - a-3, but for locations categorized as top solar. 
c-1 - c-3. Same as a-1 - a-3, but for locations categorized as top wind. 
d-1 - d-3. Same as a-1 - a-3, but for locations categorized as top.
3. 
a-1. Calculate utilization costs as a function of percent curtailment for top locations using previous optimization results of combined solar and wind.
a-2. Calculate utilization costs as a function of percent curtailment for top locations using previous optimization results of only solar. 
a-3. Calculate utilization costs as a function of percent curtailment for top locations using previous optimization results of only wind. 
b-1 - b-3. Same as a-1 - a-3, except for locations categorized as top solar. 
c-1 - c-3. Same as a-1 - a-3, except for locations categorized as top wind. 
d-1 - d-3. Same as a-1 - a-3, except for location categorized as average. 
4. 
a. Calculate the amount of solar and wind energy utilized monthly for combined solar and wind at average locations. 
b. Calculate the amount of solar and wind energy utilized monthly for combined solar and wind at top solar locations. 
c. Calculate the amount of solar and wind energy utilized monthly for combined solar and wind at top wind locations. 
d. Calculate the amount of solar and wind energy utilized monthly for combined solar and wind at top locations. 

##Generate maps of LCOA and fraction curtailment in **Map_Generation.mlx**.
It requires the mat files **SP_nocurt**, **WP_nocurt**, **SP_curt**, **WP_curt**, **SP_ramp**, **WP_ramp**, **CP_curt**, and **PinkHeatMap**, and **YellowBluemap**.
This live script has 9 sections.
1. Generate map of LCOA with solar and no curtailment. Saves map in images "Solar NoCurt Map.png".
2. Generate maps of LCOA and fraction curtailment with solar and curtailment.  Saves maps in "Solar Curt Map.png" and "Solar Fraction Curtailment Map.png". 
3. Generate maps of LCOA, fraction curtailment, and change in fraction curtailment with solar energy and ramping of the HB process. Saves maps in "Solar Ramp Map.png" "Solar Ramp Curt Map.png" and "Solar RampCurt Change Curt Map.png".
4. Generate map of LCOA with wind and no curtailment. Saves map in images "Wind NoCurt Map.png".
5. Generate maps of LCOA and fraction curtailment with wind and curtailment. Saves maps in "Wind Curt Map.png" and "Wind Fraction Curtailment Map.png". 
6. Generate maps of LCOA, fraction curtailment, and change in fraction curtailment with wind energy and ramping of the HB process. Saves maps in "Wind Ramp Map.png" "Wind Ramp Curt Map.png" and "Wind RampCurt Change Curt Map.png". 
7. Generate maps of LCOA, fraction curtailment, and fraction of energy from solar energy with combined solar and wind energy and curtailment. Saves maps in "Combined Curt Map.png" "Combined Curtailment Map.png" and "Combined Curt FractSolar Map.png",
