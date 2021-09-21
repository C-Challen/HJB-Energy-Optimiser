function HJB_main()      
%% Setup

    global T Nt Nc Nv Elec_Forecast Elec_Price Elec_Sell Therm_Forecast Therm_Price PV_Forecast Cap_E Cap_E_Min Cap_T Cap_T_Min;    

    % USER DEFINED INPUTS
    Cap_E = 4.8;                % Capacity of Electrical storage in kWh
    Cap_E_Min = Cap_E * 0.2;    % Minimum electrical capacity limit in kWh
    Cap_T = 3.5;                % Capacity of Thermal storage in kWh
    Cap_T_Min = 0;              % Minimum thermal capacity limit in kWh
    dt = 1;                     % Timestep size in minutes(T * 60)/(Nt-1); 
    T = 24;                     % Duration in hours
    Nc = 100;                   % Number of discrete SOC steps
    Nv = 100;                   % Number of discrete control steps        
    
    
    Nt = (T * 60) / dt + 1;      % Number of timesteps  
    xlswrite('Crest_Demand_Model_v2.3.1.xlsm',[T,dt],'Results','Q2:R2');
    Data = xlsread('CREST_Demand_Model_v2.3.1.xlsm','Results',['E2:L' num2str(Nt)+1]);
    Data(Nt,:) = Data(Nt-1,:);

    ES_SOC = zeros(1,Nt+1); 
    TS_SOC = zeros(1,Nt+1); 
    ES_Plan = zeros(1,Nt); 
    TS_Plan = zeros(1,Nt);  
    Cost_Plan = zeros(1,Nt);
      
%% PV generation

    PV_Forecast = Data(:,2)';
    
%% Electrical forecast 

    Elec_Forecast = Data(:,1)';
    Elec_Price = Data(:,7)';
    Elec_Sell = Data(:,8)';

%% Thermal forecast    
    
    Therm_Forecast = Data(:,6)';   
    Therm_Price = Elec_Price; % The thermal storage is the only water heater, there is no gas supply  

%% Starting state of charge

    ES_SOC(1,1) = Cap_E_Min * 0.5; % Starting at 50% SOC 
    TS_SOC(1,1) = Cap_T_Min * 0.5; 
    
%% Solver   

    close all;
    HJB();
    vpp = HJB();
    vpp.solve();
    
%% Planner

    ES_Opt = vpp.ve;
    TS_Opt = vpp.vh;
    Cost = vpp.u;

    for t = 1:Nt
       
        ES_SOC_Index = min(max(round((ES_SOC(1,t)-Cap_E_Min)/(Cap_E-Cap_E_Min)*Nc),1),Nc); % Convert HJB SOC values to percentages
        TS_SOC_Index = min(max(round(TS_SOC(1,t)/Cap_T*Nc),1),Nc);
        ES_Plan(1,t) = ES_Opt(TS_SOC_Index,ES_SOC_Index,t); % Read optimal operation by using SOC and t values as co-ordinates
        TS_Plan(1,t) = TS_Opt(TS_SOC_Index,ES_SOC_Index,t);
        Cost_Plan(1,t) = Cost(TS_SOC_Index,ES_SOC_Index,t);
        
        if ES_Plan(1,t) > 0 % If ES is discharging
            ES_SOC(1,t+1) = max(min(ES_SOC(1,t) - ES_Plan(1,t) / vpp.Eed * dt/60 - vpp.De * dt,Cap_E),Cap_E_Min); % SOC after discharge and dissipation
        else
            ES_SOC(1,t+1) = max(min(ES_SOC(1,t) - ES_Plan(1,t) * vpp.Eec * dt/60 - vpp.De * dt,Cap_E),Cap_E_Min); % SOC after charge and dissipation
        end
        
        if Therm_Forecast(1,t) < 0 % If TS is discharging
            TS_SOC(1,t+1) = max(min(TS_SOC(1,t) + Therm_Forecast(t) / vpp.Ehd * dt/60 - vpp.Dh * dt,Cap_T),Cap_T_Min); % SOC after discharge and dissipation
        else
            TS_SOC(1,t+1) = max(min(TS_SOC(1,t) - TS_Plan(1,t) * vpp.Ehd * dt/60 - vpp.Dh * dt,Cap_T),Cap_T_Min); % SOC after charge and dissipation            
        end
    end

%% Results
    
    clearvars ES_SOC_Temp t TS_SOC_Temp Data ES_Opt TS_Opt ES_SOC_Index TS_SOC_Index ans
    save('HJB_Output');
  
end