classdef HJB < handle
    % The class HJB sets up a numerical approximation of the HJB equation
    % associated to the optimal control of a virtual power plant with
    % one electrical and one hear storage device.
    %
    % Max Jensen, Chris Challen, July 2020
    % 
    % We solve the HJB equation numerically with a semi-Lagrangian
    % method following (8.48) of "Semi-Lagrangian Approximation Schemes for
    % Linear and Hamilton-Jacobi Equations" by Falcone, Ferretti, see
    % http://epubs.siam.org/doi/book/10.1137/1.9781611973051
    properties
        % Properties of the model
        Fec         % forecast of the electricity consumption
        Fhc         % forecast of the heat consumption
        Fpv         % forecast of the photovoltaic generation
        Fepb        % forecast of the electricity buying price
        Feps        % forecast of the electricity selling price
        CeM         % capacity of the electrical storage
        Cem         % minimum charge level of electrical storage
        ChM         % capacity of the heat storage
        Chm         % minimum charge level of heat storage
        Pesm        % lower electrical storage power flow constraint
        PesM        % upper electrical storage power flow constraint
        Phsm        % lower heat storage discharge power flow constraint
        PhsM        % upper heat storage discharge power flow constraint
        Ewh         % conversion efficiency for water heater
        Eec         % conversion efficiency for electric charging
        Eed         % conversion efficiency for electric discharging
        Ehc         % conversion efficiency for heat charging
        Ehd         % conversion efficiency for heat discharging
        De          % electric charge dissipation
        Dh          % heat charge dissipation
        T           % forecast period, e.g. 24 hours
        u           % value function
        ve          % optimal control ve
        vh          % optimal control vh
        % Properties of the numerical discretisation
        Nt          % number of numerical time steps
        Nc          % number of numerical charge steps
        Nv          % number of numerical control steps
        tt          % time steps of computation
        Ce          % nodes discretising electric charge
        Ch          % nodes discretising heat charge
        CeMG        % meshgrid of Ce
        ChMG        % meshgrid of Ch
        vee         % nodes discretising ve control parameter
        vhh         % nodes discretising vh control parameter
        vcc         % nodes discretising vc control parameter
        ine         % precomputed interpolation nodes: Ce coordinate
        inh         % precomputed interpolation nodes: Ch coordinate
        ineMG       % meshgrid of ine
        inhMG       % meshgrid of inh
        scost       % precomputed stationary running costs
        dt          % time step size
    end
    methods
        function obj = HJB()
            % object constructor
            % the first part sets the values for the model 
            global T Nt Nc Nv Elec_Forecast Elec_Price Elec_Sell Therm_Forecast PV_Forecast Cap_E Cap_E_Min Cap_T Cap_T_Min; 
            obj.T = T;
            obj.Nt = Nt; 
            obj.Nc    = Nc; 
            obj.Nv    = Nv; 
            obj.CeM   = Cap_E;          % kWh 
            obj.Cem   = Cap_E_Min;      % kWh
            obj.ChM   = Cap_T;          % kWh
            obj.Chm   = Cap_T_Min;      % kWh
            obj.Pesm  = - 0.750;        % kW - charging 
            obj.PesM  = 0.850;          % kW - discharging 
            obj.Phsm  = -2.8;           % kW - charging (electrical energy)
            obj.PhsM  = 35.0;           % kW - discharging (thermal energy)
            obj.Eec   = 0.865; 
            obj.Eed   = 0.88; 
            obj.Ewh   = 0.95;            
            obj.Ehc   = 0.961;                   
            obj.Ehd   = 0.961;                    
            obj.De    = 0.00016 / obj.CeM;                   
            obj.Dh    = 0.00036 / obj.ChM;                  
            % precompute quantities which can be reused across multiple forecasts
            obj.precompute();            
            % download forecasts for the model
            % remember the convention that energy produced at a point in a
            % network has a positive sign at that point. Therefore, we expect
            % the consumption forecasts usually to carry a negative sign
            obj.Fec   = Elec_Forecast(1:obj.Nt); 
            obj.Fhc   = Therm_Forecast(1:obj.Nt);
            obj.Fpv   = PV_Forecast(1:obj.Nt);
            obj.Fepb  = Elec_Price(1:obj.Nt);
            obj.Feps  = Elec_Sell(1:obj.Nt); 
        end
        function dCe = Re(obj, Pes, ce)
            % rate of change of electric charge given an applied power
            if Pes > 0
                dCe = - Pes / obj.Eed - obj.De * ce; % Discharging from electrical storage
            else
                dCe = - obj.Eec * Pes - obj.De * ce; % Charging to electrical storage
            end
        end
        function dCh = Rh(obj, Phs, ch)
            % rate of change of heat charge given an applied power
            if Phs > 0
                dCh = - Phs / obj.Ehd - obj.Dh * ch; % Discharing from thermal storage (thermal energy)
            else
                dCh = - obj.Ehc * Phs - obj.Dh * ch; % Charging to thermal storage (electrical energy)
            end
        end
        function obj = precompute(obj)
            % now precompute quantities which can be reused across multiple forecasts
            % initialise data structures
            obj.dt    = obj.T / (obj.Nt-1);
            obj.tt    = linspace(0,obj.T,obj.Nt);
            obj.Ce    = linspace(obj.Cem,obj.CeM,obj.Nc);
            obj.Ch    = linspace(obj.Chm,obj.ChM,obj.Nc);
            [obj.CeMG, obj.ChMG] = meshgrid(obj.Ce,obj.Ch); 
            obj.vee   = [linspace(obj.Pesm,0,ceil(obj.Nv * 1/2)),linspace(obj.PesM/(obj.Nv*1/2),obj.PesM,floor(obj.Nv * 1/2))];
            obj.vhh   = [linspace(obj.Phsm,-1,obj.Nv/5),linspace(-1+1/(obj.Nv/5),0,obj.Nv/5),linspace(obj.PhsM/3.5/obj.Nv*2.5,obj.PhsM/3.5,obj.Nv/2.5),linspace(obj.PhsM/7*3,obj.PhsM,obj.Nv/5)];
            obj.u     = zeros(obj.Nc,obj.Nc,obj.Nt);
            obj.ve    = zeros(obj.Nc,obj.Nc,obj.Nt);
            obj.vh    = zeros(obj.Nc,obj.Nc,obj.Nt);
            obj.ine   = zeros(obj.Nc,obj.Nv);
            obj.inh   = zeros(obj.Nc,obj.Nv);
            obj.scost = zeros(obj.Nc,obj.Nc,obj.Nv,obj.Nv);
            obj.ineMG = zeros(obj.Nc,obj.Nc,obj.Nv); % initialise mesh grids
            obj.inhMG = zeros(obj.Nc,obj.Nc,obj.Nv);
            
            for ive = 1:obj.Nv % control process ve
                for ice = 1:obj.Nc % electric charge level
                    obj.ine(ice,ive) = obj.Ce(ice) + obj.dt * obj.Re(obj.vee(ive), obj.Ce(ice));
                end
                obj.ineMG(:,:,ive) = ones(obj.Nc,1) * obj.ine(:,ive)' ; 
            end
            
            for ivh = 1:obj.Nv % control process vh
                for ich = 1:obj.Nc % heat charge level
                    obj.inh(ich,ivh) = obj.Ch(ich) + obj.dt * obj.Rh(obj.vhh(ivh), obj.Ch(ich));
                end
                obj.inhMG(:,:,ivh) = obj.inh(:,ivh) * ones(1,obj.Nc); 
            end
        end
        function cost = tcost(obj, t, ve, vh)
            % this is the part of the running cost whose terms depend on
            % time and cannot be precomputed in obj.scost
            
            % forecasted energy cost
            cost = - (obj.Fec(t) + obj.Fpv(t) + ve) - min((obj.Fhc(t) + vh) / obj.Ewh,0);

            if cost > 0
                cost = cost * obj.Fepb(t);
            else
                cost = cost * obj.Feps(t);
            end

            % battery degradation represented through a generic quadratic penalty
            cost = cost + 0.000001 * ve^2;               
           
            % selection term: all other costs being the same, choose smaller powers
            cost = cost + 0.000001 * (ve^2 + vh^2);            
        end
        function obj = solve(obj)
            % solve the HJB equation numerically with a semi-Lagrangian
            % method: (8.48) in "Semi-Lagrangian Approximation Schemes for
            % Linear and Hamilton-Jacobi Equations" by Falcone, Ferretti.
            
            ep = - max(obj.Fepb); % electricity price at final time
            
            % Now assign final time conditions (value of stored energy stored at end of T
            % divided by efficiency loss in the future discharing of this energy)
            obj.u(:,:,obj.Nt) = ep * linspace(0,obj.CeM-obj.Cem,obj.Nc)' * ones(1,obj.Nc) * obj.Eed ...
                + ep * ones(obj.Nc,1) * linspace(obj.Chm,obj.ChM,obj.Nc) * obj.Ehd;
            
            % A loop to go over all time steps backwards
            fprintf('time step (of %d):    ',obj.Nt);
            tic;
            for t = obj.Nt-1:-1:1             
                fprintf('\b\b\b\b\b %4d',t);
                % candidate solution, starting with artificially high values
                ucand = 1000 * ones(obj.Nc,obj.Nc);
                % there arrays will hold the optimal control
                vecand = zeros(obj.Nc,obj.Nc);
                vhcand = zeros(obj.Nc,obj.Nc);
                for ivh = 1:obj.Nv % control process ve
                    for ive = 1:obj.Nv % control process vh
                            % This implements the term inside on the Min in                            
                            % (8.48) for a given control
                            
                            unew = obj.dt * obj.tcost(t, obj.vee(ive), obj.vhh(ivh)) + ... % New cost
                                interp2(obj.CeMG,obj.ChMG,obj.u(:,:,t+1),obj.ineMG(:,:,ive),obj.inhMG(:,:,ivh)); % Plus old cost 
                            
                            % Gives total cost to reach the final time
                            % check were the new control achieves lower
                            % values than in the computations before
                            m = ucand > unew;
                            % updave the candidate solution and the associated controls
                            ucand(m) = unew(m);
                            vecand(m) = ive;
                            vhcand(m) = ivh;                                                      
                    end
                end
                % save found values for time step t                
                obj.u(:,:,t) = ucand;
                obj.ve(:,:,t) = obj.vee(vecand);
                obj.vh(:,:,t) = obj.vhh(vhcand);
            end
            fprintf('\n');
            toc
        end
    end
end