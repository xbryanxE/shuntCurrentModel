classdef thermo
    properties
        %T % °C
        %P % Pa
        %x % KOH fraction
        state = state; % electrolyte state

        % Numerical derivative parameters
        dT = 2.5;

        % Ionic radius
        aR = [0.133 0.110]*1e-9; % [K+ OH-]
        z = [1 1]; % Charge number

        % Physical constants
        % Gas constant
        R = 8.314462618; % J.mol^-1.K^-1
        % Faraday constant
        Fc = 96485.332123; % A.s.mol^-1
        % Avogadro number
        Nv = 6.0221367e23; % mol^-1
        % Gravity
        g = 9.8; % m.s^-2
        % Boltzmann constant
        k = 1.380658e-23; % m^2.kg.s^-2.K^-1
        % Speed of light
        c = 299792458; % m.s^-1
        % Permitivity of free space
        eps_ref = 1/(4e-7*pi*299792458^(2)); % C^2.J^(-1).m^(-1)
        % Molecular dipole moment of water
        muw = 6.138e-30; % C.m
        % water molecular polarizability
        polar_a = 1.636e-40 % C^2.J^(-1).m^(2)
        % Elementary charge
        e = 1.602176634e-19; % C
        % Debye-Huckel parameter
        OMG = 55.51; % mol.kg^-1
        % water moleculer weight
        Mw_h2o = 18.015268; % g.mol^-1
        % oxygen molecular weight
        Mw_o2 = 32.0; % g.mol^-1
        % hydrogen molecular weight
        Mw_h2 = 2.0; % g.mol^-1
        % critical density of water
        rho_crit = 322/(18.015268*1e-3); % mol.m^(-3)
        % critical temperature of water
        T_crit = 647.096; % K
        % surface tension matrix
        Asigma = [75.4787, -0.138489, -0.336392e-3, 0.475362e-6, -0.264479e-9; ...
                  -32.8890, 1.34382, -0.910138e-2, 0.396124e-4, -0.573565e-7; ...
                  614.527, -12.8736, 0.104855, -0.449076e-3, 0.651193e-6; ...
                  -1455.06, 39.8511, -0.344234, 0.144383e-2, -0.207599e-5; ...
                  1333.62, -38.3316, 0.335129, -0.137313e-2, 0.194911e-5];
        % viscosity vector
        muVector = [0.9105535967, -0.01062211683, 4.680761561e-5, ...
                    -9.209312883e-8, 6.814919843e-11];

        dataInfo = 0;

    end % properties

    methods
        %% read parameters for correlations 
        function output = readInfoDict(obj)
            data_electrolyte = readtable("Info.xlsx", 'sheet', 'electrolyte');
            data_water = readtable("Info.xlsx", 'sheet', 'water');
            data_kinetics = readtable("Info.xlsx", 'sheet', 'kinetics');
            data_custom = readtable("Info.xlsx", 'sheet', 'custom');
            data_dimensions = readtable("Info.xlsx", 'sheet', 'dimensions');
            
            names = ["electrolyte", "water", "kinetics", "custom", "dimensions"];
            values = {data_electrolyte, data_water, data_kinetics, data_custom, data_dimensions};
            DATA = dictionary(names, values);
            obj.dataInfo = DATA;
            output = obj;
        end

        %% Parameter reader
        function p = readPar(obj, sheetname, paramList, n)
            dataDict = obj.dataInfo(sheetname);
            data = dataDict{1};
            p = zeros(n, length(paramList));
            for i=1:length(paramList)
                p(:,i) = data{1:n, paramList{i}};
            end
        end

        %% Redlich Meyer parameter B for heat capacity
        function B = RedlichMeyer(~, T, rho, m)
            % coefficients
            ccoeff = [-1.27, -0.1472, 21.2, 30.2 62., 52., -31.1, -0.0655, 0.236 ... 
                       68.2, -7.51, -1.45,  0.326, 0.628, -0.0534, 0.99];
            % functions
            F = {@(x, y, z) y.^(2)./(0.647 - x);
                 @(x, y, z) y.^(3)./x.^(3);
                 @(x, y, z) y.^(4);
                 @(x, y, z) y.^(3)./log(x);
                 @(x, y, z) y.^(4).*x.^(2);
                 @(x, y, z) x.^(2);
                 @(x, y, z) y.^(2).*x.^(2).*z;
                 @(x, y, z) y.^(2)./(0.647 - x).*z;
                 @(x, y, z) y./x.*z;
                 @(x, y, z) y.^(3).*x.^(3).*z;
                 @(x, y, z) x.^(2).*z.^(1.5);
                 @(x, y, z) y.^(3).*x.*z.^(1.5);
                 @(x, y, z) y.^(3).*log(x).*z.^(1.5);
                 @(x, y, z) y.^(3)./(0.647 - x).*z.^(1.5);
                 @(x, y, z) y.^(2)./(0.647 - x).*z.^(2);
                 @(x, y, z) x.^(2).*z.^(2);};
            
            B = 0;
            for ii = 1:length(ccoeff)
            B = B + ccoeff(ii)*F{ii}(T*1e-3, rho*1e-3, m);
            end
        end 

        %% Osmotic coefficient from Debye Huckel Theory
        function A = funcA(obj, T, P, key)
            % Pure water density
            rho_pure = obj.state.rho;
            % dielectric constant of pure water
            eps_pure = obj.state.eps;
            if strcmp('osmotic', key) == true
                aDH = (2*pi*rho_pure*obj.Nv/1e21).^(1/2);
                bDH = (obj.e.^(2).*obj.c.^(2)./(eps_pure.*obj.k.*T)).^(3/2);
                A = 1/3*(aDH.*bDH);
            elseif strcmp('enthalpic', key) == true
                % Temperature vector for numerical integration
                T_trial = [T; T + obj.dT]; % K
                P_trial = [P; P]; % Pa
                % read pressure
                % Temperature derivative of pure water dielectric constant
                eps_trial = obj.waterProp(T_trial, P_trial, 'epsilon');
                d_eps = (eps_trial(2,:) - eps_trial(1,:))./obj.dT;
                % Temperature derivative of pure water density
                rho_trial = obj.waterProp(T_trial, P_trial, 'density');
                d_rho = (rho_trial(2,:) - rho_trial(1,:))./obj.dT;
                % Osmotic coefficient
                A_osm = obj.funcA(T, P, 'osmotic');
                % Enthalpic coefficient
                A = 6.*A_osm.*obj.R.*T.*(1 + T.*d_eps.*1/obj.state.eps ...
                    + T.*d_rho*1/(3*obj.state.rho));
            elseif strcmp('heatCapacity', key) == true
                T_trial = [T; T + obj.dT; T + 2*obj.dT]; % K
                P_trial = [P; P; P]; % Pa
                A_trial = obj.funcA(T_trial, P_trial, 'osmotic');
                A = (A_trial(1,:) - 2*A_trial(1,:) + A_trial(1,:))./(2*obj.dT);
            end
        end

        %% Function for Debye-Huckel force parameter Gf
        function gpar = funcG(~, x)
            gpar = 2.*(1 - (1 + x).*exp(-x))./(x.^2);
        end % g parameter
        
        %% Debye Huckel forces
        function DH = DebyeHuckel(obj, T, P, x)
            fxn = @(xn) x - xn*56.1056/(18*(1 - xn) + 56.1056*xn);
            x = fsolve(fxn, x);
            % Osmotic parameter
            A_osm = obj.funcA(T, P, 'osmotic');
            % Ionic strenght
            Ix = obj.ionicStrength(x);
            % read Debye Huckel coefficients
            paramList = {'rhox', 'alphax', 'Bx'};
            out = obj.readPar('electrolyte', paramList, 1);
            rhox = out(1);
            alphax = out(2);
            Bx = out(3);
            % DH contribution of solvent
            c1 = 2*A_osm.*Ix.^(3/2)/(1 + rhox.*Ix.^(1/2));
            c2 = (x.^(2)/4).*Bx.*exp(-alphax*Ix.^(1/2));
            DHw = c1 - c2;
            % DH contribution of solute
            eta1 = (2/rhox).*log(1 + rhox*Ix.^(1/2));
            eta2 = Ix.^(1/2).*(1 - 2*Ix)/(1 + rhox*Ix.^(1/2));
            c1 = -A_osm.*(eta1 + eta2);
            Zf = obj.funcG(alphax*Ix.^(1/2));
            c2 = (x/4).*Bx.*(Zf + (1 - x).*exp(-alphax*Ix.^(1/2)));
            DHx = c1 + c2;
            DH = [DHw; DHx];
        end
        
        %% Margules force contribution
        function S = Margules(obj, T, x)
            fxn = @(xn) x - xn*56.1056/(18*(1 - xn) + 56.1056*xn);
            x = fsolve(fxn, 0.1);
            paramList = {'U', 'W'};
            out = obj.readPar('electrolyte', paramList, 4);
            % U value
            Cxu = out(:,1);
            U = Cxu(1) + Cxu(2)./T + Cxu(3)*T + Cxu(4)./(647. - T);
            % W value
            Cxw = out(:,2);
            W = Cxw(1) + Cxw(2)./T + Cxw(3)*T + Cxw(4)./(647. - T);
            % Margules contribution for solvent
            Sw = x.^(2).*(W + (2*x - 1).*U);
            % Margules contribution for solute
            Sx = ((1 - x).^(2) - 1).*W + 2*(1 - x).^(2).*x.*U;
            S = [Sw; Sx];
        end

        %% Ionic strength
        function Is = ionicStrength(obj, x)
            % positive charge
            paramList = {'qPos', 'qNeg', 'nProtons', 'nElectrons'};
            out = obj.readPar('electrolyte', paramList, 1);
            zp = out(1);
            zn = out(2);
            np = out(3);
            nn = out(4);
            % fraction of dissociated species
            % positive moiety
            xp = (np/(np + nn))*x;
            % negative moiety
            xn = (nn/(np + nn))*x;
            Is = 1/2.*(xp*zp.^(2) + xn*zn.^(2));
        end

        %% Water properties
        function prop = waterProp(obj, T, P, key)
            if strcmp('epsilon', key) == true
                % vectors of coefficients for g
                Nh = [0.978224486826, ...
                      -0.957771379375, ...
                      0.237511794148, ...
                      0.714692244396, ... 
                      -0.298217036956, ... 
                      -0.108863472196, ...
                      0.0949327488264, ...
                      -0.00980469816509, ...
                      0.0000165167634970, ...
                      0.0000937359795772, ...
                      -0.0000123179218720e-5, ...
                      0.00196096504426];
                ih = [1, ...
                      1, ...
                      1, ...
                      2, ...
                      3, ...
                      3, ...
                      4, ...
                      5, ...
                      6, ...
                      7, ...
                      10];  
                jh = [0.25, ...
                      1, ...
                      2.5, ...
                      1.5, ...
                      1.5, ...
                      2.5, ...
                      2, ...
                      2, ...
                      5, ...
                      0.5, ...
                      10];
                % check for density
                if any(obj.state.rho_pure) == false
                    rho = obj.waterProp(T, P, 'density');
                    obj.state.rho_pure = rho;
                else
                    rho = obj.state.rho_pure;
                end
                rho_mole = rho./(obj.Mw_h2o*1e-3);
                % compute g value
                gpar = 1 + Nh(12).*(rho_mole./obj.rho_crit).*(T./228 - 1).^(-1.2);
                for ii=1:11
                    c1 = Nh(ii).*(rho_mole./obj.rho_crit).^(ih(ii)).*(obj.T_crit./T).^(jh(ii));
                    gpar = gpar + c1;
                end
                A = obj.Nv.*obj.muw^(2).*rho_mole.*gpar./(obj.eps_ref.*obj.k.*T);
                B = obj.Nv*obj.polar_a.*rho_mole./(3*obj.eps_ref);
                % composed parameters
                bpar = 1 + A + 5*B;
                cpar = 9 + 2*A + 18*B + A.^(2) + 10*A.*B + 9*B.^(2);
                apar = 4 - 4*B;
                epsilon = (bpar + sqrt(cpar))./apar;
                if any(obj.state.eps) == false
                    obj.state.eps = epsilon;
                end
                prop = epsilon; % dimensionless

            elseif strcmp('density', key) == true
                % convert Pa to MPa
                Zf = 1e-6;
                % compute density
                rho = 1./IAPWS_IF97('v_pT', P*Zf, T)';
                % add variable to the state
                if any(obj.state.rho_pure) == false
                    obj.state.rho_pure = rho;
                end
                prop = rho; % kg.m^-3

            elseif strcmp('heatCapacity', key) == true
                % convert Pa to MPa
                Zf = 1e-6;
                % enthalpy
                h = IAPWS_IF97('h_pT', P*Zf, T);
                % isobaric specific heat capacity
                cp = IAPWS_IF97('cp_ph', P*Zf, h)';
                % add value to the state
                if any(obj.state.cp_pure) == false
                    obj.state.cp_pure = cp;
                end
                prop = cp; % kJ.kg^(-1).K^(-1)
                
            elseif strcmp('vaporPressure', key) == true
                % convert MPa to Pa
                Zf = 1e6; 
                Psat = IAPWS_IF97('psat_T', T)'; % MPa
                % add value to the state
                if any(obj.state.Psat_pure) == false
                    obj.state.Psat_pure = Psat;
                end
                prop = Psat*Zf; % Pa

            elseif strcmp('Hvap', key) == true
                paramList = {'Hvap'};
                out = obj.readPar('water', paramList, 3);
                Cx = out;
                Hvap = Cx(1) + Cx(2)*T + Cx(3)*T.^(2); 
                % Add value to the state
                if any(obj.state.Hvap) == false
                    obj.state.Hvap = Hvap;
                end
                prop = Hvap; % J.kg^(-1)
            end
        end

        %% Mixture properties
        function mixProp = mixtureProp(obj, T, P, x, key)
        if strcmp('vaporPressure', key) == true
            if any(obj.state.Psat_pure) == false
                Psat_pure = obj.waterProp(T, P, 'vaporPressure');
            else
                Psat_pure = obj.state.Psat_pure;
            end
            % Margules contribution
            S = obj.Margules(T, x);
            Sw = S(1,:);
            % DH contribution
            DH = obj.DebyeHuckel(T, P, x);
            DHw = DH(1,:);
            % logarithmic activity for solvent
            aw_log = log(1 - x) + DHw + Sw;
            Psat = Psat_pure.*exp(aw_log);
            mixProp = Psat;
        elseif strcmp('density', key) == true
            paramList = {'density'};
            out = obj.readPar('electrolyte', paramList, 6);
            Cx = out;
            % convert K to C
            t = T - 273.15;
            % pure water density
            rho_pure = obj.state.rho_pure;
            % solution density
            test_rho = 1000.*ones(size(rho_pure));
            % molecular weight
            paramList = {'Mw_salt'};
            out = obj.readPar('electrolyte', paramList, 1);
            Mwx = out;
            % factor
            Zf = x/Mwx;
            % define function
            F = @(x) x - rho_pure - Cx(1).*Zf.*x - Cx(2).*Zf.*x.*t -  ...
                Cx(3).*Zf.*x.*t.^(2) - Cx(4).*(Zf.*x).^(3/2) - ...
                Cx(5).*(Zf.*x).^(3/2).*t - Cx(6).*(Zf.*x).^(3/2).*t.^(2);
           % solve
           rho = fsolve(F, test_rho);
           mixProp = rho;

        elseif strcmp('conductivity', key) == true
            paramList = {'conductivity'};
            out = obj.readPar('electrolyte', paramList, 6);
            Cx = out;
            % compute density if it has not been calculated
            if any(obj.state.rho) == false
                rho = obj.mixtureProp(T, P, x, 'density');
            else
                rho = obj.state.rho;
            end
            % salt molecular weight
            paramList = {'Mw_salt'};
            Mwx = obj.readPar('electrolyte', paramList, 1);
            % general coefficient
            M = rho.*x/(Mwx);
            % conductivty
            kappa = Cx(1)*M + Cx(2)*M.^(2) + Cx(3).*M.*T + Cx(4)*M./T + ...
                Cx(5)*M.^(3) + Cx(6)*M.^(2).*T.^(2);
            mixProp = kappa*1e2;% S.m^(-1)

        elseif strcmp('conductivityCleggPitzer', key) == true
            if any(obj.state.rho) == false
                rho = obj.mixtureProp(T, P, x, 'density');
            else
                rho = obj.state.rho;
            end
            % Margules contribution
            S = obj.Margules(T, x);
            Sw = S(2,:);
            % DH contribution
            DH = obj.DebyeHuckel(T, P, x);
            DHw = DH(2,:);
            % logarithmic activity for solvent
            Gamm_log = DHw + Sw;
            Gamm_X = exp(Gamm_log);
            
            mu_ = obj.mixtureProp(T, P, x, 'viscosity');
            u = obj.z.*obj.e./(6*pi*mu_*obj.aR); % ionic mobility
            % molecular weight
            paramList = {'Mw_salt'};
            Mwx = obj.readPar('electrolyte', paramList, 1);
            cX = rho.*x/(Mwx*1e-3); % mol/m^3
            kappa = sum(u.*obj.z.*obj.Fc)*Gamm_X;
            mixProp = kappa;


        elseif strcmp('eqPotential', key) == true
            paramList = {'potential'};
            out = obj.readPar('electrolyte', paramList, 4);
            Cx = out;
            Eref = Cx(1) + Cx(2)*T + Cx(3)*T.*log(T) + Cx(4)*T.^(2);
            % convert Pa to bar
            Zf = 1e-5;
            % vapor pressure of pure water in bar
            if any(obj.state.Psat_pure) == false
                Psat_pure = obj.waterProp(T, P, 'vaporPressure');
                obj.state.Psat_pure = Psat_pure;
            else
                Psat_pure = obj.state.Psat_pure;
            end
            % vapor pressure of the mixture
            if any(obj.state.Psat) == false
                Psat = obj.mixtureProp(T, P, x, 'vaporPressure');
                obj.state.Psat = Psat;
            else
                Psat = obj.state.Psat;
            end
            % absolute pressure
            Pabs = P*Zf - 1.013;
            % equilibrium potentoal
            c1 = ((Pabs - (Psat*Zf - 1.013)).^(1.5))./((Psat*Zf - 1.013)./(Psat_pure*Zf - 1.013));
            En = Eref + (obj.R*T/(2*obj.Fc)).*log(c1);
            mixProp = En;
        
        elseif strcmp('heatCapacity', key) == true
            % compute pure water density if it has not been computed
            if any(obj.state.rho_pure) == false
                rho_pure = obj.waterProp(T, P, 'density');
                obj.state.rho_pure = rho_pure;
            else
                rho_pure = obj.state.rho_pure;
            end
            % compute isobaric specific heat capacity
            if any(obj.state.cp_pure) == false
                cp_pure = obj.waterProp(T, P, 'heatCapacity');
                obj.state.cp_pure = cp_pure;
            else
                cp_pure = obj.state.cp_pure;
            end
            % read molecular weight of salt
            paramList = {'Mw_salt'};
            out = obj.readPar('electrolyte', paramList, 1);
            Mwx = out;
            % Compute molality
            m = x./((1 - x)*Mwx*1e-3); % mol.kg^(-1)
            % Redlich Meyer term
            B = obj.RedlichMeyer(T, rho_pure, m);
            % Debye Huckel coefficient
            Ac = obj.funcA(T, P, 'heatCapacity');
            % heat capacity
            c1 = (m.*Ac.*(m.^(1/2) + B) + cp_pure); 
            c2 = (Mwx*1e-3.*m + 1);
            cp = c1./c2; % J.K^(-1).kg^(-1)
            mixProp = cp;
        
        elseif strcmp('thermoneutral', key) == true
            % convert K to C
            t = T - 273.15;
            % get empirical coefficients
            paramList = {'Etn'};
            out = obj.readPar('electrolyte', paramList, 4);
            Cx = out;
            % convert Pa to bar
            f1 = 1e-5;
            % pure water vapor enthalpy
            if any(obj.state.Hvap) == false
                Hvap = obj.waterProp(T, P, 'Hvap');
                obj.state.Hvap = Hvap;
            else
                Hvap = obj.state.Hvap;
            end
            % vapor pressure of mixture
            if any(obj.state.Psat) == false
                Psat = obj.mixtureProp(T, P, x, 'vaporPressure')*f1;
                obj.state.Psat = Psat;
            else
                Psat = obj.state.Psat;
            end
            
            % absolute pressure
            Pabs = P*f1 - 1.013;
            % compute thermoneutral potential
            Etn = Cx(1) + Cx(2)*t + Cx(3)*t.^(2) + ...
                Cx(4)*Psat./(Pabs - Psat).*Hvap/(2*obj.Fc);
            mixProp = Etn;
        
        elseif strcmp('surfTension', key) == true
            t = T - 273.15;
            n = 5;
            sigma = 0.;
            for i=1:n
                for j=1:n
                    sigma = sigma + obj.Asigma(i,j)*t.^(j-1)*x.^(i-1)*1e-3;
                end
            end
            mixProp = sigma; % N.m^(-1)
        
        elseif strcmp('viscosity', key) == true
            n = 5;
            mu = 0;
            for i=1:n
                mu = mu + obj.muVector(i)*T.^(i-1);
            end
            mixProp = mu; % Pa.s
        end        
        end

        %% Initialize state
        function output = initState(obj, T, P, x)
            % convert bar to Pa
            P = P*1e5; % Pa
            % convert °C to K
            T = T + 273.15;
            % set temperature
            obj.state.T = T; % K
            % set pressure
            obj.state.P = P; % Pa
            % set mass fraction of salt
            obj.state.x = x;
            % pure water density
            rho_pure = obj.waterProp(T, P, 'density');
            obj.state.rho_pure = rho_pure;
            % mixture density
            rho = obj.mixtureProp(T, P, x, 'density');
            obj.state.rho = rho;
            % dielectric constant pure water
            eps = obj.waterProp(T, P, 'epsilon');
            obj.state.eps = eps;
            % saturation pressure pure water
            Psat_pure = obj.waterProp(T, P, 'vaporPressure');
            obj.state.Psat_pure = Psat_pure;
            % heat capcacity of pure water
            cp_pure = obj.waterProp(T, P, 'heatCapacity');
            obj.state.cp_pure = cp_pure;
            % enthalpy of vaporization of pure water
            Hvap = obj.waterProp(T, P, 'Hvap');
            obj.state.Hvap = Hvap;
            % saturation pressure of mixture
            Psat = obj.mixtureProp(T, P, x, 'vaporPressure');
            obj.state.Psat = Psat;
            % heat capacity of mixture
            cp = obj.mixtureProp(T, P, x, 'heatCapacity');
            obj.state.cp = cp;
            % conductivity
            kappa = obj.mixtureProp(T, P, x, 'conductivity');
            obj.state.kappa = kappa;
            % equilibrium potential
            Eref = obj.mixtureProp(T, P, x, 'eqPotential');
            obj.state.Eref = Eref;
            % thermoneutral potential
            Etn = obj.mixtureProp(T, P, x, 'thermoneutral');
            obj.state.Etn = Etn;
            % surface tension
            sigma = obj.mixtureProp(T, P, x, 'surfTension');
            obj.state.sigma = sigma;
            % dynamic viscosity
            mu = obj.mixtureProp(T, P, x, 'viscosity');
            obj.state.mu = mu;
            % conductivity from Clegg and Pitzer
            kappaClPtz = obj.mixtureProp(T, P, x, 'conductivityCleggPitzer');
            obj.state.kappaClPtz = kappaClPtz;

            % save modifications in the object
            output = obj;
        end
        % Gas molar volume
        function output = gasMolarVol(obj, T, P)
            output = obj.R*T/P;
        end
    end % methods
end
        