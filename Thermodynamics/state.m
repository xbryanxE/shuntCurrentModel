classdef state
    properties
    % pure substance density
    rho_pure = 0;
    % density of mixture
    rho = 0;
    % static dielectric constant
    eps = 0;
    % pure substance vapor pressure
    Psat_pure = 0;
    % vapor pressure as mixture
    Psat = 0;
    % Pressure
    P = 0;
    % Temperature
    T = 0;
    % mass fraction
    x = 0;
    % isobaric heat capacity of pure substance
    cp_pure = 0;
    % isobaric heat capacity of mixture
    cp = 0;
    % enthalpy of vaporization
    Hvap = 0;
    % conductivity of mixture
    kappa = 0;
    % conductivity from Clegg and Pitzer
    kappaClPtz = 0;
    % equilibrium potential
    Eref = 0;
    % thermoneutral potential
    Etn = 0;
    % surface tension
    sigma = 0;
    % dynamic viscosity
    mu = 0;
    end
end