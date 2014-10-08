%PBL TEAM 13235 test
%Show up please!!
%hi teja
%hi megan
<<<<<<< HEAD
<<<<<<< HEAD
%hi tony!
%teja_10/6_test
=======
<<<<<<< HEAD
=======
% <<<<<<< HEAD
>>>>>>> 3fad6573fd84ffb1fcd1b5df790ae057b54dd72a

function PBLtest
gasexchange2
entry
end

function entry
vper1 = [vper1O2 vper1CO2 vper1N2 vper1H2O vper1He vper1Ne vper1Ar vper1Kr vper1Xe];
vper2 = [vper2O2 vper2CO2 vper2N2 vper2H2O vper2He vper2Ne vper2Ar vper2Kr vper2Xe];
vper4 = [vper4O2 vper4CO2 vper4N2 vper4H2O vper4He vper4Ne vper4Ar vper4Kr vper4Xe];
vper5 = [vper5O2 vper5CO2 vper5N2 vper5H2O vper5He vper5Ne vper5Ar vper5Kr vper5Xe];
vper6 = [vper6O2 vper6CO2 vper6N2 vper6H2O vper6He vper6Ne vper6Ar vper6Kr vper6Xe];
vper7 = [vper7O2 vper7CO2 vper7N2 vper7H2O vper7He vper7Ne vper7Ar vper7Kr vper7Xe];
vper = [vper1; vper2; vper4; vper5; vper6; vper7];
% vper = volume percentages for each constitunet in each stream
% Note: vper values for streams 4, 5, and 6 will be equal
% vper7 based on gas exchange 1 and 2

M1 = [31.9988 44.0095 28.01348 33.00674 4.006202 20.1797 39.948 83.8 ...
    131.29];
M = [M1; M1; M1; M1; M1; M1];
% M = molar mass of each constituent in each stream

w = massfrac(vper,M);
% w = the mass fraction of each constituent in each stream

mass1 = totalmass(0.500,vper1,M);
% mass1 = the total mass of inspired air
m1 = mass1 / t_insp;
% t_insp = time of inspiration
% m1 = the mass flow rate of inspired air

deadfrac = 0.30;
m4,m5,m6 = deadspace(m1,deadfrac);
% calculates mass flow rate of streams 4, 5, and 6

m2 = m7 + m5;
% m2 depends on m7
% m7 depends on gas exchange 1 and 2
% if no gas exchange during expiration, m7 = m9 which is calculated in 
% gasexchange2

m = [m1; m2; m3; m4; m5; m6; m7];
% m = mass flow rates of all mass streams
% units? check all kg vs g especially

ms = constituentflow(m,w);
% calculates mass flow rates of constituents of a stream

% calculate m3
% can do so based on m1H2O and m6H2O

c = 1005;
Tb = 310; % body temperature
Ta = 288; % inspired air temperature
Q16 = thermal(mass1,c,Ta,Tb);
Q17 = thermal(mass1,c,Tb,Ta);
% calculates transfer of thermal energy in streams 16 and 17
end

% thermal calculates thermal energy flow rate
% Q = amount of thermal energy tranferred in J
% m = mass of inhaled/exhaled air in kg
% c = the specific heat capacity of air in J/kg/K
% Ti = the temperature of one object involved in heat transfer
% Tf = the temperature of the other object involved in heat transfer
% dt = the change in temperature

function Q = thermal(mass,c,Ti,Tf)
dT = Tf - Ti;
Q = mass * c * dT;
end

% massfrac calculates the mass fraction of constituents in a stream
% vper = the volume percentage of the constituents
% vfrac = the volume fraction of the constituents
% M = the molar mass of the constituents
% mratio = the mass ratios of the constituents
% sum_ratio = the sum of the mass ratios
% w = the mass fractions of the constituents

function w = massfrac(vper,M)
vfrac = vper ./ 100;
% calculates volume fractions from volume percentages
mratio = M .* vfrac;
% calculates mass ratios from volume fraction and molar mass
sum_ratio = sum(mratio);
% sum of mass ratios
w = mratio ./ sum_ratio;
% mass fractions calculated from mass ratios and sum
end

% constituent flow calculates the mass flow rate of each constituent in a
% stream
% m = total mass of the stream
% w = mass fractions of constituents

function ms = constituentflow(m, w)
ms = m .* w;
end

% deadspace calculates the mass flow rate of streams 4, 5, and 6
% "Entry" box is a splitter
% calculates mass flow rate based on fraction of air that goes to dead
% space

function [m4,m5,m6] = deadspace(m1,deadfrac)
m4 = m1 * deadfrac; % fraction of inlet air to deadspace
m5 = m4;
m6 = m1 * (1 - deadfrac);
end

% total mass finds the total mass of inspired air
% vtot = total volume of inspired air
% pp = partial pressure of a constituent in inspired air
% Ta = temperature of inspired air
% R = universal gas constant 
% n = moles of constituent in inspired air
% species_mass = the mass of each species in inspired air
% mass1 = the total mass of inspired air

function mass1 = totalmass(vtot,vper,M)
vfrac = vper ./ 100;
v = vfrac * vtot;
pp = 760 * vfrac;
Ta = 288;
R = 62.3637;
n = (pp .* v) / (R * Ta);
species_mass = n .* M;
mass1 = sum(species_mass);
end 

function gasexchange2
% composition of stream 8 will be calculated by gas exchange 1
% mass flow rates of stream 8 will be calculated by gas exchange 1
PO2alv = 104; %mmHg
PCO2alv = 46; %mmHg
MO2 = 31.9988;
MCO2 = 44.0095;
Tb = 310;
dO2alv = density(PO2alv,MO2,Tb);
dCO2alv = density(PCO2alv,MCO2,Tb);
% calculates densities of O2 and CO2 in alveoli <-- is this constant during
% transport across membrane??

% v14 =
% v15 =
m14 = v14 * dO2alv;
m15 = v15 * dCO2alv;
% calculate mass flow rates from volumetric flow rates and densities

m9O2 = m8O2 - m14;
m9CO2 = m8CO2 + m15;
% calculate mass flow rates of O2 and CO2 in stream 9
% m8O2 and m8CO2 calculated in gas exchange 1

m9 = m8 - m14 + m15;
% mass flow rate of stream 9
% mass flow rates of N2, H20, inert gases remains the same

w9O2 = m9O2 / m9;
w9CO2 = m9CO2 / m9;
% mass fractions of O2 and CO2 in stream 9

end

% density calculates the density of a gas
% P = pressure
% M = molar mass
% T = temperature
% R = universal gas constant
% d = density 

function d = density(P,M,T)
R = 62.3637; % if units chosen are mmHg, L, and K
d = (M * P) / (T * R);
end


=======
%hi tony!
>>>>>>> a97a06baaa8273de4bda938cba276dd71daffa4a
>>>>>>> 5570270266d8e750bdca3fcea7c461f4f2d3928f
