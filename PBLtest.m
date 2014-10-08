%PBL TEAM 13235 test
%Show up please!!
%hi teja
%hi megan

function PBLtest
gasexchange2
entry
blood
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

%PV=nRT
%R=62.36 (L*mmHg)(mol*K)

% BTP values assumed for function deadspace and function gas_exchange_1
function values_ds= deadspace(n2,nCO2,nN2,nH2O)
DSV=0.30*TV;   %30% of the total inhaled air goes to dead space
VO2=(nO2*62.36*310)/149.3;  %partial volume of O2 during both inhalation and exhalation
VCO2=(nCO2*62.36*310)/0.3;  %partial volume of CO2 during both inhalation and exhalation
VN2=(nN2*62.36*310)/563.4;  %partial volume of N2 during both inhalation and exhalation
VH2O=(nH2O*62.36*310)/47;   %partial volume of H2O during both inhalation and exhalation
values_ds=[VO2 VCO2 VN2 VH2O];  %partial volumes for deadspace unit
A: m3-m4=0;  %total mass flow rate equation
m3w3,CO2-m4w4,CO2=0;  %mass flow rate equation for CO2
m3w3,O2-m4w4,O2=0; %mass flow rate equation for O2
m3w3,N2-m4w4,N2=0; %mass flow rate equation for N2
m3w3,H2O-m4w4,H2O=0; %mass flow rate equation for H2O
end

function values_ge1= gas_exchange_1(nO2_i,nCO2_i,nO2_f,nCO2_f,nN2,nH2O)
ge1v=0.70*tv; %70% of the total inhaled air goes to the gas exchange 1 box
VO2_i=(nO2_i*62.36*310)/149.3;  %partial volume of O2 during inhalation
VCO2_i=(nCO2_i*62.36*310)/0.3; %partial volume of CO2 during inhalation
VO2_f=(nO2_f*62.36*310)/120; %partial volume of O2 during exhalation
VCO2_f=(nCO2_f*62.36*310)/27;  %partial volume of CO2 during exhalation
VN2_i=(nN2*62.36*310)/563.4; %partial volume of N2 during exhalation
VN2_f=(nN2*62.36*310)/566; %partial volume of N2 during exhalation
VH2O=(nH2O*62.36*310)/47;  %partial volume of H2O during both inhalation and exhalation
values_ge1=[VO2_i VCO2_i VO2_f VCO2_f VN2_i VN2_f VH2O];  %partial volumes for gas exchange 1 unit
Total: m5+m8+m10-m6-m7-m9=0;  %total mass flow rate equation
m10w10,CO2+m8w8,CO2-m6w6,CO2=0;  %mass flow rate equation for CO2
m5w5,O2-m9w9,O2-m7w7,O2=0;  %mass flow rate equation for O2
m5w5,N2-m6w6,N2=0;  %mass flow rate equation for N2
m5w5,H2O-m6w6,H2O=0;  %mass flow rate equation for H2O
end

%only for steady state
function blood
vO2inGE1 = 21; %ml/min/mm Hg
%vO2inGE2 = ;

dO2 =  density(1, 31.9988, 288); %1 atm, 31.9988 g/mol, 288 K
dCO2 = density(1, 44.0095, 288); %1 atm, 44.0095 g/mol, 288 K

[mCO2outGE1, mO2inGE1] = bloodGE1(vCO2outGE1, vO2inGE1, dO2, dCO2) ;
[mCO2outGE2, mO2inGE2] = bloodGE2(vCO2outGE2, vO2inGE2, dO2, dCO2);
[mCO2inSurroundings, mO2outSurroundings] = bloodSurroundings(vCO2inSurr, vO2outSurr, dO2, dCO2);
 
mCO2in = mCO2inSurroundings; %adding up all the mass flow rates of streams bringing CO2 in
mCO2out = mCO2outGE1 + mCO2outGE2; %adding up all the mass flow rates of streams bringing CO2 out
mO2in = mO2inGE1 + mO2inGE2; %adding up all the mass flow rates of streams bringing O2 in
mO2out = mO2outSurroundings; %adding up all the mass flor rates of streams bringing O2 out

blood = [mO2in mCO2in mO2out mCO2out]; %blood box = mass flow rate 
return

function [mCO2in, mO2out] = bloodGE1(vCO2, vO2, dO2, dCO2)
mCO2in = v2m(vCO2 , dCO2);
mO2out = v2m(vO2 , dO2);
end

function [mCO2in, mO2out] = bloodGE2(vCO2, vO2, dO2, dCO2)    
mCO2in = v2m(vCO2 , dCO2);
mO2out = v2m(vO2 , dO2);
end

function [mCO2out, mO2in] = bloodSurroundings(vCO2, vO2, dO2, dCO2)
mCO2out = v2m(vCO2 , dCO2);
mO2in = v2m(vO2 , dO2);
end

%this function converts volumetric flow rate to mass flor rate
function massflowrate = v2m(volumetricflowrate, density)
massflowrate = volumetricflowrate*density;
end 