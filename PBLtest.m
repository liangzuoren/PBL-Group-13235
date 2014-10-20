
function PBLtest
entry
% blood
end

%Entry box is where inhaled air is humidified and exhaled air leaves the
%body from
function entry
M = [31.9988 44.0095 28.01348 33.00674];
% M = molar mass of each constituent
TV = 0.5; % liters
H = 1.73; % height in meters of standard man
W = 68; % weight in kilograms of standard man

PP2 = [116.0 32.0 565.0 47.0];
% partial pressures in expired air
PP1 = [158.0 0.3 596.0 5.7]; 
% partial pressures in inspired air
PP7 = [100.0 40.0 573.0 47.0];
% partial pressures in alveolar air

vfrac2 = PP2 ./ 760;
vfrac1 = PP1 ./ 760;
vfrac7 = PP7 ./ 760;
% calculates the volume fractions of the constituents from the partial
% pressures
% vper4 = humid(vper_1,0.5,50);
% calculates composition after humidification in Entry unit
% vfrac6 = vfrac4;
% the compositions of stream 6 and 4 are equal because stream 4 is
% the only inlet stream and stream 6 is the only outlet stream of the dead 
% space unit, and there are no reactions
[RFin,RFex] = RF(H,W);
% RF calculates the respiration frequencies for inspiration and expiration 
% based on height and weight
tin = 30/RFin;
% RFin is number of breaths per minute, divide by 2 to get number of
% inspirations per minute and then multiply by 60s/1min to get time of
% inspiration
texp = 30/RFex;
% RFex is number of breaths per minute, divide by 2 to get number of
% expirations per minute and then multiply by 60s/1min to get time of
% expiration
tresp = tin + texp; % + t_bh
% sum of inspiration, breathold, and expiration times is the time of one 
% full respiraory cycle

insp_range = 0:0.01:tin;
exp_range = 0:0.01:texp;
resp_range = 0:0.01:tresp;
index_in = length(insp_range);
index_exp = length(exp_range);
index_resp = length(resp_range);
for i = 1:index_in
    vflow1_in(i) = -volumetricflow(RFin,TV,insp_range(i));
    % volumetric flow rate for inspiration in L/s
    % airflow into body defined as negative direction
    
    vflow4_in(i) = vflow1_in(i);
    vflow6_in(i) = 0.7 * vflow1_in(i);
    % entry unit is splitter, 30% of flow to Dead Space and 70% to Gas Exchange

    % calculate vflow3_in somehow...
     
    vflow2_in(i) = 0;
    vflow5_in(i) = 0;
    vflow7_in(i) = 0;
    %zero flow rate in these streams during inspiration
    
    Pav_in(i) = alveolarpressure(vflow1_in(i));
    % calculates alveolar pressure with respect to time during inspiration
    
end
for i = 1:index_exp
    vflow2_ex(i) = volumetricflow(RFex,TV,exp_range(i));
    % volumetric flow rate for expiration in L/s
    x = deadspace_expfrac(TV,RFex,texp);
    vflow5_ex(i) = vflow2_ex(i);
    vflow7_ex(i) = vflow5_ex(i) - x*vflow2_ex(i);
    
    Pav_ex(i) = alveolarpressure(vflow2_ex(i));
    % calculates alveolar pressure with respect to time during expiration
    
    vflow1_ex(i) = 0;
    vflow4_ex(i) = 0;
    vflow6_ex(i) = 0;
    % vflow3_ex = 0;
    %zero flow rate in these streams during expiration
    
end

Pav = [Pav_in Pav_ex];

plot(resp_range,Pav)
title('Alveolar Pressure Over One Respiratory Cycle')
xlabel('Time (s)')
ylabel('Alveolar Pressure (mmHg)')
% graph of alveolar pressure over one full respiratory cycle

vflow1 = [vflow1_in vflow1_ex];
vflow2 = [vflow2_in vflow2_ex];
% vflow3 = [vflow3_in vflow3_ex];
vflow4 = [vflow4_in vflow4_ex];
vflow5 = [vflow5_in vflow5_ex];
vflow6 = [vflow6_in vflow6_ex];
vflow7 = [vflow7_in vflow7_ex];
% combines functions for inspiration and expiration

% calculate Pav_O2 and Pav_CO2 from Pav


for i = 1:index_resp
    vflow8(i) = 21*60/100*Pav(i);
    % calculates O2 transport with respect to time
    vflow9(i) = 425*60/100*Pav(i);
    % calculates CO2 transport with respect to time
    vflow10(i) = vflow8(i);
    vflow11(i) = vflow9(i);
    % solves volumetric flow conservation equations for blood unit to
    % calculate O2 and CO2 transport rate out of the blood to tissues
end

vflow = [vflow1; vflow2; vflow4; vflow5; vflow6; vflow7; vflow8; ...
    vflow9; vflow10; vflow11];

t_start_4 = traveltime(RFin,TV,index_in,insp_range)

figure
plot(resp_range,vflow1,resp_range,vflow2,resp_range,vflow4,resp_range,...
    vflow5,resp_range,vflow6,resp_range,vflow7)
title('Volumetric Flow Rates of Streams 1, 2, 4, 5, 6, and 7')
xlabel('Time (s)')
ylabel('Volumetric Flow Rate (L/s)')
% figure
% plot(resp_range,vflow8,resp_range,vflow10,resp_range,vflow9, ...
    % resp_range,vflow11)
% title('Volumetric Flow Rate for O2 and CO2 Diffusion')



% w = massfrac(vper,M);
% % w = the mass fraction of each constituent in each stream

% mass1 = totalmass(TV,vper1,M);
% mass1 = the total mass of inspired air

% vs = constituent_volume(V,w);
% % calculates volume of constituents in a unit

% calculate m3
% can do so based on m1H2O and m6H2O

% c = 1005;
% Tb = 310; % body temperature
% Ta = 288; % inspired air temperature
% Q12 = thermal(mass1,c,Ta,Tb)
% Q13 = thermal(mass1,c,Tb,Ta)
% % calculates transfer of thermal energy in streams 12 and 13
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


function vs = constituent_volume(V, w)
vs = V .* w;
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
v = vfrac .* vtot;
pp = 760 .* vfrac;
Ta = 288;
R = 62.3637;
n = (pp .* v) / (R * Ta);
species_mass = n .* M;
mass1 = sum(species_mass);
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
%{
%PV=nRT
%R=62.36 (L*mmHg)(mol*K)

<<<<<<< HEAD
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
%}

% % BTP values assumed for function deadspace and function gas_exchange_1
% function values_ds= deadspace(n2,nCO2,nN2,nH2O)
% DSV=0.30*TV;   %30% of the total inhaled air goes to dead space
% VO2=(nO2*62.36*310)/149.3;  %partial volume of O2 during both inhalation and exhalation
% VCO2=(nCO2*62.36*310)/0.3;  %partial volume of CO2 during both inhalation and exhalation
% VN2=(nN2*62.36*310)/563.4;  %partial volume of N2 during both inhalation and exhalation
% VH2O=(nH2O*62.36*310)/47;   %partial volume of H2O during both inhalation and exhalation
% values_ds=[VO2 VCO2 VN2 VH2O];  %partial volumes for deadspace unit
% A: m3-m4=0;  %total mass flow rate equation
% m3w3,CO2-m4w4,CO2=0;  %mass flow rate equation for CO2
% m3w3,O2-m4w4,O2=0; %mass flow rate equation for O2
% m3w3,N2-m4w4,N2=0; %mass flow rate equation for N2
% m3w3,H2O-m4w4,H2O=0; %mass flow rate equation for H2O
% end

% function values_ge1= gas_exchange_1(nO2_i,nCO2_i,nO2_f,nCO2_f,nN2,nH2O)
% ge1v=0.70*tv; %70% of the total inhaled air goes to the gas exchange 1 box
% VO2_i=(nO2_i*62.36*310)/149.3;  %partial volume of O2 during inhalation
% VCO2_i=(nCO2_i*62.36*310)/0.3; %partial volume of CO2 during inhalation
% VO2_f=(nO2_f*62.36*310)/120; %partial volume of O2 during exhalation
% VCO2_f=(nCO2_f*62.36*310)/27;  %partial volume of CO2 during exhalation
% VN2_i=(nN2*62.36*310)/563.4; %partial volume of N2 during exhalation
% VN2_f=(nN2*62.36*310)/566; %partial volume of N2 during exhalation
% VH2O=(nH2O*62.36*310)/47;  %partial volume of H2O during both inhalation and exhalation
% values_ge1=[VO2_i VCO2_i VO2_f VCO2_f VN2_i VN2_f VH2O];  %partial volumes for gas exchange 1 unit
% Total: m5+m8+m10-m6-m7-m9=0;  %total mass flow rate equation
% m10w10,CO2+m8w8,CO2-m6w6,CO2=0;  %mass flow rate equation for CO2
% m5w5,O2-m9w9,O2-m7w7,O2=0;  %mass flow rate equation for O2
% m5w5,N2-m6w6,N2=0;  %mass flow rate equation for N2
% m5w5,H2O-m6w6,H2O=0;  %mass flow rate equation for H2O
% end


%only for steady state
function blood
%vO2inGE2 = 21; %ml/min/mm Hg
%vO2inGE1 = vO2inGE2/0.7*0.3;
%vCO2outGE2 = 451;
%vCO2outGE1 = vCO2outGE2/0.7*0.3;
%vO2outblood = vO2inGE1 + vO2inGE2;
%vCO2inblood = vCO2outGE1 + vCO2outGE2;
tin = 30/RFin;
% RFin is number of breaths per minute, divide by 2 to get number of
% inspirations per minute and then multiply by 60s/1min to get time of
% inspiration
texp = 30/RFex;
% RFex is number of breaths per minute, divide by 2 to get number of
% exspirations per minute and then multiply by 60s/1min to get time of
% exspiration
tresp = tin + texp;
% calculates length of time per inspiration and expiration from breathing
% frequencies
PaO2i = 100;%100 mmHg partial pressure in beginning for O2
PaCO2i = 40;%40 mmHg partial pressure in beginning for CO2
Patm = 760; 
PH2O = 47;%presures taken from diagram in assumed values
PN2 = 573;

DifCapO2 = 21; %mL/min/mmHg
DifCapCO2 = 400; %mL/min/mmHg

dO2 =  density(1, 31.9988, 310); %1 atm, 31.9988 g/mol, 310 K(body temperature)
dCO2 = density(1, 44.0095, 310); %1 atm, 44.0095 g/mol, 310 K(body temperature)


for x=0:1 %modeling 1 breath
    t=0;
    PaO2alveoli = PaO2i;
    PaCO2alveoli = PaCO2i;
    trange = 1;
    tstep = 0.01;

    while t>tin;
        
        difRateO2 = DifCapO2 * (PaO2alveoli - PaO2capillary);
        difRateCO2 = DifCapCO2 * (PaCO2alveoli - PaCO2capillary);
        mO2 = dO2*difRateO2;
        mCO2 = dCO2*difRateCO2;
        nO2 = mO2/31.9988;
        nCO2 = mCO2/44.0095;
        %PV=nRT, P = nRT/V
        PaO2diff = nO2*tstep*62.36367*310/3; %gives pressure in mmHg
        PaCO2diff = nCO2*tstep*62.36367*310/3;
        PaO2alveoli = PaO2alveoli - PaO2diff;
        PaCO2alveoli = PaCO2alveoli + PaCO2diff;
        PaO2capillary = PaO2capillary - PaO2diff;
        PaCO2capillary = PaCO2capillary - PaCO2diff;
        
        %mTotal = residualVolume + vflow6_in(trange)*dInspiredAir; %volumetric flow rate * density?
        %PaO2alveoli = (PaO2alveoli * mTotal- mO2 )/ mTotal;%PaO2 = FiO2*(Patm - PH2O) - (PaCO2/RQ)
        %PaO2capillary = 40;
        %PaCO2alveoli = (PaCO2alveoli * mTotal - mCO2 )/ mTotal;
        %PaCO2capillary = 46;
        
        t=t+tstep;
        trange = trange+1;
    end
    
    trange = 1;
    
    while t>texp
        
        difRateO2 = DifCapO2 * (PaO2alveoli - PaO2capillary);
        difRateCO2 = DifCapCO2 * (PaCO2alveoli - PaCO2capillary);
        mO2 = dO2*difRateO2;
        mCO2 = dCO2*difRateCO2;
        %mTotal = residualVolume + vflow7_ex(trange)*dAir;
        PaO2alveoli = (PaO2alveoli * mTotal- mO2 )/ mTotal;
        PaO2capillary = 40;
        PaCO2alveoli = (PaCO2alveoli * mTotal - mCO2 )/ mTotal;
        PaCO2capillary = 46;
        
        t=t+tstep;
        trange = trange+1;
    end
    
end


%{
[mCO2outGE1, mO2inGE1] = bloodGE1(vCO2outGE1, vO2inGE1, dO2, dCO2) ;
[mCO2outGE2, mO2inGE2] = bloodGE2(vCO2outGE2, vO2inGE2, dO2, dCO2);
[mCO2inSurroundings, mO2outSurroundings] = bloodSurroundings(vCO2inSurr, vO2outSurr, dO2, dCO2);
 
mCO2in = mCO2inSurroundings; %adding up all the mass flow rates of streams bringing CO2 in
mCO2out = mCO2outGE1 + mCO2outGE2; %adding up all the mass flow rates of streams bringing CO2 out
mO2in = mO2inGE1 + mO2inGE2; %adding up all the mass flow rates of streams bringing O2 in
mO2out = mO2outSurroundings; %adding up all the mass flow rates of streams bringing O2 out

mBlood = [mO2in mCO2in mO2out mCO2out]; %blood box = mass flow rate 

end



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

%this function converts volumetric flow rate to mass flow rate
function massflowrate = v2m(volumetricflowrate, density)
massflowrate = volumetricflowrate*density;
%}
end

%calcuate mass flow rateswith regards to time
%model diffusion over time
%partial pressure change over time
%
%
%humid calculates volume percentage in humidfied air

%humidity function calculates the volume percentages of the air constiuents


%humid function calculates the volume percentages of the air constiuents

%after is has been humidified to 100%
%h_i= initial relative humidity of air when inspired
%m_i=initial amount of water in air when inspired (mg)
%m_f=final amount of water in air after humdidification (mg)
%water_added=total amount of water added to air when it is humidified (mg)
%mass_tot1=total mass of gas before humidified (mg)
%mass_tot2= total mass of gas after humidified (mg)
%w_i= initial mass fractions of gas constiuents before humidified
%mass_i=amount of each consituent in inhaled air (mg)
%mass_f=amount of each consituent in humidified air (mg)
%w_f= mass fractions of gas constituents after humidified
%vper_1=volume percentages of air in stream before it is humidified
%vper_h= volume percentages of air after it is humidified
%vtot= volume of air inhaled- 0.5 L in this model

function vper_h=humid(h_i,vtot,vper1,M)
m_i=(0.18*h_i)*vtot;
%at 24 degrees Celsius, amount of water in air before humidified (mg), at 50% humidity
m_f=22;
%at 37 degrees Celsius, amount of water in air after humidified (mg), at 100% humidity
water_added=m_f-m_i;
density=1184; 
%density of inhaled air (mg/L)
mass_tot1=density*vtot;
%total mass of inhaled air (mg)
mass_tot2=mass_tot1+water_added;
%total mass of humidifed air (mg)
w_i=massfrac(vper1,M);   
%finds mass fractions of all constituents of inhaled air in Stream 1
mass_i=w_i.*mass_tot1;  
%finds amount of each consituent in inhaled air (mg)
mass_f=mass_i;
%creates variable to hold mass amounts of all consituents after humidified
mass_f(4)=mass_f(4)+ water_added;  
%amount of water in inspired gas after it has been humidified (mg)
w_f= mass_f ./ mass_tot2; 
%finds mass fractions of all constituents of air after humidified
vper_h=volumepercentage(mass_f,M);
end

%volumepercentage function calculates the volume percentages of constituents in a stream 
%based on the amounts of each constituent (in mg)
%mass_f=mass of all constitiuents (mg)
%M= molar masses of all constituents (g/mol)
%v_fracf=volumetric fractions of all gas constituents after humidified
%vper_f= volume percentages of all gas constituents after humidified
function vper_f=volumepercentage(mass_f,M)
nratio=(mass_f/1000)./M;  
%calculates number of moles of all gas constituents after humidified
nratio_sum= sum(nratio);  
%calculates sum of all moles of gas after humidified
v_fracf= nratio./nratio_sum;  
%calculates volumetric fractions (which are equal to mole fractions)
vper_f=vfrac_f.*100;  
%finds volume percentages of gas after humidified
end

% calculates respiration frequency in breaths per minute
% H = height in meters
% W = weight in kilograms
% RFin = respiration frequency for inspiration
% RFex = respiration frequency for expiration

function [RFin, RFex] = RF(H,W)
RFin = 55.55 - 32.86*H + 0.2602*W;
RFex = 77.03 - 45.42*H + 0.2373*W;
end

% calculates volumetric flow rate of inspiration or expiration
% RF = breathing frequency
% TV = tidal volume
% t = time

function vflow = volumetricflow(RF,TV,t)
vflow = (pi*RF*TV)/60 * sin(pi*RF*t/30);
end

% calculates the pressure in the alveoli with respect to time
% Raw = airway resistance
% Pb = baromatic pressure
% vflow = the volumetric flow rate of air (either inspiratory or expiratory
% flow rate, depending on the time)

function Pav = alveolarpressure(vflow1)
Raw = 1.41 / 1.36; % 1.36 is conversion factor to convert cmH2O to mmHg
Pb = 760; % mmHg
Pav = vflow1 * Raw + Pb;
end


% deadspace_expfraction calculates the fraction x of the expired air flow
% rate that comes from the deadspace
% deadvol = the volume of air the dead space contains
% Integral of the (x * [equation for expiratory flow rate]) with bounds 0 
% to texp equals deadvol
% Solves this equation for x, the fraction of the expired air flow rate
% that comes from the deadspace

function x = deadspace_expfrac(TV,RFex,texp)
deadvol = 0.3 * TV;
x = 2*deadvol/TV * (1 / (1 - cos(pi*RFex*texp/30)));
end
 
%calculates the volume compositions of each constituent in the gas in each
%unit
function partial_vol= volume(vper1)
pp=zeros(4,8);
%creates variable to hold the partial pressures for all gas constituents in
%each unit during respiration
sum_vper1=sum(vper1);
vper1_frac=vper1./sum_vper1;
pp(1,:)=[vper1_frac*760];%partial pressures for entry box during inspiration (mmHg)
pp(2,:)=[ppf_O2 ppf_CO2 ppf_N2 ppf_H2O];
%partial pressures for entry box during expiration (mmHg)
pp(3,:)=[ppi_O2 ppi_CO2 ppi_N2 ppi_H2O];
%partial pressures for dead space box during inspiration (mmHg)
pp(4,:)=[ppf_O2 ppf_CO2 ppf_N2 ppf_H2O];
%partial pressures for dead space box during expiration (mmHg)
pp(5,:)=[ppi_O2 ppi_CO2 ppi_N2 ppi_H2O];
%partial pressures for gas exchange box during inspiration (mmHg)
pp(6,:)=[ppf_O2 ppf_CO2 ppf_N2 ppf_H2O];
%partial pressures for gas exchange box during expiration (mmHg)
pp(7,:)=[ppi_O2 ppi_CO2 ppi_N2 ppi_H2O];
%partial pressures for cappilaries box during inspiration (mmHg)
pp(8,:)=[ppf_O2 ppf_CO2 ppf_N2 ppf_H2O];
%partial pressures for cappilaries box during expiration (mmHg)
vol=[29.6 150 65 3000];
%total volume of each unit (mL)
sum_pressures=zeros(8,1);
pp_frac=zeros(8,4);
%finds the partial pressure fractions
for i=1:8
    sum_pressures(i)=sum(pp(i,:));
end
for i=1:8
    for j=1:4
        pp_frac(i,j)=pp(i,j)/sum_pressures(i);
    end
end
partial_vol=volume_calcs(vol,pp_frac);

end

%performs the actual calculations to find partial volume of each gas
%constituent in each unit
function partial_vol= volume_calcs(vol,pp_frac)
for i=1:8
        for j=1:4
             partial_vol(i,j)=pp_frac(i,j)*vol(i);
        end
    
end
end

% travel time calculates the time of air to travel to different locations
% in the lungs
% t_start 4 = the time required for air to travel to the dead space unit
% region, which will be when flow begins in stream 4
% t1 = the time to travel through the extrathoracic region
% t2 = the time required to travel though the trachea
function t_start_4 = traveltime(RFin,TV,index_in,insp_range)
for i = 1:index_in
    int_vflow(i) = -0.5*TV*cos(pi*RFin*insp_range(i)/30)+0.5*TV;
    if abs(int_vflow(i) - 50/1000) < 0.01
        t1 = insp_range(i);
        % the integral of the inspiratory volumetric flow rate from 0 to t1
        % equals the volume of the extrathoracic volume
    end
end
for i = 1:index_in
    if abs(int_vflow(i)*1000/(pi*(1.539/2)^2 - 10.26)) < 0.01
        t2 = insp_range(i);
        % int_vflow(i)*1000/(pi*(1.539/2)^2 - 10.26 is the integral of the 
        % linear velocity
        % this integral from 0 to t2 equals the length of the trachea
    end
end
t_start_4 = t1 + t2; 
end

