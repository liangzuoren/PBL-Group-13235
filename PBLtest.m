
function PBLtest
t_start_6 = entry;
blood(t_start_6)
end

%Entry unit is where inhaled air is humidified and warmed and exhaled air 
%leaves the body from
%physiologically equivalent to the nose, pharynx, larynx, and trachea
function t_start_6 = entry
M = [31.9988 44.0095 28.01348 33.00674];
% M = molar mass of each constituent
TV = 0.5; % liters
H = 1.73; % height in meters of standard man
W = 68; % weight in kilograms of standard man
vper1=[20.95 0.033 78.08 0.937];
%volume percentages of inspired air for O2, CO2, N2, H2O
PP2 = [116.0 32.0 565.0 47.0];
% partial pressures in expired air
PP1 = partial_pressure_calcs(vper1);
% partial pressures in inspired air 50% humidity
PP7 = [100.0 40.0 573.0 47.0];
% partial pressures in alveolar air

vfrac2 = PP2 ./ 760;
vfrac1 = PP1 ./ 760;
vfrac7 = PP7 ./ 760;
% calculates the volume fractions of the constituents from the partial
% pressures

vper_h = humid(TV,vper1,M);
vfrac4 = vper_h./100;
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
tresp = tin + texp;
% sum of inspiration, breathold, and expiration times is the time of one 
% full respiratory cycle

insp_range = 0:0.01:tin;
exp_range = 0:0.01:texp;
resp_range = 0:0.01:tresp;
index_in = length(insp_range);
index_exp = length(exp_range);
index_resp = length(resp_range);

[t_start_4,t_start_6] = traveltime(RFin,TV,index_in,insp_range);
[t_delay_2,t_delay_7] = traveltime(RFex,TV,index_exp,exp_range);
% t_start_4 = time for air to reach dead space during inspiration
% t_start_6 = time for air to reach alveoli during inspiration
% t_delay_2 = time for air to travel from exit of dead space unit 
% (connection of main bronchi to trachea) to the exit of the mouth
% t_delay_7 = time for air to travel from alveolit to exit of mouth during
% expiration
tot_t = t_start_6 + tin + t_delay_7 + texp;
% tot_t = the total time of one respiration cycle
t_start_7 = t_start_6 + tin;
t_start_5 = t_start_6 + tin + t_delay_7 - t_delay_2;
t_start_2 = t_start_6 + tin + t_delay_7;
% finds the times at which air flow begins in streams 2, 5, and 7

range_1a = 0:0.001:tin-0.001;
range_1b = tin:0.001:tot_t;
range_4a = 0:0.001:t_start_4-0.001;
range_4b = t_start_4:0.001:tin+t_start_4-0.001;
range_4c = tin+t_start_4:0.001:tot_t;
range_6a = 0:0.001:t_start_6-0.001;
range_6b = t_start_6:0.001:t_start_7-0.001;
range_6c = t_start_7:0.001:tot_t;
range_7a = 0:0.001:t_start_7-0.001;
range_7b = t_start_7:0.001:t_start_7+texp-0.001;
range_7c = t_start_7+texp:0.001:tot_t;
range_5a = 0:0.001:t_start_5-0.001;
range_5b = t_start_5:0.001:t_start_5+texp-0.001;
range_5c = t_start_5+texp:0.001:tot_t;
range_2a = 0:0.001:t_start_2-0.001;
range_2b = t_start_2:0.001:tot_t;
overall_range = 0:0.001:tot_t;
% sets up ranges for significant sections of time for each stream
% divided into time ranges of zero flow rate and time range of nonzero flow
% rate for each stream

index_1a = length(range_1a);
index_1b = length(range_1b);
index_4a = length(range_4a);
index_4b = length(range_4b);
index_4c = length(range_4c);
index_6a = length(range_6a);
index_6b = length(range_6b);
index_6c = length(range_6c);
index_7a = length(range_7a);
index_7b = length(range_7b);
index_7c = length(range_7c);
index_5a = length(range_5a);
index_5b = length(range_5b);
index_5c = length(range_5c);
index_2a = length(range_2a);
index_2b = length(range_2b);
% indices used in following for loops

for i = 1:index_1a
    vflow1_a(i) = -volumetricflow(RFin,TV,range_1a(i));
    % volumetric flow rate for inspiration in stream 1 in L/s
    % airflow into body defined as negative direction
end
for i = 1:index_1b
    vflow1_b(i) = 0;
    % includes zero flow rates over time interval when there is no flow 
    % rate in stream 1
end
% calculates the flow rates of stream 1 over time intervals 1a and 1b

for i = 1:index_4a
    vflow4_a(i) = 0;
end
for i = 1:index_4b
    vflow4_b(i) = -volumetricflow(RFin,TV,range_4b(i)-t_start_4);
end
for i = 1:index_4c
    vflow4_c(i) = 0;
end
% calculates the flow rates of stream 4 over time intervals 4a,4b,4c

for i = 1:index_6a
    vflow6_a(i) = 0;
end
for i = 1:index_6b
    vflow6_b(i) = - 0.7 * volumetricflow(RFin,TV,range_6b(i)-t_start_6);
end
for i = 1:index_6c
    vflow6_c(i) = 0;
end
% calculates the flow rates of stream 6 over time intervals 6a,6b,6c

for i = 1:index_7a
    vflow7_a(i) = 0;
end
for i = 1:index_7b
    vflow7_b(i) = 0.7 * volumetricflow(RFex,TV,range_7b(i)-t_start_7);
end
for i = 1:index_7c
    vflow7_c(i) = 0;
end
% calculates the flow rates of stream 7 over time intervals 7a,7b,7c

for i = 1:index_5a
    vflow5_a(i) = 0;
end
for i = 1:index_5b
    vflow5_b(i) = volumetricflow(RFex,TV,range_5b(i)-t_start_5);
end
for i = 1:index_5c
    vflow5_c(i) = 0;
end
% % calculates the flow rates of stream 5 over time intervals 5a,5b,5c

for i = 1:index_2a
    vflow2_a(i) = 0;
end
for i = 1:index_2b
    vflow2_b(i) = volumetricflow(RFex,TV,range_2b(i)-t_start_2);
end
% calculates the flow rates of stream 2 over time intervals 2a,2b,2c

vflow1 = [vflow1_a vflow1_b];
vflow4 = [vflow4_a vflow4_b vflow4_c];
vflow6 = [vflow6_a vflow6_b vflow6_c];
vflow7 = [vflow7_a vflow7_b vflow7_c];
vflow5 = [vflow5_a vflow5_b vflow5_c];
vflow2 = [vflow2_a vflow2_b];
% combines functions for entire respiratory cycle

plot(overall_range,vflow1,overall_range,vflow4,overall_range,vflow6,...
    overall_range,vflow7,overall_range,vflow5,overall_range,vflow2)

title('Volumetric Flow Rates of Streams 1, 2, 4, 5, 6, and 7')
xlabel('Time (s)')
ylabel('Volumetric Flow Rate (L/s)')
% legend('1','2','4','5','6','7')
% RV4 = 0;
% composition(vfrac4,overall_range,RFin,TV,RV4,t_start_4)

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

%function partial_pressure_calc finds the partial pressure of inspired air
%at 24 degrees Celsius
%water vapor pressure of water at 24 degrees Celsius is 22.377 mmHg
%PP1=partial pressure values for inspired air (in mmHg)
function PP1= partial_pressure_calcs(vper1)
fracs=vper1./100;
PP1=(760-22.377).*fracs;
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
% vfrac_temp = the volume fraction of the constituents
% M = the molar mass of the constituents
% mratio = the mass ratios of the constituents
% sum_ratio = the sum of the mass ratios
% w = the mass fractions of the constituents

function w = massfrac(vper,M)
vfrac_temp = vper ./ 100;
% calculates volume fractions from volume percentages
mratio = M .* vfrac_temp;
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
% TV = total volume of inspired air
% pp = partial pressure of a constituent in inspired air
% Ta = temperature of inspired air
% R = universal gas constant 
% n = moles of constituent in inspired air
% species_mass = the mass of each species in inspired air
% mass1 = the total mass of inspired air

function mass1 = totalmass(TV,vper,M)
vfrac = vper ./ 100;
v = vfrac .* TV;
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


%PV=nRT
%R=62.36 (L*mmHg)(mol*K)
%{
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
%}
%only for steady state

function blood(timeToAlveoli)
%only for steady state

%vO2inGE2 = 21; %ml/min/mm Hg
%vO2inGE1 = vO2inGE2/0.7*0.3;
%vCO2outGE2 = 451;
%vCO2outGE1 = vCO2outGE2/0.7*0.3;
%vO2outblood = vO2inGE1 + vO2inGE2;
%vCO2inblood = vCO2outGE1 + vCO2outGE2;
H = 1.73; % height in meters of standard man
W = 68; % weight in kilograms of standard man
[RFin,RFex] = RF(H,W);
tin = 30/RFin;
% RFin is number of breaths per minute, divide by 2 to get number of
% inspirations per minute and then multiply by 60s/1min to get time of
% inspiration
texp = 30/RFex;
% RFex is number of breaths per minute, divide by 2 to get number of
% exspirations per minute and then multiply by 60s/1min to get time of
% exspiration
tresp = tin + texp +timeToAlveoli;
% calculates length of time per inspiration and expiration from breathing
% frequencies
PaO2i = 100;%100 mmHg partial pressure in beginning for O2
PaCO2i = 40;%40 mmHg partial pressure in beginning for CO2
Patm = 760; 
PH2O = 47;%presures taken from diagram in assumed values
PN2 = 573;
Ptotal = PaO2i+PaCO2i+PH2O+PN2;
TV = 0.5; %tidal volume = 500 mL O2 perfect at 0.9
DifCapO2 = 21; %mL/min/mmHg
DifCapCO2 = 400; %mL/min/mmHg
xMax = 10;
dO2 =  density(760, 31.9988, 310); %1 atm, 31.9988 g/mol, 310 K(body temperature)
dCO2 = density(760, 44.0095, 310); %1 atm, 44.0095 g/mol, 310 K(body temperature)
trange = 1;
tstep = 0.01;
tsize = length(0:tstep:tresp);
vTotalInitial = 3000;

PaO2alveoliOverTime = zeros(1,tsize*xMax);
PaCO2alveoliOverTime = zeros(1,tsize*xMax);
PH2OalveoliOverTime = zeros(1, tsize*xMax);
PN2alveoliOverTime = zeros(1, tsize*xMax);
PaO2alveoli = PaO2i;
PaCO2alveoli = PaCO2i;
vTotal = vTotalInitial;
PaO2capillary = 40;
PaCO2capillary = 46;

mO2overTime = zeros(1,tsize*xMax);
mCO2overTime = zeros(1,tsize*xMax);
vPercentOverTime = zeros(4, tsize*xMax);
PTotalOverTime = zeros(1, tsize*xMax);


for x=1:xMax %modeling 1 breath
    t=tin+timeToAlveoli;
    mO2inTotal = 0;
    mCO2inTotal = 0;
    Pinspiredair = [158.0 0.3 5.7 596.0];
    %PaO2alveoli = PaO2i;
    %PaCO2alveoli = PaCO2i;
    
    while t<tin+timeToAlveoli+texp
        %mO2out = ((pi*RFex*TV/60)*sin(pi*RFex*t/30)-0.3*mO2inTotal)*dO2/1000*PaO2alveoli/Ptotal/tstep*60;
        %mCO2out = ((pi*RFex*TV/60)*sin(pi*RFex*t/30)-0.3*mCO2inTotal)*dCO2/1000*PaCO2alveoli/Ptotal/tstep*60;
        
        mO2out = ((pi*RFex*TV/60)*sin(pi*RFex*t/30)*0.7)*dO2*PaO2alveoli/Ptotal*tstep;
        mCO2out = ((pi*RFex*TV/60)*sin(pi*RFex*t/30)*0.7)*dCO2*PaCO2alveoli/Ptotal*tstep;
        mH2Oout = ((pi*RFex*TV/60)*sin(pi*RFex*t/30)*0.7)*dCO2*PH2O/Ptotal*tstep;
        mN2out = ((pi*RFex*TV/60)*sin(pi*RFex*t/30)*0.7)*dCO2*PN2/Ptotal*tstep;
        if mO2out > 0
        mO2out = mO2out*-1;
        end
        if mCO2out > 0
        mCO2out = mCO2out*-1;
        end
        if mH2Oout > 0 
        mH2Oout = mH2Oout*-1;
        end
        if mN2out > 0
        mN2out = mN2out*-1;
        end
        mN2out
        
        %difRateO2 = -DifCapO2 * (PaO2alveoli - PaO2capillary)*dO2/1000*tstep/60;
        %difRateCO2 = -DifCapCO2 * (PaCO2alveoli - PaCO2capillary)*dCO2/1000*tstep/60;
        
        if PaO2alveoli < PaO2capillary
        difRateO2 = DifCapO2 * (PaO2capillary - PaO2alveoli) *dO2/1000*tstep/60;
        end
        if PaO2alveoli > PaO2capillary
        difRateO2 = -DifCapO2 * (PaO2alveoli - PaO2capillary) *dO2/1000*tstep/60;
        end
        if PaCO2alveoli < PaCO2capillary
        difRateCO2 = DifCapCO2 * (PaCO2capillary - PaCO2alveoli) *dCO2/1000*tstep/60;
        end
        if PaCO2alveoli > PaCO2capillary
        difRateCO2 = -DifCapCO2 * (PaCO2alveoli - PaCO2capillary) *dCO2/1000*tstep/60;
        end
        
        
        mO2 = difRateO2 + mO2out;
        mCO2 = difRateCO2 + mCO2out;
        %{
        mO2alveoli = mO2alveoli + mO2;
        
        mCO2alveoli = mCO2alveoli - mCO2;
        
        nO2 = mO2alveoli/31.9988;
        nCO2 = mCO2alveoli/44.0095;
        
        PaO2alveoli = nO2*tstep*62.36367*1000*310/3;
        PaCO2alveoli = nCO2*tstep*62.36367*1000*310/3;
        %}
        nO2 = mO2/31.9988;
        nCO2 = mCO2/44.0095;
        nN2 = mN2out/28.014;
        nH2O = mH2Oout/18.01528;
        %PV=nRT, P = nRT/V
        
        PaO2diff = nO2*62.36367*1000*310/vTotal; %gives pressure in mmHg
        PaCO2diff = nCO2*62.36367*1000*310/vTotal;
        PN2diff = nN2*62.36367*1000*310/vTotal;
        PH2Odiff = nH2O*62.36367*1000*310/vTotal;
        
        nTotal = PaO2alveoli*vTotal/310/62.36367/1000;
        
        PaO2alveoli = PaO2alveoli + PaO2diff;
        PaCO2alveoli = PaCO2alveoli + PaCO2diff;
        PN2 = PN2 + PN2diff;
        PH2O = PH2O + PH2Odiff;
        
        vTotal = nTotal*62.36367*1000*310/PaO2alveoli;
        
        Ptotal = PaO2alveoli + PaCO2alveoli + PH2O + PN2;
        PTotalOverTime(trange) = Ptotal;
        PaO2alveoliOverTime(trange)= PaO2alveoli;
        PaCO2alveoliOverTime(trange)= PaCO2alveoli;
        PN2alveoliOverTime(trange) = PN2;
        PH2OalveoliOverTime(trange) = PH2O;
        mO2overTime(trange) = mO2out;
        mCO2overTime(trange) = mCO2out;
        t=t+tstep;
        trange = trange+1;
    end 
    
    while t<texp + timeToAlveoli*2 + tin
        if PaO2alveoli < PaO2capillary
        difRateO2 = DifCapO2 * (PaO2capillary - PaO2alveoli) *dO2/1000*tstep/60;
        end
        if PaO2alveoli > PaO2capillary
        difRateO2 = -DifCapO2 * (PaO2alveoli - PaO2capillary) *dO2/1000*tstep/60;
        end
        if PaCO2alveoli < PaCO2capillary
        difRateCO2 = DifCapCO2 * (PaCO2capillary - PaCO2alveoli) *dCO2/1000*tstep/60;
        end
        if PaCO2alveoli > PaCO2capillary
        difRateCO2 = -DifCapCO2 * (PaCO2alveoli - PaCO2capillary) *dCO2/1000*tstep/60;
        end
        mO2 = difRateO2;
        mCO2 = difRateCO2;
        %{
        mO2alveoli = mO2alveoli + mO2;
        mCO2alveoli = mCO2alveoli - mCO2;
        
        nO2 = mO2alveoli/31.9988;
        nCO2 = mCO2alveoli/44.0095;
        
        PaO2alveoli = nO2*tstep*62.36367*1000*310/3;
        PaCO2alveoli = nCO2*tstep*62.36367*1000*310/3;
        %}
        nO2 = mO2/31.9988;
        nCO2 = mCO2/44.0095;
        %PV=nRT, P = nRT/V
        
        PaO2diff = nO2*62.36367*1000*310/vTotal; %gives pressure in mmHg
        PaCO2diff = nCO2*62.36367*1000*310/vTotal;
        
        nTotal = PaO2alveoli*vTotal/310/62.36367/1000;
        
        PaO2alveoli = PaO2alveoli + PaO2diff;
        PaCO2alveoli = PaCO2alveoli + PaCO2diff;
        
        vTotal = nTotal*62.36367*1000*310/PaO2alveoli;
        
        Ptotal = PaO2alveoli + PaCO2alveoli + PH2O + PN2;
        PTotalOverTime(trange) = Ptotal;
        PaO2alveoliOverTime(trange)= PaO2alveoli;
        PaCO2alveoliOverTime(trange)= PaCO2alveoli;
        PN2alveoliOverTime(trange) = PN2;
        PH2OalveoliOverTime(trange) = PH2O;
        
        mO2overTime(trange) = 0;
        mCO2overTime(trange) = 0;
        
        t=t+tstep;
        trange = trange+1;
    end 
    
    while t<(tin*2+texp+timeToAlveoli*2)
        
        mO2in = -0.7*(pi*RFin*TV/60)*sin(pi*RFin*t/30)*0.218*dO2*tstep; %insert partial pressures/total pressure or mass frac for humidified air here
        mCO2in = -0.7*(pi*RFin*TV/60)*sin(pi*RFin*t/30)*0.0003*dCO2*tstep;
        mH2Oin = -0.7*((pi*RFin*TV/60)*sin(pi*RFin*t/30))*dCO2*0.0697*tstep;
        mN2in = -0.7*((pi*RFin*TV/60)*sin(pi*RFin*t/30))*dCO2*0.78*tstep;
        
        if mO2in < 0
        mO2in = mO2in*-1;
        end
        if mCO2in < 0
        mCO2in = mCO2in*-1;
        end
        if mH2Oin < 0 
        mH2Oin = mH2Oin*-1;
        end
        if mN2in < 0
        mN2in = mN2in*-1;
        end
        
        
        if PaO2alveoli < PaO2capillary
        difRateO2 = DifCapO2 * (PaO2capillary - PaO2alveoli) *dO2/1000*tstep/60;
        end
        if PaO2alveoli > PaO2capillary
        difRateO2 = -DifCapO2 * (PaO2alveoli - PaO2capillary) *dO2/1000*tstep/60;
        end
        if PaCO2alveoli < PaCO2capillary
        difRateCO2 = DifCapCO2 * (PaCO2capillary - PaCO2alveoli) *dCO2/1000*tstep/60;
        end
        if PaCO2alveoli > PaCO2capillary
        difRateCO2 = -DifCapCO2 * (PaCO2alveoli - PaCO2capillary) *dCO2/1000*tstep/60;
        end
        mO2 = difRateO2 + mO2in;
        mCO2 = difRateCO2 + mCO2in;
        %{
        mO2alveoli = mO2alveoli + mO2;
        mCO2alveoli = mCO2alveoli - mCO2;
        
        nO2 = mO2alveoli/31.9988;
        nCO2 = mCO2alveoli/44.0095;
        
        PaO2alveoli = nO2*tstep*62.36367*1000*310/3;
        PaCO2alveoli = nCO2*tstep*62.36367*1000*310/3;
        %}
        
        nO2 = mO2/31.9988;
        nCO2 = mCO2/44.0095;
        nN2 = mN2in/28.014;
        nH2O = mH2Oin/18.01528;
        %PV=nRT, P = nRT/V
        
        PaO2diff = nO2*62.36367*1000*310/vTotal; %gives pressure in mmHg
        PaCO2diff = nCO2*62.36367*1000*310/vTotal;
        PN2diff = nN2*62.36367*1000*310/vTotal;
        PH2Odiff = nH2O*62.36367*1000*310/vTotal;
        
        nTotal = PaO2alveoli*vTotal/310/62.36367/1000;
        
        PaO2alveoli = PaO2alveoli + PaO2diff;
        PaCO2alveoli = PaCO2alveoli + PaCO2diff;
        PN2 = PN2 + PN2diff;
        PH2O = PH2O + PH2Odiff;
       
        vTotal = nTotal*62.36367*1000*310/PaO2alveoli;
       
        Ptotal = PaO2alveoli + PaCO2alveoli + PH2O + PN2;
        PTotalOverTime(trange) = Ptotal;
        mO2overTime(trange) = mO2in;
        mCO2overTime(trange) = mCO2in;
        PaO2alveoliOverTime(trange)= PaO2alveoli;
        PaCO2alveoliOverTime(trange)= PaCO2alveoli;
        PN2alveoliOverTime(trange) = PN2;
        PH2OalveoliOverTime(trange) = PH2O;
        
        %mTotal = residualVolume + vflow6_in(trange)*dInspiredAir; %volumetric flow rate * density?
        %PaO2alveoli = (PaO2alveoli * mTotal- mO2 )/ mTotal;%PaO2 = FiO2*(Patm - PH2O) - (PaCO2/RQ)
        %PaCO2alveoli = (PaCO2alveoli * mTotal - mCO2 )/ mTotal;
        mO2inTotal = mO2inTotal + mO2in;
        mCO2inTotal = mCO2inTotal + mCO2in;
        t=t+tstep;
        trange = trange+1;
    end
end    


resp_range = 0:0.01:length(PaO2alveoliOverTime)*tstep-tstep;
length(resp_range);
vPercentOverTime(1,:)= (PaO2alveoliOverTime./PTotalOverTime*100); 
vPercentOverTime(2,:) = (PaCO2alveoliOverTime./PTotalOverTime*100);
vPercentOverTime(3,:) = (PH2OalveoliOverTime./PTotalOverTime*100);
vPercentOverTime(4,:) = (PN2alveoliOverTime./PTotalOverTime*100);
figure
plot(resp_range, PaO2alveoliOverTime);
title('Partial Pressure of O2 in Alveoli Over Time')
xlabel('Time (s)')
ylabel('Partial Pressure (mmHg)')
figure
plot(resp_range, PaCO2alveoliOverTime)
title('Partial Pressure of CO2 in Alveoli Over Time')
xlabel('Time (s)')
ylabel('Partial Pressure (mmHg)')
figure
plot(resp_range, mO2overTime);
title('Mass Flow Rate of O2 In and Out of Alveoli Over Time')
xlabel('Time (s)')
ylabel('Mass Flow Rate (g/s)')
figure
plot(resp_range,mCO2overTime);
title('Mass Flow Rate of CO2 In and Out of Alveoli Over Time')
xlabel('Time (s)')
ylabel('Mass Flow Rate (g/s)')
figure
plot(resp_range, vPercentOverTime(1,:));
title('Volume Percentage of O2 in Alveoli Over Time')
xlabel('Time (s)')
ylabel('Volume Percentage')
figure
plot(resp_range, vPercentOverTime(2,:));
title('Volume Percentage of CO2 in Alveoli Over Time')
xlabel('Time (s)')
ylabel('Volume Percentage')
figure
plot(resp_range, vPercentOverTime(3,:));
title('Volume Percentage of H2O in Alveoli Over Time')
xlabel('Time (s)')
ylabel('Volume Percentage')
figure
plot(resp_range, vPercentOverTime(4,:));
title('Volume Percentage of N2 in Alveoli Over Time')
xlabel('Time (s)')
ylabel('Volume Percentage')

end

%calcuate mass flow rateswith regards to time
%model diffusion over time
%partial pressure change over time
%
%
%humid calculates volume percentage in humidfied air

%humidity function calculates the volume percentages of the air constiuents


%humid function calculates the volume percentages of the air constiuents
%after is has been humidified from 50% (when inspired) to 100%

%m_i=initial amount of water in air when inspired (g)
%m_f=final amount of water in air after humdidification (g)
%water_added=total amount of water added to air when it is humidified (g)
%mass_tot1=total mass of gas before humidified (g)
%mass_tot2= total mass of gas after humidified (g)
%w_i= initial mass fractions of gas constiuents before humidified
%mass_i=amount of each consituent in inhaled air (g)
%mass_f=amount of each consituent in humidified air (g)
%w_f= mass fractions of gas constituents after humidified
%vper_1=volume percentages of air in stream before it is humidified
%vper_h= volume percentages of air after it is humidified
%VT= volume of air inhaled- 0.5 L in this model

function vper_h=humid(vtot,vper1,M)
m_i=.0045; 
%at 24 degrees Celsius, amount of water in air before humidified (g), at
%50% humidity, for inspired air
m_f=.022;
%at 37 degrees Celsius, amount of water in air after humidified (g), at 100% humidity
water_added=m_f-m_i;
density=1.184; 
%density of inhaled air (g/L)
mass_tot1=density*vtot;
%total mass of inhaled air (g)
mass_tot2=mass_tot1+water_added;
%total mass of humidifed air (g)
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

%volumepercentage function calculates the volume percentages of 
%constituents in a stream based on the mass fractions of each constituent
%mass_f=mass of all constitiuents (mg)
%M= molar masses of all constituents (g/mol)
%v_fracf=volumetric fractions of all gas constituents after humidified
%vper_f= volume percentages of all gas constituents after humidified
function vper_f=volumepercentage(mass_f,M)
nratio=(mass_f)./M;  
%calculates number of moles of all gas constituents after humidified
nratio_sum= sum(nratio);  
%calculates sum of all moles of gas after humidified
vfrac_f= nratio./nratio_sum;  
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
 
%calculates the volume compositions of each constituent in the gas in each
%unit
function partial_vol= volume(PP1)
pp=zeros(4,8);
%creates variable to hold the partial pressures for all gas constituents in
%each unit during respiration
pp(1,:)=[PP1];%partial pressures for entry box during inspiration (mmHg)
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
vol=[29.6 29.6 150 150 65 65 3000 3000];
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

% traveltime calculates the time of air to travel to different locations
% in the lungs
% t1 = the time to travel through the extrathoracic region
% t2 = the time required to travel though the trachea
% t_start_4 = the time required for air to first reach the dead space unit,
% (the time required to reach the generation 1 bronchi), which will be when 
% flow begins in stream 4
% t_start_6 = the time required for air to first reach the alveoli
% (generation 17) 
function [t_start_4,t_start_6] = traveltime(RF,TV,index,range)
d = [1.539 1.043 0.71 0.479 0.385 0.299 0.239 0.197 0.159 0.132 0.111 ...
    0.093 0.081 0.07 0.063 0.056 0.051 0.046 0.043 0.04 0.038 0.037 ...
    0.035 0.035];
% d = the diameters of each generations of tubes in cm, based on scaled 
% Weibel A model
r = d ./ 2;
% r = the radii of each generation of tubes in cm
A = pi*r.^2;
% A = the cross-sectional area of each generation of tubes in cm^2
l = [10.26 4.07 1.624 0.65 1.086 0.915 0.769 0.65 0.547 0.462 0.393 ...
    0.333 0.282 0.231 0.197 0.171 0.141 0.121 0.1 0.085 0.071 0.06 ...
    0.05 0.043];
% l = the length of each generation of tubes in centimeters
% the trachea is generation 0, which corresponds to d(1), r(1), and l(1)
for i = 1:index
    int_vflow(i) = antiderivative(RF,TV,range(i))- antiderivative(RF,TV,0);
    if abs(int_vflow(i) - 50/1000) < 0.01
        t1 = range(i);
        % the integral of the inspiratory volumetric flow rate from 0 to t1
        % equals the volume of the extrathoracic volume
    end
end
for i = 1:index
    x = int_vflow(i)*1000/A(1);
    if abs(int_vflow(i)*1000/A(1) - l(1)) < 0.5
        t2 = range(i);
        % int_vflow(i)*1000/(pi*(r(1))^2 is the integral of the 
        % linear velocity
        % this integral from 0 to t2 equals the length of the trachea l(1)
    end
end

t_start_4 = t1 + t2; 
% the time to first reach the dead space unit (time to reach generation 1 
% brochi) is the sum of the time to travel through the extrathoracic region 
% and the time to travel through the trachea

for i = 2:17
    for j =1:index
        if abs(int_vflow(j)*1000/A(i) - l(i)) < 0.1
            t_dead_gen(i) = range(j);
        end
    end
end

t_dead_sum = sum(t_dead_gen);
% the first alveoli appear at generation 17
% the time to travel from the exit of the trachea to the first alveoli is 
% the sum of the time to travel though the generation 1 through 16 airways

t_start_6 = t_start_4 + t_dead_sum;
% t_start_6 = the total time required for air to reach the first alveoli
end

%This function finds the volume compositions of gas constituents in the
%streams where volume composition does not change over time-Stream 1, 4, 6
%vflow= contains all the flow rates of the diagram
%vfrac= contains all the volume fractions
function partial_vol= volume(vflow, vfrac)
partial_vol=zeros(6,length(vflow));
for i=1:length(vflow)
    partial_vol=vfrac.*vflow;
end
end

% RV = residual volume that remains in unit
% only aplies for stream 4 right now
function composition(vfrac,overall_range,RF,TV,RV,t_start_4)
V = zeros(4,length(overall_range));
V(:,1:t_start_4) = RV .* vfrac';
for i = 1:length(vfrac)
    for j = t_start_4:length(overall_range)
        V(i,j) = V(i,j-1) + antiderivative(RF,TV,j)-antiderivative(RF,TV,j-1);
        for k = 1:length(overall_range)
            Vs = V(i,:);
            Vtot(j) = sum(Vs);
        end
        vfrac(i,j) = V(i,j)/Vtot(j);
    end
end
plot(overall_range,vfrac(1,:))
title('Partial Pressure of Oxygen in Stream 4')
xlabel('Time (s)')
ylabel('PO2 (mmHg)')
end

function output = antiderivative(RF,TV,t)
output = -0.5*TV*cos(pi*RF*t/30);
end




