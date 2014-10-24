% The entry unit is where inhaled air is humidified and warmed and exhaled 
% air leaves the body from
% The entry unit is physiologically equivalent to the nose, pharynx, 
% larynx, and trachea
% The dead space unit is the area where air travels through the unit but no 
% gas exchange, heat transfer, or humidification occurs
% The dead space unit is physiologically equivalent to the generation 1 
% through 16 bronci
% The gas exchange unit is where gas exchange occurs
% The gas exchnage unit is physiologically equivalent to the alveoli,
% including the alveoli of respiratory bronchioles, alveolar ducts, and
% alveolar sacs

function PBLtest2
TV = 0.5; % liters
H = 1.73; % height in meters of standard man
W = 68; % weight in kilograms of standard man

[RFin,RFex,tin,texp,t_start_4,t_start_6,tot_t,t_start_7,t_start_5,...
    t_start_2,t_delay_7,index_in,insp_range] = time(TV,H,W);

[vPercentOverTime] = blood(t_start_6,tin,t_delay_7,texp,tot_t);

airflow(RFin,RFex,tin,texp,t_start_4,t_start_6,tot_t,t_start_7,t_start_5,...
    t_start_2,TV,index_in,insp_range,vPercentOverTime)
end

% time calculates time intervals of significance to the respiratory cycle
% time also calculates the volume of the entry unit and the dead sapce unit
% Volume_entry= the volume of the entry unit in mL
% Volume_deadspace = the volume of the dead space unit in mL

function [RFin,RFex,tin,texp,t_start_4,t_start_6,tot_t,t_start_7,...
    t_start_5,t_start_2,t_delay_7,index_in,insp_range] = time(TV,H,W)

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

insp_range = 0:0.001:tin;
exp_range = 0:0.001:texp;
index_in = length(insp_range);
index_exp = length(exp_range);

[t_start_4,t_start_6,~,~] = traveltime(RFin,TV,index_in,insp_range);
[t_delay_2,t_delay_7,Volume_deadspace,Volume_entry] = traveltime...
    (RFex,TV,index_exp,exp_range);
% t_start_4 = time for air to reach dead space during inspiration
% t_start_6 = time for air to reach alveoli during inspiration
% t_delay_2 = time for air to travel from exit of dead space unit 
% (connection of main bronchi to trachea) to the exit of the mouth
% t_delay_7 = time for air to travel from alveoli to exit of mouth during
% expiration
tot_t = t_start_6 + tin + t_delay_7 + texp;
% tot_t = the total time of one respiration cycle
t_start_7 = t_start_6 + tin;
t_start_5 = t_start_6 + tin + t_delay_7 - t_delay_2;
t_start_2 = t_start_6 + tin + t_delay_7;
% finds the times at which air flow begins in streams 2, 5, and 7

Volume_deadspace = Volume_deadspace
Volume_entry = Volume_entry

end

% The airflow function calculates: 
% The overall volumetric flow rate of streams 1, 2, 3, 4, 5, 6, and 7 over
% one full respiratory cycle
% The volumetric flow rate of each constituent in each of these streams 
% over one full respiratory cycle 
% The partial pressures of each constituent in the entry and dead space 
% units over one full respiratory cycle
% The volume percentages of each constituent in the entry and dead space 
% units over one full respiratory cycle 
% Returns graphs for each of these calulations

function airflow(RFin,RFex,tin,texp,t_start_4,t_start_6,tot_t,t_start_7,t_start_5,...
    t_start_2,TV,index_in,insp_range,vPercentOverTime)


vfrac_ex_O2 = vPercentOverTime(1,:) ./ 100;
vfrac_ex_CO2 = vPercentOverTime(2,:) ./ 100;
vfrac_ex_H2O = vPercentOverTime(3,:) ./ 100;
vfrac_ex_N2 = vPercentOverTime(4,:) ./ 100;

M = [31.9988 44.0095 28.01348 33.00674];
% M = molar mass of each constituent
M_air = sum(M);
R = 62.3637;
Tb = 310; % body temperature
air_density=1.184; % g/L

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
% calculates composition after humidification in Entry unit
vfrac6 = vfrac4;
% the compositions of stream 6 and 4 are equal because stream 4 is
% the only inlet stream and stream 6 is the only outlet stream of the dead 
% space unit, and there are no reactions
vfrac3 = [0 0 0 1];

range_1a = 0:0.001:tin;
range_1b = tin+0.001:0.001:tot_t;
range_4a = 0:0.001:t_start_4-0.001;
range_4b = t_start_4:0.001:tin+t_start_4;
range_4c = tin+t_start_4+0.001:0.001:tot_t;
range_6a = 0:0.001:t_start_6-0.001;
range_6b = t_start_6:0.001:t_start_7;
range_6c = t_start_7+0.001:0.001:tot_t;
range_7a = 0:0.001:t_start_7-0.001;
range_7b = t_start_7:0.001:t_start_7+texp;
range_7c = t_start_7+texp+0.001:0.001:tot_t;
range_5a = 0:0.001:t_start_5-0.001;
range_5b = t_start_5:0.001:t_start_5+texp;
range_5c = t_start_5+texp+0.001:0.001:tot_t;
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
index_overall = length(overall_range);
% indices used in following for loops

for i = 1:index_1a
    vflow1_a(i) = volumetricflow(RFin,TV,range_1a(i));
    % volumetric flow rate for inspiration in stream 1 in L/s
    % airflow into body defined as positive direction
    intflow1_a(i) = - antiderivative(RFin,TV,range_1a(i));
end
for i = 1:index_1b
    vflow1_b(i) = 0;
    % includes zero flow rates over time interval when there is no flow 
    % rate in stream 1
    intflow1_b(i) = - antiderivative(RFin,TV,range_1a(index_1a));
end

% calculates the flow rates of stream 1 over time intervals 1a and 1b

for i = 1:index_4a
    vflow4_a(i) = 0;
    intflow4_a(i) = - antiderivative(RFin,TV,range_4b(1)-t_start_4);
end
for i = 1:index_4b
    vflow4_b(i) = volumetricflow(RFin,TV,range_4b(i)-t_start_4);
    intflow4_b(i) = - antiderivative(RFin,TV,range_4b(i)-t_start_4);
end
for i = 1:index_4c
    vflow4_c(i) = 0;
    intflow4_c(i) = - antiderivative(RFin,TV,range_4b(index_4b)-t_start_4);
end
% calculates the flow rates of stream 4 over time intervals 4a,4b,4c

for i = 1:index_6a
    vflow6_a(i) = 0;
    intflow6_a(i) = 0.7 * - antiderivative(RFin,TV,range_6b(1)-t_start_6);
end
for i = 1:index_6b
    vflow6_b(i) = 0.7 * volumetricflow(RFin,TV,range_6b(i)-t_start_6);
    intflow6_b(i) = - 0.7 * antiderivative(RFin,TV,range_6b(i)-t_start_6);
end
for i = 1:index_6c
    vflow6_c(i) = 0;
    intflow6_c(i) = - 0.7 * antiderivative(RFin,TV,range_6b(index_6b)-t_start_6);
end
% calculates the flow rates of stream 6 over time intervals 6a,6b,6c

for i = 1:index_7a
    vflow7_a(i) = 0;
    intflow7_a(i) = 0.7 * antiderivative(RFex,TV,range_7b(1)-t_start_7);
end
for i = 1:index_7b
    vflow7_b(i) = - 0.7 * volumetricflow(RFex,TV,range_7b(i)-t_start_7);
    intflow7_b(i) = 0.7 * antiderivative(RFex,TV,range_7b(i)-t_start_7);
end
for i = 1:index_7c
    vflow7_c(i) = 0;
    intflow7_c(i) = 0.7 * antiderivative(RFex,TV,range_7b(index_7b)-t_start_7);
end
% calculates the flow rates of stream 7 over time intervals 7a,7b,7c

for i = 1:index_5a
    vflow5_a(i) = 0;
    intflow5_a(i) = antiderivative(RFex,TV,range_5b(1)-t_start_5);
end
for i = 1:index_5b
    vflow5_b(i) = - volumetricflow(RFex,TV,range_5b(i)-t_start_5);
    intflow5_b(i) = antiderivative(RFex,TV,range_5b(i)-t_start_5);
end
for i = 1:index_5c
    vflow5_c(i) = 0;
    intflow5_c(i) = antiderivative(RFex,TV,range_5b(index_5b)-t_start_5);
end
% calculates the flow rates of stream 5 over time intervals 5a,5b,5c

for i = 1:index_2a
    vflow2_a(i) = 0;
    intflow2_a(i) = antiderivative(RFex,TV,range_2b(1)-t_start_2);
end
for i = 1:index_2b
    vflow2_b(i) = - volumetricflow(RFex,TV,range_2b(i)-t_start_2);
    intflow2_b(i) = antiderivative(RFex,TV,range_2b(i)-t_start_2);
end
% calculates the flow rates of stream 2 over time intervals 2a,2b,2c

vflow1 = [vflow1_a vflow1_b];
vflow4 = [vflow4_a vflow4_b vflow4_c];
vflow6 = [vflow6_a vflow6_b vflow6_c];
vflow7 = [vflow7_a vflow7_b vflow7_c];
vflow5 = [vflow5_a vflow5_b vflow5_c];
vflow2 = [vflow2_a vflow2_b];
% combines the sections of the volumetric flow rate function for the entire 
% respiratory cycle

intflow1 = [intflow1_a intflow1_b];
intflow4 = [intflow4_a intflow4_b intflow4_c];
intflow6 = [intflow6_a intflow6_b intflow6_c];
intflow7 = [intflow7_a intflow7_b intflow7_c];
intflow5 = [intflow5_a intflow5_b intflow5_c];
intflow2 = [intflow2_a intflow2_b];
% combines the sections of the integrals of the volumetric flow rate 
% function for the entire respiratory cycle

for i = 1:length(vfrac1)
    for j = 1:index_overall
        vflows1(i,j) = vfrac1(i)' * vflow1(j);
    end
end
% calculates the volumetric flow rates of constituents in stream 1

for i = 1:length(vfrac4)
    for j = 1:index_overall
        vflows4(i,j) = vfrac4(i)' * vflow4(j);
    end
end
% calculates the volumetric flow rates of constituents in stream 4

for i = 1:length(vfrac6)
    for j = 1:index_overall
        vflows6(i,j) = vfrac6(i)' * vflow6(j);
    end
end
% calculates the volumetric flow rates of constituents in stream 6

for j = 1:index_overall
    vflow2O2(j) = vfrac_ex_O2(j) * vflow2(j);
    vflow2CO2(j) = vfrac_ex_CO2(j) * vflow2(j);
    vflow2N2(j) = vfrac_ex_N2(j) * vflow2(j);
    vflow2H2O(j) = vfrac_ex_H2O(j) * vflow2(j);
end

% calculates the volumetric flow rates of constituents in stream 2 


for j = 1:index_overall
    vflow5O2(j) = vfrac_ex_O2(j) * vflow5(j);
    vflow5CO2(j) = vfrac_ex_CO2(j) * vflow5(j);
    vflow5N2(j) = vfrac_ex_N2(j) * vflow5(j);
    vflow5H2O(j) = vfrac_ex_H2O(j) * vflow5(j);
end

% calculates the volumetric flow rates of constituents in stream 5

for j = 1:index_overall
    vflow7O2(j) = vfrac_ex_O2(j) * vflow7(j);
    vflow7CO2(j) = vfrac_ex_CO2(j) * vflow7(j);
    vflow7N2(j) = vfrac_ex_N2(j) * vflow7(j);
    vflow7H2O(j) = vfrac_ex_H2O(j) * vflow7(j);
end
% calculates the volumetric flow rates of constituents in stream 7

figure
plot(overall_range,vflows1(1,:),overall_range,vflows4(1,:),overall_range,vflows6(1,:),overall_range,vflow2O2,overall_range,vflow5O2,overall_range,vflow7O2)
title('Volumetric Flow Rates of Oxygen Through Airways')

figure
plot(overall_range,vflows1(2,:),overall_range,vflows4(2,:),overall_range,vflows6(2,:),overall_range,vflow2CO2,overall_range,vflow5CO2,overall_range,vflow7CO2)
title('Volumetric Flow Rates of Carbon Dioxide Through Airways')

figure
plot(overall_range,vflows1(3,:),overall_range,vflows4(3,:),overall_range,vflows6(3,:),overall_range,vflow2N2,overall_range,vflow5N2,overall_range,vflow7N2)
title('Volumetric Flow Rates of Nitrogen Through Airways')

figure
plot(overall_range,vflows1(4,:),overall_range,vflows4(4,:),overall_range,vflows6(4,:),overall_range,vflow2H2O,overall_range,vflow5H2O,overall_range,vflow7H2O)
title('Volumetric Flow Rates of Water Vapor Through Airways')

RVentry = 0;
[VO2entry,VCO2entry,VN2entry,VH2Oentry,Vtotentry, ...
    vperO2entry,vperCO2entry,vperN2entry,...
    vperH2Oentry] = composition(vPercentOverTime,vfrac1,intflow1,intflow4,intflow5,...
    intflow2,RVentry,index_overall,M_air);
% calculates the volumes and partial pressures of constituents in the entry
% box over one full respiratory cycle

RVds = 0;
[VO2ds,VCO2ds,VN2ds,VH2Ods,Vtotds,PO2ds,PCO2ds,PN2ds,PH2Ods,vperO2ds,...
    vperCO2ds,vperN2ds,vperH2Ods] = composition(vPercentOverTime,vfrac4,intflow4,...
    intflow6,intflow7,intflow5,RVds,length(overall_range),M_air);

 
nflownet_entry = molrate(vflow1,-vflow4,zeros(length(vflow1)),vfrac1,vfrac4,vfrac3,air_density,M_air,index_overall);
% doesn't take water flowing in stream 3 into account
nflowO2entry = nflownet_entry(1,:);
nflowCO2entry = nflownet_entry(2,:);
nflowN2entry = nflownet_entry(3,:);
nflowH2Oentry = nflownet_entry(4,:);

PO2entry = zeros(1,length(nflowO2entry));
PCO2entry = zeros(1,length(nflowCO2entry));
PN2entry = zeros(1,length(nflowN2entry));
PH2Oentry = zeros(1,length(nflowH2Oentry));


for i=2:length(nflowO2entry)
   if VO2entry(i) > 0;
       PO2diff = (nflowO2entry(i)*R*Tb*0.001)/(VO2entry(i)*10);
       PO2entry(i) = PO2entry(i-1) + PO2diff; 
   else
       PO2entry(i) = PO2entry(i-1);
   end    
       
   if VCO2entry(i)>0;
   PCO2diff = (nflowCO2entry(i)*R*Tb*0.001)/(VCO2entry(i)*10);
   PCO2entry(i) = PCO2entry(i-1) + PCO2diff;
   else
   PCO2entry(i) = PCO2entry(i-1);
   end
   
   if VN2entry(i)>0;
   PN2diff = (nflowN2entry(i)*R*Tb*0.001)/(VN2entry(i)*10);
   PN2entry(i) = PN2entry(i-1) + PN2diff;
   else
   PN2entry(i) =  PN2entry(i-1);
   end
   
   if VH2Oentry(i)>0;
   PH2Odiff = (nflowH2Oentry(i)*R*Tb*0.001)/(VH2Oentry(i)*10);
   PH2Oentry(i) = PH2Oentry(i-1) + PH2Odiff; 
   else
   PH2Oentry(i) = PH2Oentry(i-1);
   end
end
   
%PO2entry =(nflowO2entry*R*Tbss*0.001)./VO2entry;
%PCO2entry =(nflowCO2entry*R*Tb*0.001)./VCO2entry;
%PN2entry =(nflowNs2entry*R*Tb*0.001)./VN2entry;
%PH2Oentry =(nflowH2Oentry*R*Tb*0.001)./VH2Oentry;


figure
plot(overall_range,PO2entry)
title('Pressure of Oxygen in Entry Unit')

figure
plot(overall_range,PCO2entry)
title('Pressure of Carbon Dioxide in Entry Unit')

figure
plot(overall_range,PN2entry)
title('Pressure of Nitrogen in Entry Unit')

figure
plot(overall_range,PH2Oentry)
title('Pressure of Water in Entry Unit')

figure
plot(overall_range,vflow1,overall_range,vflow4,overall_range,vflow6,...
    overall_range,vflow7,overall_range,vflow5,overall_range,vflow2)

title('Volumetric Flow Rates of Streams 1, 2, 4, 5, 6, and 7')
xlabel('Time (s)')
ylabel('Volumetric Flow Rate (L/s)')
% legend('1','2','4','5','6','7')

figure
plot(overall_range,vflows1(1,:),overall_range,vflows1(2,:),...
    overall_range,vflows1(3,:),overall_range,vflows1(4,:),overall_range,vflow1)
title('Volumetric Flow Rate of Constituents in Stream 1')
xlabel('Time (s)')
ylabel('Volumetric Flow Rate (L/s)')

figure
plot(overall_range,vflows4(1,:),overall_range,vflows4(2,:),...
    overall_range,vflows4(3,:),overall_range,vflows4(4,:),overall_range,vflow4)
title('Volumetric Flow Rate of Constituents in Stream 4')
xlabel('Time (s)')
ylabel('Volumetric Flow Rate (L/s)')

figure
plot(overall_range,vflows6(1,:),overall_range,vflows6(2,:),...
    overall_range,vflows6(3,:),overall_range,vflows6(4,:),overall_range,vflow6)
title('Volumetric Flow Rate of Constituents in Stream 6')
xlabel('Time (s)')
ylabel('Volumetric Flow Rate (L/s)')




% volcont1 = zeros(1,index_overall);
% volcont_temp1(1) = 0;
% volcont2 = zeros(1,index_overall);
% volcont_temp2(1) = 0;
% % volcont4 = zeros(1,index_overall);
% % volcont_temp4(1) = 0;
% for i = 2:index_overall
%     volcont1(i) = volcont_temp1(i-1) + intflow1(i) - intflow1(i-1);
%     volcont_temp1(i) = volcont1(i);
% %     volcont4(i) = volcont_temp4(i-1) + intflow4(i) - intflow4(i-1);
% %     volcont_temp4(i) = volcont4(i);
%     volcont2(i) = volcont_temp2(i-1) + intflow2(i) - intflow2(i-1);
%     volcont_temp2(i) = volcont2(i);
% end



figure
plot(overall_range,PO2entry)
title('Partial Pressure of Oxygen in Entry Unit')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
% plots the partial pressures of constituents in the entry box over time
figure
plot(overall_range,VO2entry,overall_range,VCO2entry,overall_range,...
    VN2entry,overall_range,VH2Oentry,overall_range,Vtotentry)
title('Volume of Constituents in Entry Unit')
xlabel('Time (s)')
ylabel(' Volume (L)')
% plots the volumes of constituents in the entry box over time
figure
plot(overall_range,vperO2entry,overall_range,vperCO2entry,overall_range,...
    vperN2entry,overall_range,vperH2Oentry)
title('Volume Percentages of Constituents in Entry Unit')
xlabel('Time (s)')
ylabel(' Volume Percent (%)')

% NEED TO ADD HUMIDIFICATION OCCURING OVER TIME IN ENTRY



figure
plot(overall_range,PO2ds,overall_range,PCO2ds,overall_range,...
    PN2ds,overall_range,PH2Ods)
title('Partial Pressure of Constituents in Dead Space Unit')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
% plots the partial pressures of constituents in the dead space unit over 
% time
figure
plot(overall_range,VO2ds,overall_range,VCO2ds,overall_range,...
    VN2ds,overall_range,VH2Ods,overall_range,Vtotds)
title('Volume of Constituents in Dead Space Unit')
xlabel('Time (s)')
ylabel(' Volume (L)')
% plots the volumes of constituents in the dead space unit over time

figure
plot(overall_range,vperO2ds,overall_range,vperCO2ds,overall_range,...
    vperN2ds,overall_range,vperH2Ods)
title('Volume Percentages of Constituents in Dead Space')
xlabel('Time (s)')
ylabel(' Volume Percent (%)')

mass1 = totalmass(TV,vfrac1,M);
% mass1 = the total mass of inspired air

c = 1005;
Tb = 310; % body temperature
Ta = 288; % inspired air temperature
Q12 = thermal(mass1,c,Ta,Tb);
% calculates transfer of thermal energy in stream 12 to inspired air
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

% total mass finds the total mass of inspired air
% TV = total volume of inspired air
% pp = partial pressure of a constituent in inspired air
% Ta = temperature of inspired air
% R = universal gas constant 
% n = moles of constituent in inspired air
% species_mass = the mass of each species in inspired air
% mass1 = the total mass of inspired air

function mass1 = totalmass(TV,vfrac,M)
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

function vPercentOverTime = blood(t_start_6,tin,t_delay_7,texp,tot_t)
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


tresp = tot_t
% calculates length of time per inspiration and expiration from breathing
% frequencies
PaO2i = 100;%100 mmHg partial pressure in beginning for O2
PaCO2i = 40;%40 mmHg partial pressure in beginning for CO2
PH2O = 47;%presures taken from diagram in assumed values
PN2 = 573;
Ptotal = PaO2i+PaCO2i+PH2O+PN2;
vTotalInitial = 3000;
vTotal = vTotalInitial;
nN2i = 573*vTotal/62.36367/1000/310;
nH2Oi = 47*vTotal/62.36367/1000/310;
TV = 0.5; %tidal volume = 500 mL O2 perfect at 0.9
DifCapO2 = 21; %mL/min/mmHg
DifCapCO2 = 400; %mL/min/mmHg
xMax = 3;
dO2 =  density(760, 31.9988, 310); %1 atm, 31.9988 g/mol, 310 K(body temperature)
dCO2 = density(760, 44.0095, 310); %1 atm, 44.0095 g/mol, 310 K(body temperature)
dN2 = density(760, 28.014, 310);
dH2O = density(760, 18.01528 ,310);
trange = 1;
tstep = 0.001;
tsize = length(0:tstep:tresp)

PaO2alveoliOverTime = zeros(1,tsize*xMax);
PaCO2alveoliOverTime = zeros(1,tsize*xMax);
PN2alveoliOverTime = zeros(1, tsize*xMax);
PH2OalveoliOverTime = zeros(1, tsize*xMax);
nO2alveoliOverTime = zeros(1,tsize*xMax);
nCO2alveoliOverTime = zeros(1,tsize*xMax);
nH2OalveoliOverTime = zeros(1, tsize*xMax);
nN2alveoliOverTime = zeros(1, tsize*xMax);
PaO2alveoli = PaO2i;
PaCO2alveoli = PaCO2i;
nN2 = nN2i;
nH2O = nH2Oi;
PaO2capillary = 60; %Partial pressure are really good at 65
PaCO2capillary = 46;
vO2rateoverTime = zeros(1,tsize*xMax);
vCO2rateoverTime = zeros(1,tsize*xMax);
vN2rateoverTime = zeros(1,tsize*xMax);
vH2OrateoverTime = zeros(1,tsize*xMax);
vO2overTime = zeros(1,tsize*xMax);
vCO2overTime = zeros(1,tsize*xMax);
vN2overTime = zeros(1,tsize*xMax);
vH2OoverTime = zeros(1,tsize*xMax);
mO2overTime = zeros(1,tsize*xMax);
mCO2overTime = zeros(1,tsize*xMax);
mN2overTime = zeros(1,tsize*xMax);
mH2OoverTime = zeros(1,tsize*xMax);
vPercentOverTime = zeros(4, tsize*xMax+xMax*2);
PTotalOverTime = zeros(1, tsize*xMax);
nTotalOverTime = zeros(1, tsize*xMax);
DiffusionRateO2OverTime = zeros(1,tsize*xMax);
DiffusionRateCO2OverTime = zeros(1, tsize*xMax);

for x=1:xMax %modeling 1 breath
    %t = t_start_6+tin;
    t = 0;
    mO2inTotal = 0;
    mCO2inTotal = 0;
    
    
    while t<texp
        
        vO2out = volumetricflow(RFex, TV, t)*0.7*tstep*PaO2alveoli/Ptotal;
        vCO2out = volumetricflow(RFex, TV, t)*0.7*tstep*PaCO2alveoli/Ptotal;
        vH2Oout = volumetricflow(RFex, TV, t)*0.7*tstep*0.0344;
        vN2out = volumetricflow(RFex, TV, t)*0.7*tstep*0.7611;
        
        if vO2out > 0
        vO2out = vO2out*-1;
        end
        if vCO2out > 0
        vCO2out = vCO2out*-1;
        end
        if vH2Oout > 0 
        vH2Oout = vH2Oout*-1;
        end
        if vN2out > 0
        vN2out = vN2out*-1;
        end
        
        mO2out = vO2out*dO2;%0.15263;%0.2042;
        mCO2out = vCO2out*dCO2;%0.0421;%*0.0003;%
        mH2Oout = vH2Oout*dH2O;%PH2O/Ptotal;%0.0344;%0.0618;%
        mN2out = vN2out*dN2;%PN2/Ptotal;%0.7611;%0.74342;%
        
        
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
        nO2 = mO2/31.9988;
        nCO2 = mCO2/44.0095;
        nN2diff = mN2out/28.014;
        nH2Odiff = mH2Oout/18.01528;
        
        %PV=nRT, P = nRT/V
        PaO2diff = nO2*62.36367*1000*310/vTotal; %gives pressure in mmHg
        PaCO2diff = nCO2*62.36367*1000*310/vTotal;
        
        PaO2alveoli = PaO2alveoli + PaO2diff;
        PaCO2alveoli = PaCO2alveoli + PaCO2diff;
        
        nN2 = nN2 + nN2diff;
        nH2O = nH2O + nH2Odiff;
        %PN2diff = nN2*62.36367*1000*310/vTotal;
        %PH2Odiff = nH2O*62.36367*1000*310/vTotal;
        
        nO2Total = PaO2alveoli*vTotal/310/62.36367/1000;
        nCO2Total = PaCO2alveoli*vTotal/310/62.36367/1000;
        nGasTotal = nO2Total + nN2 + nH2O + nCO2Total;
        %PN2 = PN2 + PN2diff;
        %PH2O = PH2O + PH2Odiff;
        
        
        %vTotal = nO2Total*62.36367*1000*310/PaO2alveoli;
        
        PH2O = nH2O*62.36367*1000*310/vTotal;
        PN2 = nN2*62.36367*1000*310/vTotal;
        Ptotal = PaO2alveoli + PaCO2alveoli + PH2O + PN2;
        DiffusionRateO2OverTime(trange) = difRateO2;
        DiffusionRateCO2OverTime(trange) = difRateCO2;
        PTotalOverTime(trange) = Ptotal;
        nTotalOverTime(trange) = nGasTotal;
        PaO2alveoliOverTime(trange) = PaO2alveoli;
        PaCO2alveoliOverTime(trange) = PaCO2alveoli;
        PN2alveoliOverTime(trange) = PN2;
        PH2OalveoliOverTime(trange) = PH2O;
        nO2alveoliOverTime(trange)= nO2Total;
        nCO2alveoliOverTime(trange)= nCO2Total;
        nN2alveoliOverTime(trange) = nN2;
        nH2OalveoliOverTime(trange) = nH2O;
        vO2a = nO2Total*62.36367*1000*310/PaO2alveoli;
        vCO2a = nCO2Total*62.36367*1000*310/PaCO2alveoli;
        vN2a = nN2*62.36367*1000*310/PN2;
        vH2Oa = nH2O*62.36367*1000*310/PH2O;
        vO2overTime(trange) = vO2a;
        vCO2overTime(trange) = vCO2a;
        vN2overTime(trange) = vN2a;
        vH2OoverTime(trange) = vH2Oa;
        vO2rateoverTime(trange) = vO2out;
        vCO2rateoverTime(trange) = vCO2out;
        vN2rateoverTime(trange) = vN2out;
        vH2OrateoverTime(trange) = vH2Oout;
        mO2overTime(trange) = mO2out;
        mCO2overTime(trange) = mCO2out;
        mN2overTime(trange) = mN2out;
        mH2OoverTime(trange) = mH2Oout;
        t=t+tstep;
        trange = trange+1;
    end 
    
    t = 0;
    while t<t_start_6+t_delay_7
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
        nO2 = mO2/31.9988;
        nCO2 = mCO2/44.0095;
        %PV=nRT, P = nRT/V
        
        PaO2diff = nO2*62.36367*1000*310/vTotal; %gives pressure in mmHg
        PaCO2diff = nCO2*62.36367*1000*310/vTotal;
        
        nO2Total = PaO2alveoli*vTotal/310/62.36367/1000;
        nCO2Total = PaCO2alveoli*vTotal/310/62.36367/1000;
        nGasTotal = nO2Total + nN2 + nH2O + nCO2Total;
        
        %vTotal = nO2Total*62.36367*1000*310/PaO2alveoli;
        PaO2alveoli = PaO2alveoli + PaO2diff;
        PaCO2alveoli = PaCO2alveoli + PaCO2diff;
        
        
        Ptotal = PaO2alveoli + PaCO2alveoli + PH2O + PN2;
        PTotalOverTime(trange) = Ptotal;
        nTotalOverTime(trange) = nGasTotal;
        DiffusionRateO2OverTime(trange) = difRateO2;
        DiffusionRateCO2OverTime(trange) = difRateCO2;
        PaO2alveoliOverTime(trange) = PaO2alveoli;
        PaCO2alveoliOverTime(trange) = PaCO2alveoli;
        PN2alveoliOverTime(trange) = PN2;
        PH2OalveoliOverTime(trange) = PH2O;
        nO2alveoliOverTime(trange)= nO2Total;
        nCO2alveoliOverTime(trange)= nCO2Total;
        nN2alveoliOverTime(trange) = nN2;
        nH2OalveoliOverTime(trange) = nH2O;
        vO2a = nO2Total*62.36367*1000*310/PaO2alveoli;
        vCO2a = nCO2Total*62.36367*1000*310/PaCO2alveoli;
        vN2a = nN2*62.36367*1000*310/PN2;
        vH2Oa = nH2O*62.36367*1000*310/PH2O;
        vO2overTime(trange) = vO2a;
        vCO2overTime(trange) = vCO2a;
        vN2overTime(trange) = vN2a;
        vH2OoverTime(trange) = vH2Oa;
        vO2rateoverTime(trange) = 0;
        vCO2rateoverTime(trange) = 0;
        vN2rateoverTime(trange) = 0;
        vH2OrateoverTime(trange) = 0;
        
        mO2overTime(trange) = 0;
        mCO2overTime(trange) = 0;
        mN2overTime(trange) = 0;
        mH2OoverTime(trange) = 0;
        
        t=t+tstep;
        trange = trange+1;
    end
    
    t=0;
    
    while t<tin
        vO2in = volumetricflow(RFin, TV, t)*0.7*0.2042*tstep;
        vCO2in = volumetricflow(RFin, TV, t)*0.7*0.0003*tstep;
        vH2Oin = volumetricflow(RFin, TV, t)*0.7*0.0344*tstep;
        vN2in =volumetricflow(RFin, TV, t)*0.7*0.7611*tstep;
        
        if vO2in < 0
        vO2in = vO2in*-1;
        end
        if vCO2in < 0
        vCO2in = vCO2in*-1;
        end
        if vH2Oin < 0 
        vH2Oin = vH2Oin*-1;
        end
        if vN2in < 0
        vN2in = vN2in*-1;
        end
        
        mO2in = vO2in*dO2;
        mCO2in = vCO2in*dCO2;
        mH2Oin = vH2Oin*dH2O;
        mN2in = vN2in*dN2;
        %{
        mO2in = -0.7*(pi*RFin*TV/60)*sin(pi*RFin*t/30)*0.218*dO2*tstep; %insert partial pressures/total pressure or mass frac for humidified air here
        mCO2in = -0.7*(pi*RFin*TV/60)*sin(pi*RFin*t/30)*0.0003*dCO2*tstep;
        mH2Oin = -0.7*((pi*RFin*TV/60)*sin(pi*RFin*t/30))*dH2O*0.0697*tstep;
        mN2in = -0.7*((pi*RFin*TV/60)*sin(pi*RFin*t/30))*dN2*0.78*tstep;
        %}
        
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
        
        nO2 = mO2/31.9988;
        nCO2 = mCO2/44.0095;
        nN2diff = mN2in/28.014;
        nH2Odiff = mH2Oin/18.01528;
        
        %PV=nRT, P = nRT/V
        
        PaO2diff = nO2*62.36367*1000*310/vTotal; %gives pressure in mmHg
        PaCO2diff = nCO2*62.36367*1000*310/vTotal;
        
        nN2 = nN2 + nN2diff;
        nH2O = nH2O + nH2Odiff;
        
        nO2Total = PaO2alveoli*vTotal/310/62.36367/1000;
        nCO2Total = PaCO2alveoli*vTotal/310/62.36367/1000;
        nGasTotal = nO2Total + nN2 + nH2O + nCO2Total;
        
        PaO2alveoli = PaO2alveoli + PaO2diff;
        PaCO2alveoli = PaCO2alveoli + PaCO2diff;
        
        PH2O = nH2O*62.36367*1000*310/vTotal;
        PN2 = nN2*62.36367*1000*310/vTotal;
        Ptotal = PaO2alveoli + PaCO2alveoli + PH2O + PN2;
        DiffusionRateO2OverTime(trange) = difRateO2;
        DiffusionRateCO2OverTime(trange) = difRateCO2;
        PTotalOverTime(trange) = Ptotal;
        nTotalOverTime(trange) = nGasTotal;
        PaO2alveoliOverTime(trange) = PaO2alveoli;
        PaCO2alveoliOverTime(trange) = PaCO2alveoli;
        PN2alveoliOverTime(trange) = PN2;
        PH2OalveoliOverTime(trange) = PH2O;
        vO2a = nO2Total*62.36367*1000*310/PaO2alveoli;
        vCO2a = nCO2Total*62.36367*1000*310/PaCO2alveoli;
        vN2a = nN2*62.36367*1000*310/PN2;
        vH2Oa = nH2O*62.36367*1000*310/PH2O;
        vO2overTime(trange) = vO2a;
        vCO2overTime(trange) = vCO2a;
        vN2overTime(trange) = vN2a;
        vH2OoverTime(trange) = vH2Oa;
        vO2rateoverTime(trange) = vO2in;
        vCO2rateoverTime(trange) = vCO2in;
        vN2rateoverTime(trange) = vN2in;
        vH2OrateoverTime(trange) = vH2Oin;
        mO2overTime(trange) = mO2in;
        mCO2overTime(trange) = mCO2in;
        mN2overTime(trange) = mN2in;
        mH2OoverTime(trange) = mH2Oin;
        nO2alveoliOverTime(trange)= nO2Total;
        nCO2alveoliOverTime(trange)= nCO2Total;
        nN2alveoliOverTime(trange) = nN2;
        nH2OalveoliOverTime(trange) = nH2O;
        
        mO2inTotal = mO2inTotal + mO2in;
        mCO2inTotal = mCO2inTotal + mCO2in;
        t=t+tstep;
        trange = trange+1;
        
    end
end    

resp_range = 0:tstep:length(nO2alveoliOverTime)*tstep-tstep;
vPercentOverTime(1,:)= (nO2alveoliOverTime./nTotalOverTime*100); 
vPercentOverTime(2,:) = (nCO2alveoliOverTime./nTotalOverTime*100);
vPercentOverTime(3,:) = (nH2OalveoliOverTime./nTotalOverTime*100);
vPercentOverTime(4,:) = (nN2alveoliOverTime./nTotalOverTime*100);
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
plot(resp_range, PN2alveoliOverTime)
title('Partial Pressure of N2 in Alveoli Over Time')
xlabel('Time (s)')
ylabel('Partial Pressure (mmHg)')
figure
plot(resp_range, PH2OalveoliOverTime)
title('Partial Pressure of H2O in Alveoli Over Time')
xlabel('Time (s)')
ylabel('Partial Pressure (mmHg)')
figure
plot(resp_range, vO2rateoverTime/tstep);
title('Volumetric Flow Rate of O2 in Streams 6 and 7')
xlabel('Time (s)')
ylabel('Volumetric Flow Rate (L/s)')
figure
plot(resp_range,vCO2rateoverTime/tstep);
title('Volumetric Flow Rate of CO2 in Streams 6 and 7')
xlabel('Time (s)')
ylabel('Volumetric Flow Rate (L/s)')
figure
plot(resp_range, vN2rateoverTime()/tstep);
title('Volumetric Flow Rate of N2 in Streams 6 and 7')
xlabel('Time (s)')
ylabel('Volumetric Flow Rate(L/s)')
figure
plot(resp_range, vH2OrateoverTime()/tstep);
title('Volumetric Flow Rate of H2O in Streams 6 and 7')
xlabel('Time (s)')
ylabel('Volumetric Flow Rate(L/s)')
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
figure
plot(resp_range, vO2overTime); 
title('Volume of O2 in Alveoli Over Time')
xlabel('Time (s)')
ylabel('Volume')
figure
plot(resp_range, vCO2overTime); 
title('Volume of CO2 in Alveoli Over Time')
xlabel('Time (s)')
ylabel('Volume')
figure
plot(resp_range, vN2overTime); 
title('Volume of N2 in Alveoli Over Time')
xlabel('Time (s)')
ylabel('Volume')
figure
plot(resp_range, vH2OoverTime); 
title('Volume of H2O in Alveoli Over Time')
xlabel('Time (s)')
ylabel('Volume')

%{
figure
plot(resp_range, PTotalOverTime);
title('Total Pressure Over Time')
%}
figure
plot(resp_range, DiffusionRateO2OverTime()/tstep);
title('Diffusion Rate of O2 Across Respiratory Membrane Over Time (Stream 8)')
xlabel('Time (s)')
ylabel('Mass Flow Rate(g/s)')
figure
plot(resp_range, DiffusionRateCO2OverTime()/tstep)
title('Diffusion Rate of CO2 Across Respiratory Membrane Over Time (Stream 9)')
xlabel('Time (s)')
ylabel('Mass Flow Rate(g/s)')
figure
%{
plot(resp_range, nTotalOverTime)
title('Total Moles of Gas in Alveoli Over Time')
xlabel('Time (s)')
ylabel('Total Moles')
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
function [t_start_4,t_start_6,Volume_deadspace,Volume_entry] = traveltime(RF,TV,index,range)
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
for j = 1:index
    int_vflow(j) = antiderivative(RF,TV,range(j))- antiderivative(RF,TV,0);
    if abs(int_vflow(j) - 50/1000) < 0.01
        t1 = range(j);
        % the integral of the inspiratory volumetric flow rate from 0 to t1
        % equals the volume of the extrathoracic volume
    end
end
for j = 1:index
    if abs(int_vflow(j)*1000/A(1) - l(1)) < 0.05
        t2 = range(j);
        % int_vflow(i)*1000/(pi*(r(1))^2 is the integral of the 
        % linear velocity
        % this integral from 0 to t2 equals the length of the trachea l(1)
    end
end

Volume_entry = 50 + pi * r(1)^2 * l(1);

t_start_4 = t1 + t2; 
% the time to first reach the dead space unit (time to reach generation 1 
% brochi) is the sum of the time to travel through the extrathoracic region 
% and the time to travel through the trachea

for i = 2:15
    for j =1:index
        if abs(int_vflow(j)*1000/A(i) - l(i)) < 0.01
            t_dead_gen(j) = range(j);
        end
    end
end

t_dead_sum = sum(t_dead_gen);
% the first alveoli appear at generation 15
% the time to travel from the exit of the trachea to the first alveoli is 
% the sum of the time to travel though the generation 1 through 14 airways

t_start_6 = t_start_4 + t_dead_sum;
% t_start_6 = the total time required for air to reach the first alveoli

for i = 2:15
    Volume(i) = pi * r(i)^2 * l(i) * 2^(i-1);
end
Volume_deadspace = sum(Volume);
% calculates the volume of the dead space unit
end

%This function finds the volume compositions of gas constituents in the
%streams where volume composition does not change over time-Stream 1, 4, 6
%vflow= contains all the flow rates of the diagram
%vfrac= contains all the volume fractions
function partial_vol= const_volume(vflow, vfrac)
partial_vol=zeros(6,length(vflow));
for i=1:length(vflow)
    partial_vol=vfrac.*vflow;
end
end

% intflowin = integral of volumetric flow rate of stream going in to unit
% intflowout = integral of volumetric flow rate of stream going out of unit
% **TO-DO**: For expiration streams, change vfrac(i) to vfracexp(i,j) <-- which we get
% from Tony's function
function [VO2,VCO2,VN2,VH2O,Vtot,PO2,PCO2,PN2,PH2O,vperO2,vperCO2,...
    vperN2,vperH2O] = composition(vPercentOverTime,vfrac,intflowin_in,intflowout_in,...
    intflowin_exp,intflowout_exp,RV,index,M_air)
vFractionOverTime(4,:) = vPercentOverTime(3,:) ./ 100;
vFractionOverTime(3,:) = vPercentOverTime(4,:) ./ 100;
vFractionOverTime(1,:) = vPercentOverTime(1,:) ./ 100;
vFractionOverTime(2,:) = vPercentOverTime(2,:) ./ 100;
Tb = 310;
R = 62.3637; 
V = zeros(length(vfrac),index);
V(:,1) = RV .* vfrac';
Vtot = zeros(1,index);
Vtot(1) = sum(RV .* vfrac');
V_temp = zeros(1,index);
for i = 1:length(vfrac)
    for j = 2:index
        V(i,j) = V(i,j-1) - (vfrac(i)*intflowin_in(j)-vfrac(i)*...
            intflowin_in(j-1)) + (vfrac(i)*intflowout_in(j)-vfrac(i)*...
            intflowout_in(j-1)) + (vFractionOverTime(i,j)*intflowin_exp(j)-vFractionOverTime(i,j)*...
            intflowin_exp(j-1)) - (vFractionOverTime(i,j)*intflowout_exp(j)-vFractionOverTime(i,j)*...
            intflowout_exp(j-1));
        % calculates the volume of each constituent i in the unit over time 
    end
end

for i = 1:length(vfrac)
    for j = 1:index
        Vs = V(:,j);
        Vtot(j) = sum(Vs);
    end
end
% calculates the total volume of each constituent in the unit over time

for i = 1:length(vfrac)
    for j = 2:index
        if abs(Vtot(j)) < 0.0001
        vper_s(i,j) = 0;
        % TO-DO: Why isn't this working?
        else
            vper_s(i,j) = V(i,j) ./ Vtot(j) .* 100;
        end
    end
end
% calculates the volume percent of each constituent in the unit over time

vperO2 = vper_s(1,:);
vperCO2 = vper_s(2,:);
vperN2 = vper_s(3,:);
vperH2O = vper_s(4,:);
VO2 = V(1,:);
VCO2 = V(2,:);
VN2 = V(3,:);
VH2O = V(4,:);

density=1.184; 
%density of inhaled air (g/L)
ntot = Vtot .* density ./ M_air;
Ptot = ntot .* R .* Tb ./ Vtot;
for i = 1:index
    if abs(Vtot(i)) < 0.0001
        PO2(i) = 0;
        PCO2(i) = 0;
        PN2(i) = 0;
        PH2O(i) = 0;
    else
        PO2(i) = VO2(i) ./ Vtot(i) .* Ptot(i);
        PCO2(i) = VCO2(i) ./ Vtot(i) .* Ptot(i);
        PN2(i) = VN2(i) ./ Vtot(i) .* Ptot(i);
        PH2O(i) = VH2O(i) ./ Vtot(i) .* Ptot(i);
    end
end
end

% calculates the antiderivative of the volumetric flow rate equation
function output = antiderivative(RF,TV,t)
output = -0.5*TV*cos(pi*RF*t/30);
end

function nflownet = molrate(vflowA,vflowB,vflowC,vfracA,vfracB,vfracC,air_density,M_air,index_overall)
nflowA = vflowA*air_density/M_air;
nflowB = vflowB*air_density/M_air;
nflowC = vflowC*air_density/M_air;
%nflowA = vflowA .* air_density ./ M_air;
%nflowB = vflowB .* air_density ./ M_air;
%nflowC = vflowC .* air_density ./ M_air;
for i = 1:length(vfracA)
   nflows_A(i,:) = nflowA * vfracA(i);
   nflows_B(i,:) = nflowB * vfracB(i);
   %nflows_C(i,:) = nflowC * vfracC(i);  
end
nflownet = nflows_A + nflows_B; %+ nflows_C;


end