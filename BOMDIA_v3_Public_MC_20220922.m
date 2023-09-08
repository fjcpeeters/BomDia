% =============================================================
% BOM DIA version 3: An algorithm for age modelling of marine hemi-pelagic
% sediments using CaCO3 wt.percent as a proxy for sedimentation rate
%
%--------------------------------- Version 3 ------------------------------
% ------------------------------ August 18th 2022 ---------------------------
% ------------------ With Monte Carlo uncertainty estimation --------------


% ==============================================================
% Dr. Frank Peeters
% Department of Earth Sciences
% Faculty of Science, Vrije Universiteit
% De Boelelaan 1085, 1081 HV Amsterdam,
% The Netherlands. email: f.j.c.peeters@vu.nl
% =============================================================

% Two input files are required: BDv3_input_CaCO3 and BDv3_input_ACP
% BDv2_input_CaCO3:
%   1st column equal spaced depth (must start at 0 cm!)
%   2nd column CaCO3wt% (or other proxy for sed rate (between 0 - 100%)
% BDv3_input_ACP:
%   1st column: depth in core 
%   2nd column lower boundary age (use kyr)
%   3rd column median age (use kyr)
%   4th column upper boundary age (use kyr)

clear, clc % clear workspace and command window
close all

%% SECTION 1: SELECT AND LOAD DATA

fprintf('For GeoB3910 - Jaeschke etal. 2007 Paleoc_V22 ------- press 1 \n')
fprintf('For 64PE304_C80 - van der Lubbe etal. 2014 ---------- press 2 \n')
flag_sample = input('ENTER YOUR SELECTION: \n');

fprintf(' How many MC simulations do you want? Choose between 10000 or 500000')
nsim = input(' ENTER YOUR SELECTION: \n');

fprintf(' Choose number of best age models choose 30 or not more than 100')
nbest = input(' ENTER YOUR SELECTION: \n');

fprintf('Add additional uncertainty to ACPs? e.g. for 2 percent enter 2')
addun = input(' ENTER YOUR SELECTION: \n');

if addun == 0
    disp('no additional uncertainty was added')
elseif addun > 10
    disp('you have added to much uncertainty')
    disp ('choose between 0 and 10')
    disp('The script terminates. Press run to restart')
    return
else
   x = ['you have added ', num2str(addun),' percent uncertainty'];
   disp(x)
end

if flag_sample == 1
    load BDv3_input_CaCO3_GeoB3910; % record must start at 0 cm!
    BDv3_input_CaCO3 = BDv3_input_CaCO3_GeoB3910;
    load BDv3_input_ACP_GeoB3910 % load depth and age of ACP's
    BDv3_input_ACP = BDv3_input_ACP_GeoB3910;
    fname = 'GeoB3910-2 Jaeschke et al. (2007)';
    clear BDv3_input_ACP_GeoB3910 BDv3_input_CaCO3_GeoB3910
elseif flag_sample == 2
    load BDv3_input_CaCO3_64PE304_C80; % record must start at 0 cm!
    BDv3_input_CaCO3 = BDv3_input_CaCO3_64PE304_C80;
    load BDv3_input_ACP_64PE304_C80 % load depth and age of ACP's
    BDv3_input_ACP = BDv3_input_ACP_64PE304_C80;
    fname = '64PE304-80 van der Lubbe et al.,(2014)';
    clear BDv3_input_ACP_64PE304_C80 BDv3_input_CaCO3_64PE304_C80
else
    fprintf ('please enter a number from the list above only')
    flag_sample = input('ENTER YOUR SELECTION: \n');
end
tic

%% SECTION 2: ASSIGN INPUT DATA TO VARIABLES

ACP_depth = BDv3_input_ACP(:,1); % store depth of ACP
depth = BDv3_input_ACP(:,1);
agemin= BDv3_input_ACP(:,2);
agemed = BDv3_input_ACP(:,3);
agemax= BDv3_input_ACP(:,4);
nacp = length(depth);

ageminorig = agemin;
agemedorig = agemed;
agemaxorig = agemax;

% increase the ACP uncertainty by some percentage of the absolute value
if addun > 0 % if added uncertainty is lager than zero then
    agemin = agemin - (addun/100).*agemin; 
    agemax = agemax + (addun/100).*agemax; 
end

%% SECTION 3: VARIABLE DECLARATION

dZ = BDv3_input_CaCO3(2,1)-BDv3_input_CaCO3(1,1); % determine depth spacing
indx = ACP_depth*(1/dZ)+1; % index numbers for ACPdepth in Z and CaCO3
Z = BDv3_input_CaCO3(indx(1):indx(end),1); %  store depth array
Z = Z(indx(1):indx(end)); % store depth values up to the depth of last ACP
CaCO3 = BDv3_input_CaCO3(indx(1):indx(end),2); % store CaCO3 array
CaCO3_marker = zeros(length(ACP_depth),1);
ns = length(ACP_depth)-1; % The number of segments
ni = length(CaCO3); % Length of the CaCO3 record

% Define structure variables containing output for segments and agemodelMC
sgmnt = struct;
AgemodelMC = struct;

AgemodelMC.all = zeros(ni,nsim);
AgemodelMC.carMC = zeros(ni,nsim);
Agemodel = zeros(ni,1); % this is a temporary variable
SR = zeros(ni-1,1); % Sedimentation rate
SRMC = zeros(ni-1,nsim); % Sedimentation rate
Age_for_SR = zeros(ni-1,1); % Age for sedimentation rate
Z_for_SR = zeros(ni-1,1); % middepth for SR

CaCO3original = CaCO3;
% CaCO3 = CaCO3.^0.8;

%% SECTION 4: Monte Carlo ACP CREATION (EXCLUDING AGE REVERSALS)
ACPrandomset = zeros(nacp,nsim);
iter = 1;
itertotal = 0;
rand('seed',0)
while iter <= nsim && itertotal <9999999 % safe while loop?
    ACPrandomset(:,iter) = agemin(:)+(agemax(:) - agemin(:)).*rand(nacp,1);
    Diff_ACP_age = diff(ACPrandomset(:,iter));
    contains_reversal = any(Diff_ACP_age<=0);
        if contains_reversal == 1
            iter = iter+0;
            itertotal = itertotal+1;
            X = [num2str(itertotal),' reversal'];
            disp(X)
        else
            iter = iter+1;
            itertotal = itertotal+1;
            X = [num2str(itertotal),' no reversal'];
            disp(X)                
        end
end

%% SECTION 5: CALCULATE CAR VALUES FOR THE SEGMENTS

for k =1:nsim
    ACP_age = ACPrandomset(:,k); % for median ACP age
for s = 1:ns
    age_top = ACP_age(s); age_bottom = ACP_age(s+1);
    CaCO3temp = CaCO3(indx(s):indx(s+1));
    [AM,SR,car] = getcar(CaCO3temp,dZ,age_top,age_bottom);
    sgmnt.car(s) = car;
    Agemodel(indx(s):indx(s+1)-1) = AM(1:end-1);
end

Agemodel(end) = ACP_age(end);

for i = 1:ni-1
    SR(i) = dZ./(Agemodel(i+1)-Agemodel(i));
    Age_for_SR(i) = (Agemodel(i)+Agemodel(i+1))/2;
    Z_for_SR(i) = (Z(i)+Z(i+1))/2;
end

for s = 1:ns+1
    CaCO3_marker(s) = CaCO3original(indx(s)); % generate markers to show ACP's
end

for s = 1:ns
    car(indx(s):indx(s+1))= sgmnt.car(s);
end

AgemodelMC.all(:,k) = Agemodel(:);
SRMC(:,k) = SR;
AgemodelMC.carMC(:,k) = car';
sgmnt.carMC(:,k) = sgmnt.car;
end

sgmnt.carMCstd = std(sgmnt.carMC(1:end,:));

% locate the indices of Monte Carlo age models with lowest variability in
% the downcore CAR parameter
[carMCloweststd,indxcarMCstd] = mink(sgmnt.carMCstd,nbest); 
sgmnt.carMCbest=sgmnt.carMC(:,indxcarMCstd);
AgemodelMC.nsim = nsim; 
AgemodelMC.nbest = AgemodelMC.all(:,indxcarMCstd);
AgemodelMC.nbeststd = std(AgemodelMC.nbest,0,2);
AgemodelMC.LB5 = prctile(AgemodelMC.nbest,5,2);
AgemodelMC.UB95 = prctile(AgemodelMC.nbest,95,2);
AgemodelMC.median = median(AgemodelMC.nbest,2);
AgemodelMC.mean = mean(AgemodelMC.nbest,2);
for i =1:ni-1
AgemodelMC.SRmedian(i) = dZ/(AgemodelMC.median(i+1)-AgemodelMC.median(i));
AgemodelMC.SRmean(i) = dZ/(AgemodelMC.mean(i+1)-AgemodelMC.mean(i));
end

AgemodelMC.SRmedian = AgemodelMC.SRmedian';
AgemodelMC.SRmean = AgemodelMC.SRmean';

midpointCaCO3values = zeros(ni-1,1);

for i = 1:ni-1
    midpointCaCO3values(i) = (CaCO3original(i)+CaCO3original(i+1))/2;
end

AgemodelMC.carmean = AgemodelMC.SRmean.*(midpointCaCO3values/100);

toc

%% SECTION 6: CLEAN UP AND SAVE OUTPUTFILE

% clear unsed variables in workspace
clear i s ns car Agemodel
% save output to file
save('AgemodelBomDiaMC.mat','-struct','AgemodelMC');


%% SECTION 7: PLOT THE RESULTS

figure(1)
set(0, 'DefaultLineLineWidth', 1)
set(0,'defaultAxesFontSize',10)

% PLOT 1 and 2
subplot(3,3,1:2)
plot(Z, CaCO3original,'Color',[0 0.4470 0.7410], 'LineWidth', 1.5)
hold on
plot(ACP_depth,CaCO3_marker,'ro')
hold off
% set(gca,'Ydir','reverse')
ylabel ('CaCO_3 [wt%]')
title (fname),
xlabel ('Depth [cm]')

% PLOT 3
if nbest>1
        subplot(3,3,3)
    for k=2:nbest
        plot(Z,AgemodelMC.all(:,indxcarMCstd(:)),'Color',[0 0.4470 0.7410])
        hold on
    end
        % plot(Z,AgemodelMC.all(:,indxcarMCstd(1)),'-r', 'LineWidth', 2)
        plot(ACP_depth, agemedorig,'+k')
        plot(ACP_depth, ageminorig,'^k')
        plot(ACP_depth, agemaxorig,'vk')
        hold off
 else
        subplot(3,3,3)
        plot(Z,AgemodelMC.all(:,indxcarMCstd(1)),'-r', 'LineWidth', 2)
        hold on
        plot(ACP_depth, agemedorig,'+k')
        plot(ACP_depth, ageminorig,'^k')
        plot(ACP_depth, agemaxorig,'vk')
        hold off
end

if addun == 0
    X = ('Age model');
else
    X = [num2str(addun), ' percent additional uncertainty added to the ACPs'];
end
% title(X),
xlabel ('Depth [cm]')
ylabel ('Age [kyr BP]')

% PLOT 4 and 5
subplot(3,3,4:5)
plot(Z,AgemodelMC.carMC(:,indxcarMCstd(:)),'Color',[0 0.4470 0.7410])
hold on
plot(Z_for_SR,AgemodelMC.carmean,'r','LineWidth', 2)
hold off
ylim ([0 inf])
xlabel ('Depth [cm]'), ylabel ('CAR [cm/kyr]')

% PLOT 6
subplot(3,3,6)
h1 = histogram(sgmnt.carMC(:,indxcarMCstd(:)));
h1.FaceColor = [0 0.4470 0.7410];
h1.EdgeColor = 'k';
% h1.BinWidth = 0.25;
% title('CAR distribution')
xlabel ('Segment CAR value [cm/kyr]')
ylabel ('Frequency')
xlim ([0 inf])

% PLOT 7 and 8
subplot(3,3,7:8)
plot(Z_for_SR,SRMC(:,indxcarMCstd(:)),'Color',[0 0.4470 0.7410])
hold on
plot(Z_for_SR,AgemodelMC.SRmean,'-r','LineWidth', 1.5)
hold off
ylabel ('SAR [cm/kyr]')
xlabel ('Depth [cm]')
ylim([0 inf])
% title (['Sedim. accum. rate (SAR) for best ',num2str(nbest),' out of ',...
%    num2str(nsim),' MC solutions'])

% PLOT 9
subplot(3,3,9)
h2 = histogram(SRMC(:,indxcarMCstd(:)));
h2.FaceColor = [0 0.4470 0.7410];
h2.EdgeColor = 'k';
h2.BinWidth = 2;
xlabel ('SAR [cm/kyr]')
ylabel('Frequency')
ylim([0 inf])
% title('SAR distribution')

%%
figure(2)
if nbest>1
    for k=2:nbest
        plot(Z,AgemodelMC.all(:,indxcarMCstd(k)),'Color',[0 0.4470 0.7410])
        hold on
    end
    plot(Z,AgemodelMC.all(:,indxcarMCstd(1)),'-r', 'LineWidth', 2)
    plot(ACP_depth, agemedorig,'+k')
    plot(ACP_depth, ageminorig,'^k')
    plot(ACP_depth, agemaxorig,'vk')
    hold off
else
    plot(Z,AgemodelMC.all(:,indxcarMCstd(1)),'-r', 'LineWidth', 2)
    hold on
    plot(ACP_depth, agemedorig,'+k')
    plot(ACP_depth, ageminorig,'^k')
    plot(ACP_depth, agemaxorig,'vk')
    hold off
end
title(['Best ', num2str(nbest),' (out of ', num2str(nsim),')', ' MC age depth models'])
xlabel ('Depth [cm]'), ylabel ('Age [kyr BP]')

midpointCaCO3values = zeros(ni-1,1);
ni = length(CaCO3); % Length of the CaCO3 record

for i = 1:ni-1
    midpointCaCO3values(i) = (CaCO3original(i)+CaCO3original(i+1))/2;
end

%%
figure(3)
plot(midpointCaCO3values, SRMC(:,indxcarMCstd(:)),'+','Color',[0 0.4470 0.7410])
title('Scatter plot of SAR versus CaCO3 for best age model')
xlabel ('CaCO3 wt%'), ylabel ('SAR for best agemodel')

x = ([' The lowest standard deviation of CAR values is ',...
    num2str(carMCloweststd(1))]);
disp(x)

%%
figure(4)
subplot(6,1,1:3)
[lineh, bandsh] = fanChart(Z, AgemodelMC.nbest, 'median', 5:5:95, ...
    'alpha', .2, 'colormap', {'shadesOfColor', [0 0 .8]});
hold on
    plot(ACP_depth, agemedorig,'+k')
    plot(ACP_depth, ageminorig,'^k')
    plot(ACP_depth, agemaxorig,'vk')
    plot(Z,AgemodelMC.mean,'r','LineWidth', 1.5)
hold off

grid on
title(fname)
xlabel ('Depth [cm]'),
ylabel ('Age [kyr BP]')

subplot(6,1,4)
plot(Z,2.*AgemodelMC.nbeststd,'k', 'LineWidth', 1.5)
grid on
ylabel('two sigma uncertainty [kyr]')
xlabel('Depth [cm]')

subplot(6,1,5)
plot(Z_for_SR,AgemodelMC.SRmean,'k', 'LineWidth', 1.5)
grid on
xlabel('Depth [cm]')
ylabel('SAR [cm/kyr]')

subplot(6,1,6)
plot(Z,CaCO3, 'LineWidth', 1.5)
set(gca,'Ydir','reverse')
grid on
ylabel('CaCO_3 wt%')
xlabel('Depth [cm]')
 
%% DISPLAY INFO BOMDIA OUTPUT
disp (' ')
disp ('================================================================= ')
disp ('INFO BOMDIA OUTPUT')
disp ('================================================================= ')
disp (' ')
disp('The BomDia age-depth model output is stored in "AgemodelMC"')
disp (' ')
disp('AgemodelMC.all contains all MC produced age-depth models.')
disp (' ')
disp('AgemodelMC.carMC contains all CAR values for all MC produced age-depth models.')
disp (' ')
disp('AgemodelMC.nism contains number of MC simulations.')
disp (' ')
disp('AgemodelMC.nbest contains nbest age-depth models.')
disp (' ')
disp('AgemodelMC.nbeststd contains the 1 s.d. uncertainty.')
disp (' ')
disp('AgemodelMC.LB5 5% percentile lower age boundary.')
disp (' ')
disp('AgemodelMC.UB95 95% percentile upper age boundary.')
disp (' ')
disp('AgemodelMC.median contains the median age-depth model.')
disp (' ')
disp('AgemodelMC.mean contains the mean age-depth model.')
disp (' ')
disp('AgemodelMC.SRmedian contains the sediment accumulation rate for the median age-depth model.')
disp (' ')
disp('AgemodelMC.SRmean contains the sediment accumulation rate for the mean age-depth model.')