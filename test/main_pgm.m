% main_pgm.m
%
% Main program to test L1c radiance frequency calibration and Doppler correction

addpath ../src
addpath ../static

% Level 1c file name
fn = 'AIRS.2017.11.09.226.L1C.AIRS_Rad.v6.1.2.0.G17314105227.hdf';

% Frequency Calibrated Radiances
tic
radiances_cal = cal_l1c_freqs_and_doppler(fn);
toc

disp('Done calibrating radiances')
disp('Now doing plotting, etc.')

%==================== Rest of this File is Validation/Plots ===================
% Granule length
nobs = 90*135;
nxtrack = 90;
natrack = 135;
nobs = nxtrack*natrack;
nchan = 2645;

% Get existing L1c radiances for comparison
junk = hdfread(fn, 'radiances');
junk = permute(junk,[3 2 1]);
radiances  = reshape(junk,nchan,nobs)';
clear junk;

[xtr, atr] = meshgrid(1:nxtrack,1:natrack);
lxtr = reshape(xtr',nobs,1);
latr = reshape(atr',nobs,1);

% 90 x 135 variables
junk = hdfread(fn,'Latitude');
Latitude = reshape(junk',nobs,1);

junk = hdfread(fn,'Longitude');
Longitude = reshape(junk',nobs,1);

junk = hdfread(fn,'satzen');
satzen = reshape(junk',nobs,1);

radiances_cal  = reshape(radiances_cal,nchan,nobs)';

% Convert
[nx ny ] = size(radiances);

% Get fl1c 
load fl1c

btobs = NaN(nx,ny);   
btobs_cal = NaN(nx,ny);
for i = 1:nx
 btobs(i,:) = rad2bt(fl1c,radiances(i,:));
 btobs_cal(i,:) = rad2bt(fl1c,radiances_cal(i,:));
end

% Pick 2 channels on different sides of a line
ch1 = 300;
ch2 = 302;

figure
plot(btobs(:,ch1)-btobs_cal(:,ch1))
hold on;
plot(btobs(:,ch2)-btobs_cal(:,ch2))
grid;
xlabel('Scene Linear Index')
ylabel('L1c - L1c Cal in K')
hl = legend([num2str(fl1c(ch1),5) ' cm-1'],[num2str(fl1c(ch2),5) ' cm-1']);

figure
plot(lxtr,btobs(:,ch1)-btobs_cal(:,ch1),'.')
hold on;
plot(lxtr,btobs(:,ch2)-btobs_cal(:,ch2),'.')
grid;
xlabel('Xtrack');
ylabel('L1c - L1c Cal in K')
xlim([0 91])
hl = legend([num2str(fl1c(ch1),5) ' cm-1'],[num2str(fl1c(ch2),5) ' cm-1']);

figure
plot(Latitude,btobs(:,ch1)-btobs_cal(:,ch1),'.')
hold on;
plot(Latitude,btobs(:,ch2)-btobs_cal(:,ch2),'.')
grid;
ylabel('L1c - L1c Cal in K')
xlabel('Latitude');
hl = legend([num2str(fl1c(ch1),5) ' cm-1'],[num2str(fl1c(ch2),5) ' cm-1']);

figure
scatter(Longitude,Latitude,30,btobs(:,ch1)-btobs_cal(:,ch1),'filled')
grid;box on;
colorbar
xlabel('Longitude');
ylabel('Latitude');
title([num2str(fl1c(ch1),5) ' cm-1  (L1c - L1c cal) in K'])

figure
scatter(Longitude,Latitude,30,btobs(:,ch1),'filled')
grid;box on;
colorbar
xlabel('Longitude');
ylabel('Latitude');
title([num2str(fl1c(ch1),5) ' cm-1  L1c B(T)'])

figure
scatter(Longitude,Latitude,30,btobs(:,ch2)-btobs_cal(:,ch2),'filled')
grid;box on;
colorbar
xlabel('Longitude');
ylabel('Latitude');
title([num2str(fl1c(ch2),5) ' cm-1  (L1c - L1c cal) in K'])

figure
scatter(Longitude,Latitude,30,btobs(:,ch2),'filled')
grid;box on;
colorbar
xlabel('Longitude');
ylabel('Latitude');
title([num2str(fl1c(ch2),5) ' cm-1  L1c B(T)'])

figure
plot(fl1c,nanmean(btobs)-nanmean(btobs_cal))
grid
xlabel('Wavenumber')
ylabel('L1c - L1c cal in K')
xlim([630 2700])
title('Granule Mean B(T) of (L1c - L1c cal)')

