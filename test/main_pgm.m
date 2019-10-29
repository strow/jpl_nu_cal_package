% main_pgm.m
%
% Main program to test L1c radiance frequency calibration and Doppler correction

addpath ../src
addpath ../static

% Level 1c file name
% Granule below is descending with equator passing, failed on earlier bad get_opi.m function
%fn = '/asl/data/airs/L1C/2002/249/AIRS.2002.09.06.176.L1C.AIRS_Rad.v6.1.2.0.G16315005300.hdf';
%fn = '/asl/data/airs/L1C/2002/249/AIRS.2002.09.06.160.L1C.AIRS_Rad.v6.1.2.0.G16315004354.hdf';

%fn = '/asl/data/airs/L1C/2019/001/AIRS.2019.01.01.009.L1C.AIRS_Rad.v6.1.2.0.G19001103631.hdf';
fn = 'AIRS.2017.11.09.226.L1C.AIRS_Rad.v6.1.2.0.G17314105227.hdf';
%fn = 'AIRS.2019.01.01.001.L1C.AIRS_Rad.v6.1.2.0.G19001103741.hdf';
%fn = 'shortgran.hdf';
% Frequency Calibrated Radiances
tic
[radiances_cal, yoff_scan_line] = cal_l1c_freqs_and_doppler(fn);
toc

disp('Done calibrating radiances')
disp('Now doing plotting, etc.')

%==================== Rest of this File is Validation/Plots ===================
% Granule length
num_scanlines = cell2mat(hdfread(fn,'num_scanlines'));
natrack = num_scanlines;
nxtrack = 90;
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

% 90 x natrack variables
junk = hdfread(fn,'Latitude');
Latitude = reshape(junk',nobs,1);

junk = hdfread(fn,'Longitude');
Longitude = reshape(junk',nobs,1);

junk = hdfread(fn,'satzen');
satzen = reshape(junk',nobs,1);

junk = hdfread(fn,'state');
state = reshape(junk',nobs,1);

radiances_cal  = reshape(radiances_cal,nchan,nobs)';

% Convert
[nx ny ] = size(radiances);

% Get fl1c 
load fl1c

% btobs = NaN(nx,ny);   
% btobs_cal = NaN(nx,ny);
% for i = 1:nx
%  btobs(i,:) = rad2bt(fl1c,radiances(i,:));
%  btobs_cal(i,:) = rad2bt(fl1c,radiances_cal(i,:));
% end
fl1c_mat = repmat(fl1c',nobs,1);

btobs = complex_rad2bt(fl1c_mat,radiances);
btobs_cal = complex_rad2bt(fl1c_mat,radiances_cal);

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

% % For L. Strow environment only
% addpath /asl/matlib/plotutils/
% figure(1);
% aslprint('btdiff_vs_index')
% 
% figure(2)
% aslprint('btdiff_vs_xtrack')
% 
% figure(3);
% aslprint('btdiff_vs_latitude')
% 
% figure(4);
% aslprint('ch1_btdiff_map',1)
% 
% figure(5)
% aslprint('ch1_bt_map',1)
% 
% figure(6)
% aslprint('ch2_btdiff_map',1)
% 
% figure(7)
% aslprint('ch2_bt_map',1)
% 
% figure(8)
% aslprint('granule_mean_btdiff')
