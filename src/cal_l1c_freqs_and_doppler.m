function [radiances_nucal_scan, yoff_scan_line] = cal_l1c_freqs_and_doppler(fn);
%
% Returns L1c radiances (90x135) corrected for frequency shift due to the 
% instrument drifts and the Doppler effect.  The returned frequencies are
% sampled to the official L1c frequency scale.

%======== Load in data that will be presistent ==============
% The official l1c channels, we interpolate to this scale
persistent fl1c
if isempty(fl1c)
   load fl1c
end

% Indices for swtiching from l1c to l1b and back
persistent   l1b_ind_in_l1c  l1c_ind_for_l1b
if isempty(l1b_ind_in_l1c ) | isempty(l1c_ind_for_l1b)
   load indices_of_l1b_in_l1c
end
%======================== Read L1c File Vars Needed =====================

%% No error trapping yet for partial granules, not sure if this is an issue
%% If need to trap, start by getting "state" variable below, and then only
%% read variables from 1:n0, etc etc.  Maybe needed for "bad" granules too?
%% L Strow can implement if you give me some bad granules for testing
% junk = hdfread(fn, 'state');
% state = reshape( double(junk'), 1,nobs);
% i0=find( state == 0);  % Indices of "good" FOVs
% n0=length(i0);

% Check to see if this is a short granule
num_scanlines = cell2mat(hdfread(fn,'num_scanlines'));

% Granule length
nobs = 90*num_scanlines;
nxtrack = 90;
natrack = num_scanlines;
nobs = nxtrack*natrack;
nchan = 2645;

[xtr, atr] = meshgrid(1:nxtrack,1:natrack);
lxtr = reshape(xtr',nobs,1);
latr = reshape(atr',nobs,1);

% Per scan line variables
junk = cell2mat(hdfread(fn,'scan_node_type'));
junk = (junk == 68); % snode = 0 is asc, = 1 is desc
% junk(130:135) = 1;  % Testing
junk  = repmat(junk(:),1,90);
scan_node_type = reshape(junk',nobs,1);

junk = cell2mat(hdfread(fn,'sat_lat'));
junk = repmat(junk(:),1,90);
sat_lat = reshape(junk',nobs,1);

% 90 x 135 variables
junk = hdfread(fn,'Latitude');
Latitude = reshape(junk',nobs,1);

junk = hdfread(fn,'Longitude');
Longitude = reshape(junk',nobs,1);

junk = hdfread(fn,'Time');
rtime = reshape( junk',nobs,1);
rtime = rtime  + 12784 * 86400 + 27;  % TAI time

junk = hdfread(fn,'satzen');
satzen = reshape(junk',nobs,1);

junk = hdfread(fn,'satazi');
satazi = reshape(junk',nobs,1);

mtime = tai2dtime(rtime);

junk = hdfread(fn, 'radiances');
junk = permute(junk,[3 2 1]);
radiances  = reshape(junk,nchan,nobs)';

clear junk
%=========================== Get Orbit Phase ==========================
opi = get_opi(Latitude,scan_node_type);
%=========================== Get Obs Frequency ==========================
% Get indices into yoff matrix which handles 1:180 for opi
% No interpolation of orbit phase, just use closest of 180 phases in table
[c,ia,ib] = unique(opi);
% Will not worry about ab changing during a granule 
ab_time = get_ab_state(nanmean(mtime));
yoff = get_yoff(nanmean(mtime));

% Get mean yoff per scanline for Evan
for i=1:natrack
   k = find((lxtr == 45 | lxtr == 46) & latr == i);
   scan_line_opi(i) = nanmean(opi(k));
end
yoff_scan_line = yoff(:,round(scan_line_opi));

% Get freq, the actual (computed) frequencies for the observation at mtime
for i=1:length(c)
   % Only do gmodel on unique orbit phases, fill these back to all scenes   
   [f_lm,freq(i,:),m_lm,module] = gmodel(155.1325,yoff(:,opi(ia(i))),ab_time);
end
% Fill in freqs for all scenes
tmp_freqall = freq(ib,:);

% % Only shift true l1b channels in l1c
% freqall = freqall(:,l1b_ind_in_l1c);
%=========================== Get Doppler Shift ==========================
dnu_ppm = doppler_jpl(scan_node_type,lxtr,satzen,satazi,sat_lat);
dnu     = dnu_ppm*1E-6.*tmp_freqall;
tmp_freqall = tmp_freqall + dnu;

[nx ny ] = size(radiances);

% Convert tmp_freqall to L1c channels, fill with standard L1c freqs first to fill fake channels
freqall = zeros(nx,ny);
freqall = freqall + fl1c';
% Now fill in frequencies for l1b channels from grating model
freqall(:,l1c_ind_for_l1b) = tmp_freqall(:,l1b_ind_in_l1c);

% rad2bt slows way down (or memory allocation) when mixing in any complex radiances
% The following codes separates real and complex BT creation.  If a granule has only
% 2000 negative radiances, the slow down in the commented out loop below is 100X!!
% 
% Removed code due to slow rad2bt on negative radiances
% l1b_btobs = NaN(nx,ny);
% profile on;nx = 1000;
% for i = 1:nx
%   l1b_btobs(i,:) = rad2bt(freqall(i,:),radiances(i,:));
% end

% Fast rad2bt for comples, use indices k, ki, kj for reverse operation (complex_bt2rad)
[l1b_btobs k kr] = complex_rad2bt(freqall,radiances);

% Lot's of programming to avoid mixed real/imag returns from jpl_shift.m
% Some kind of Matlab performance "bug"
if isreal(l1b_btobs)
   tmp_btobs = NaN(nx,ny);
   for i=1:nx
      tmp_btobs(i,:) = jpl_shift(l1b_btobs(i,:),freqall(i,:),fl1c);
   end
else
   tmp_btobs = complex(NaN(nx,ny));   
% Cut this call up into all real and imag
   for i= 1:nx
      r_or_i(i) = isreal(l1b_btobs(i,:));
   end
% First do imaginary
   indi = find( r_or_i == 0 );
   tmpi = complex(NaN(length(indi),ny));
   j = 1;
   for i=indi
      tmpi(j,:) = jpl_shift(l1b_btobs(i,:),freqall(i,:),fl1c);
      j = j + 1;
   end
% Now do real
   indr = find( r_or_i == 1 );
   tmpr = NaN(length(indr),ny);
   j = 1;
   for i=indr
      tmpr(j,:) = jpl_shift(l1b_btobs(i,:),freqall(i,:),fl1c);
      j = j + 1;
   end
   tmp_btobs(indi,:) = tmpi;
   tmp_btobs(indr,:) = tmpr;
end

% Matlab will preserve negative radiances, says testing, do NOT ask for real part of tmp_btobs
% Plus special call to avoid slow downs, same issue as for rad2bt
radiances_nucal = complex_bt2rad(fl1c,tmp_btobs,k,kr);

radiances_nucal_scan = reshape(radiances_nucal,90,num_scanlines,2645);
radiances_nucal_scan = permute(radiances_nucal_scan,[3 1 2]);
