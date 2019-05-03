function [bt, k, kr] = complex_rad2bt(fa,ra);

[nx ny ] = size(ra);

% Get negative radiances in both logical and index form
k  = ra <  0;
kr = ra >= 0;

% Pre-allocate so can use k, kr
bt     = NaN(nx,ny);

bt(kr) = rad2bt(fa(kr),ra(kr));

% Next rad2bt will be complex, so set that up
bt  = complex(bt);

% Complex bt calcs from negative radiances
bt(k) = rad2bt(fa(k),ra(k));
