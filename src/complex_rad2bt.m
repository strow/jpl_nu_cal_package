function [bt k ki kj] = complex_rad2bt(fa,ra);

[nx ny ] = size(ra);

% Get negative radiances in both logical and index form
k = ra < 0;
[ki,kj] = find( ra < 0);

% Just do the small number of negative radiances
if length(ki) > 0
   x = NaN(1,length(ki));
   for i=1:length(ki)
      x(i) = rad2bt(fa(ki(i),kj(i)),ra(ki(i),kj(i)));
   end
end;
ra(ki,kj) = NaN;   % Get rid of negative radiances for next loop

% Now do the positive radiances
bt = NaN(nx,ny);
for i = 1:nx
   bt(i,:) = rad2bt(fa(i,:),ra(i,:));
end
bt = complex(bt);

% If there are any negative radiances, their bt in in x
if length(ki) > 0
   bt(k) = x;
end
