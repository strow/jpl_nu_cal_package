function ra = complex_bt2rad(fa,bt,k,ki,kj);

% Here assume fa is a vector (fl1c)

[nx ny ] = size(bt);

% Any negative radiances?
if length(ki) > 0
   x = NaN(1,length(ki));
   for i=1:length(ki)
      x(i) = bt2rad(fa(kj(i)),bt(ki(i),kj(i)));
   end
end;
bt(ki,kj) = NaN;   % Get rid of negative radiances for next loop

% Now do the positive radiances
ra = NaN(nx,ny);
for i = 1:nx
   ra(i,:) = bt2rad(fa,bt(i,:));
end
ra = complex(ra);

% If there were any negative radiances, their bt is in x
if length(ki) > 0
   ra(k) = x;
end
