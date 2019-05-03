function ra = complex_bt2rad(fl1c,bt,k,kr);

% Here assume fa is a vector (fl1c)

[nx ny ] = size(bt);

% Convert fa to matrix
fa = ones(nx,ny);
for i = 1:nx
   fa(i,:) = fl1c;
end

ra = NaN(nx,ny);

ra(kr) = bt2rad(fa(kr),bt(kr));

ra(k) = bt2rad(fa(k),bt(k));
