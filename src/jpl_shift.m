%
% basic JPL shift function from the L1C ATBD
%

function Tb_resamp = jpl_shift(Tb_in, v_in, v_nom, d1);

Tb_in = Tb_in(:);
v_in = v_in(:);
v_nom = v_nom(:);

a = d1.a; 
b = d1.b;

dv = v_nom - v_in;

% Just uncomment lines under Yibo or L Strow to pick interpolation approach
%
%% Yibo's approach
% Tb_spline = interp1(v_in, Tb_in, v_nom, 'spline');
% 
% Tb_resamp = Tb_in + (a .* (Tb_spline - Tb_in) ./ dv + b) .* dv;

%% L Strow's approach (requires Matlab curve fitting toolbox),
%% uses spline slope for interpolation rather than the noisy slope from data
pp = csape(v_in, Tb_in);
fprime = fnder(pp);
btderiv = fnval(fprime,v_nom);
Tb_resamp = Tb_in + (a.*btderiv + b).*dv;
