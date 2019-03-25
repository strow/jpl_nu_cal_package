%
% Basic JPL shift function from the L1C ATBD
%
function Tb_resamp = jpl_shift(Tb_in, v_in, v_nom, d1);

Tb_in = Tb_in(:);
v_in = v_in(:);
v_nom = v_nom(:);

a = d1.a; 
b = d1.b;

dv = v_nom - v_in;

Tb_spline = interp1(v_in, Tb_in, v_nom, 'spline');
Tb_resamp = Tb_in + (a .* (Tb_spline - Tb_in) ./ dv + b) .* dv;
