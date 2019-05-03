%
% Basic JPL shift function from the L1C ATBD
%
function [Tb_resamp] = jpl_shift(Tb_in, v_in, v_nom);

persistent a b
if isempty(a) | isempty(b)
   load('umbc_shift_1c','a','b');
end

Tb_in = Tb_in(:);
v_in = v_in(:);
v_nom = v_nom(:);

dv = v_nom - v_in;
Tb_spline = interp1(v_in, Tb_in, v_nom, 'pchip');
Tb_resamp = Tb_in + (a .* (Tb_spline - Tb_in) ./ dv + b) .* dv;

