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
Tb_spline = interp1(v_in, Tb_in, v_nom, 'spline');

% Trap fill channels since dv == 0
kf = dv == 0;
k  = dv ~= 0;

% Put back input for fill channels.
% Remove line below so don't divide by dv(k)  (OLD)
% Tb_resamp(k) = Tb_in(k) + (a(k) .* (Tb_spline(k) - Tb_in(k)) ./ dv(k) + b(k)) .* dv(k);

% % At JPL for testing in March 25, 2019: WRONG
Tb_resamp(k) = Tb_in(k) + a(k) .* (Tb_spline(k) - Tb_in(k))  + b(k) .* dv(k);
Tb_resamp(kf) = Tb_in(kf);

