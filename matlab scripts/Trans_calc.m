function Trans_calc
  close all; clear all; clc;
  x_step = 0.283; %monolayer width
  x1 = input('first barrier start coordinate, ML: ');
  x2 = input('first barrier end coordinate, ML: ');
  x3 = input('second barrier start coordinate, ML: ');
  x4 = input('second barrier end coordinate, ML: ');
  U1 = input('barrier height, eV: ');
  U0 = input('Eg value, eV: ');
  U1 = U1 + U0;
  E_steps = input('Energy discretisation, steps: ');
  E = linspace(U0 * 0.9,U1 * 1.1,E_steps); %widening energies specter to U0 - 0.1*U0 .. U1+0.1*U1
  x1 = x1 * x_step;
  x2 = x2 * x_step;
  x3 = x3 * x_step;
  x4 = x4 * x_step;
  A5_vals = A5_calc(x1, x2, x3, x4, U0, U1, E);
  figure;
  plot(E,A5_vals, 'b');
  hold on;
  title('transition coefficient vs electrone energy');
  xlabel('electrone energy, eV');
  ylabel('transition coefficient');
  grid on;
  hold off;
  clear all;
end

function A5 = A5_calc(x1, x2, x3, x4, U0, U1, E)
%consts
  nm2m = 1e-9;
  m0 = 0.911e-30;
  hs = 1.046e-34;
  ev = 1.6e-19;
%parameters conversion
  x1 = x1 * nm2m;
  x2 = x2 * nm2m;
  x3 = x3 * nm2m;
  x4 = x4 * nm2m;
  U0 = ev * U0;
  U1 = ev * U1;
  E = ev * E;
%calculation
  k1 = sqrt(2 * m0 * (E - U0) / hs^2);
  k2 = sqrt(2 * m0 * (U1 - E) / hs^2);
  res = 2 .* k2 .* exp(k2 .* x4) ./ ((k1 - 1i .* k1) .* exp(1i .* k1 .* x4)) .* 2 .* k1;
  denom = (exp(k2 .* x3 + 1i .* k1 .* x2) .* (k2 + 1i .* k1).^2 - (exp(k2 .* (2 .* x4 - x3) + 1i .* k1 .* x2) - exp(k2 .* (2 .* x4 - x3) + 1i .* k1 .* (2 .* x3 - x2)) + exp(k2 .* x3 + 1i .* k1 .* (x3 - x2))) .* (k2 - 1i .* k1).^2 .* (k1 + k2)) ./ (4i .* k1 .* k2 .* exp(k2 .* x2 .* 1i .* k1 .* x3));
  denom = denom +((k2 + 1i * k1).^2 .* (k2 - 1i .* k1) .* (exp(k2 .* (x3 - x4) + 1i .* k1 .* x2) + exp(k2 .* (x4 - x3) + 1i .* k1 .* (2 .* x3 - x2)) - exp(k2 .* x3 - 1i .* k1 .* (2* x3 - x2))) - exp(k2 .* (x4 - x3) + 1i .* k1 .* x2) .* (k2 - 1i .* k1).^3 .* (k1 - k2)) ./ (4i .* k1 .* k2 .* exp(-k2 .* (x4 + x2) + 1i .* k1 .* x3) .* (k2 + 1i .* k1));
  res = res ./ denom;
  A5 = abs(res).^2;
end