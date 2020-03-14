function [time, Q, Qa] = extProbTime(rgFnc, lambda, X0, tF)
% This function estimate the extinction probability as function of time.
%  rgFnc IN: reproductive generating function for each state (N,1)
%  lambda IN: hazard rate (1/expected lifespan) (N,1)
%  X0    IN: initial condition (N,1)
%  tF    IN: final time
%  time OUT: time (1,M)
%  Q    OUT: extinction probability as function of time (N,M)
%  Qa   OUT: asymptotic extinction probability

% initialization
S0 = double(X0==0);
% derivative of ft
dftdt = @(t, ft) lambda.*(rgFnc(ft) - ft);

% numerical integration
[time,ft] = ode45(dftdt, [0 tF], S0);
time = time'; ft = ft';
Q = ft.^repmat(X0, 1, length(time));

% asymptotic result
tol = 1e-10;
[DFDY,FAC] = numjac(@(t, x) rgFnc(x), 0, S0, rgFnc(S0) , tol, [], 0);
if any(abs(abs(eig(DFDY))-1) < tol) % unitary eigenvalue in the map
    Qa = prod(S0.^X0); % probably here it should be the minimum of Q along the orbit (periodic)
else
    dfun = @(ft) (rgFnc(ft) - ft);
    Qa = prod(fsolve(dfun,Q(:,end)));
end

% [DFDY,FAC] = numjac(@(t, x) rgFnc(x), 0, ones(size(S0)), rgFnc(ones(size(S0))) , tol, [], 0);
% eig(DFDY)

return