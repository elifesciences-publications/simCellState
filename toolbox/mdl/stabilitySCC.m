function dstab = stabilitySCC(eigJ, varargin)
%  dstab = stabilitySCC(eigJ, varargin)
% This function measure the stability of a SCC based on the eigenvalues of
% the jacobian. The output is the absolute value of the eigenvalue with
% maximum module with sign. This means that:
% - if all eigenvalues have negative real part, it returns the module of the largest one.
% - if all are negative and one is positive, it returns the module of the positive one.
% - if there is one (or more) zero eigenvalue and all negative, it returns zero.
% - if there is one (or more) zero eigenvalue, one (or more) positive, it returns the largest positive one.
% Optionally a shift can be passed as input (to force the SCC being
% subcritical).

if nargin == 2
    shift = varargin{1};
    if isempty(varargin{1})
        shift = 0;
    end
else
    shift = 0;
end
% output
dstab = abs(max(real(eigJ))-shift);