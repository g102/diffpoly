function dy = diffpoly(x, y, w, p)
% dy = DIFFPOLY(x, y, w, p)
%   Compute the derivative of y(x) by computing the derivative of the the
%   fitting polynomial of degree p on a stencil of width w
%   At the boundaries of the domain, the stencil has always the same size
%   but it is not centered
%
%   Usage:
%     dy = DIFFPOLY(x, y, w, p)
%   
%   Inputs:
%     x: points where the function y(x) has been computed
%     y: values of the function y(x) in the points specified in 'x'
%     w: width of the stencil used to compute the polynomial (must be odd)
%     p: degree of the interpolating polynomial
%   All the input values are mandatory, none are optional
%
%   Output:
%     dy: estimated derivative y'(x) computed in the points in 'x'
%   
%   Limitations:
%     the stencil is assumed to be symmetric around each point
%     x and y must be vectors, with x increasing
%     y has no NaNs
%
%  Written by:
%     Stefano Gambuzza (s.gambuzza AT soton DOT ac DOT uk)

% Algorithm:
%   for x0 in x:
%       define a neighbourhood of x0 of radius hw
%       if the distance of x0 from the boundary is less than hw:
%           make the neighb. not centered around x0 but keep width = w
%       compute the best fitting polynomial of degree p in the neighb.
%       estimate dydx(x0) as the derivative of the poynomial computed

%% check on inputs
narginchk(4, 4);

% errors
if ~isvector(x) || ~isvector(y)
	error('diffpoly only works on vectors');
end
if length(x) ~= length(y)
	error('x and y must have the same number of elements');
end
if ~mod(w, 2)
	error('Stencil width must be odd');
end

% warnings
if floor(w) ~= w
	w = floor(w);
	warning('Stencil width is not integer, rounded down to %d', w);
end
if floor(p) ~= p
	p = floor(p);
	warning('Polynomial degree is not integer, rounded down to %d', p);
end
if p <= 0
	dy = zeros(size(x));
	warning('Degree of the polynomial lower than 1, returning zero derivative');
	return
end
if any(diff(x) < 0)
	error('x must be an increasing vector of points');
end
if p > w
	error('P is higher than the number of points in the stencil');
elseif p == w
	warning('P is equal to the number of points in the stencil');
end
if ~iscolumn(y)
	bReturnRowVector = true;
else
	bReturnRowVector = false;
end

%% redefine polyder to only do what we want (faster than Matlab's original)
polyder = @(x) x(1:end-1) .* (length(x)-1:-1:1);

%% preliminary computations
% ensure x and y are column vectors
x = x(:);
y = y(:);
dy = zeros(size(y));

% ensure y has no NaNs
if any(isnan(y))
	warning('y has NaN values, will replace with linear interp.');
	y = fillmissing(y, 'linear');
end

% shorthand
Nx = length(x);
hw = (w-1)/2;

% build single large Vandermonde matrix
V = ones(Nx, p+1);
for jj = 2 : p+1
	% can skip jj = 1 as that's all ones
	V(:, jj) = x.^(jj-1);
end

%% compute dy
% since the interpolating polynomial is the same in the first 'hw' points,
% only compute it once, then reuse it for all points with the same stencil
bLeftBoundComputed = false;
bRightBoundComputed = false;
for jj = 1:length(y)
	if jj <= hw
		% close to the left boundary
		i1 = 1;
		i2 = 2*hw + 1;
		if ~bLeftBoundComputed
			dp = polyder(flip(transpose(V(i1:i2, :)\y(i1:i2))));
			bLeftBoundComputed = true;
		end
	elseif jj >= Nx - hw
		% close to the right boundary
		i2 = Nx;
		i1 = i2 - 2*hw + 1;
		if ~bRightBoundComputed
			dp = polyder(flip(transpose(V(i1:i2, :)\y(i1:i2))));
			bLeftBoundComputed = true;
		end
	else
		% in the middle
		i1 = jj - hw;
		i2 = jj + hw;
		dp = polyder(flip(transpose(V(i1:i2, :)\y(i1:i2))));
	end
	dy(jj) = polyval(dp, x(jj));
end

%% return row vector if input was given as row vector
if bReturnRowVector
	dy = dy.';
end
