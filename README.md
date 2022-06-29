# diffpoly-matlab

Compute the derivative of y(x) by computing the derivative of the the fitting polynomial of degree POLYDEG on a stencil of half-width w (total width 2*w + 1).
At the boundaries of the domain, the stencil has always the same size but it is not centered.

## Usage
`dy = diffpoly(x, y, w, p)`
 
 ## Inputs
  * x: points where the function y(x) has been computed
  * y: values of the function y(x)
  * w: half-width of the stencil to be used to compute the approximating polynomial
  * p: degree of the interpolating polynomial
  
  All inputs are mandatory, none are optional

## Output
 * dy: estimate of dy/dx computed in the points in `x`

## Limitations
 * The stencil is always symmetric (from -w to w); this doesn't mean the points in the stencil need to be equally spaced
 * `x` and `y` must be vectors, with `x` increasing
 * `y` has no NaNs
 
## Algorithm
```
  for x0 in x:
      define a neighbourhood of x0 of radius w
      if the distance of x0 from the boundary is less than w:
          make the neighb. not centered around x0 but keep width = w
      compute the best fitting polynomial of degree p in the neighb.
      estimate dydx(x0) as the derivative of the poynomial computed
 ```
