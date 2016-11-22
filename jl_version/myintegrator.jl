//Here are some other parameters used.

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define ABS(a) ((a) < 0 ? -(a) : (a))

type Vector
  x::Real
  y::Real
  z::Real
end

type PhaseState
  x::Real
  y::Real
  z::Real
  xd::Real
  yd::Real
  zd::Real
end
