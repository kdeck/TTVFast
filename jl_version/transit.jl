# here are the definitions of the CalcTransit and CalcRV structure.
type CalcTransit
  planet::Integer 
  epoch::Integer
  time::Real
  rsky::Real
  vsky::Real
end

type CalRV
  time::Real
  RV::Real
end
