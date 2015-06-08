# Demo how to call TTV fast from Julia
if !isdefined(:TTVFAST_INITIALIZED)
  include("TTVFast.jl")
  using TTVFast
end

# Integration parameters
inflag = 0
t_start=0.0
t_stop=4*365.2425
dt = 0.02

# initializing array of input parameters to TTVFast
# should make into nice function
nplanets=2;  # hardwired for now
  duration = t_stop-t_start
  p = Array(Cdouble,2+nplanets*7);
  p[1] = 4pi/365.2425^2; # G
  p[2] = 1.0; # Mstar
  num_events=0;
  for i in 1:nplanets
    p[2+7*(i-1)+1] = 0.0001*rand();   # Planet Mass
    p[2+7*(i-1)+2] = 3.0^i; # Period
    p[2+7*(i-1)+3] = rand();  # Eccentricity
    p[2+7*(i-1)+4] = pi/2;    # Inclination
    p[2+7*(i-1)+5] = 2pi*rand(); # Longitude of Ascending Node
    p[2+7*(i-1)+6] = 2pi*rand(); # Argument of Pericenter
    p[2+7*(i-1)+7] = 2pi*rand(); # Mean Anomaly 
    num_events += iceil( duration/p[2+7*(i-1)+2] +1); # /* large enough to fit all the transits calculated by the code*/
  end

ttvfast_input = ttvfast_inputs_type(p, t_start=t_start, t_stop=t_stop, dt=dt)

incl_rvs = true  # Make sure first rv_time is _after_ t_start+dt/2 or else won't get any RV outputs
if incl_rvs
  num_rvs = 100
  rv_times = linspace(t_start+0.501*dt,t_stop, num_rvs)
  ttvfast_output = ttvfast_outputs_type(num_events ,rv_times)
else
  ttvfast_output = ttvfast_outputs_type(num_events)
end

println(STDERR, "# About to call TTVFast")
tic();
ttvfast!(ttvfast_input,ttvfast_output)
toc();

println("# You can inspect transit times with calls like: ")
println("  get_event(ttvfast_output,1)")
println("transit[1] = ",get_event(ttvfast_output,1))
if incl_rvs
  println("  get_rv(ttvfast_output,1)")
  println("rv[1] = ", get_rv(ttvfast_output,1))
end

println("# When you're done looking at outputs, remember to run...")
println("  free_ttvfast_outputs(ttvfast_output) ")

