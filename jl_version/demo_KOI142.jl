# Demo how to call TTV fast from Julia
if !isdefined(:TTVFAST_INITIALIZED)
  include("TTVFast.jl")
  using TTVFast
end

path = "../c_version/"
setup_filename = "setup_file_astro_cartesian"
RV_in_filename = "RV_file"

# Read Integration parameters from file
setup_param = split(readall(string(path,setup_filename)))
inputfilename = string(path,setup_param[1])
t_start=parse(setup_param[2])
dt = parse(setup_param[3])
t_stop = parse(setup_param[4])
nplanets = parse(setup_param[5]) 
inflag = parse(setup_param[5]) 

# initializing array of input parameters to TTVFast using data from file
p = Array(Cdouble,2+nplanets*7);
param_str = split(readall(inputfilename))
@assert(length(p) == length(param_str))
{p[i] = parse(param_str[i]) for i in 1:length(p) }
num_events=5000  # hardwired same as TTVFast

ttvfast_input = ttvfast_inputs_type(p, t_start=t_start, t_stop=t_stop, dt=dt, inflag=inflag)

incl_rvs = true  # Make sure first rv_time is _after_ t_start+dt/2 or else won't get any RV outputs
if incl_rvs
  rvfile_str = split(readall(string(path,RV_in_filename)))
  rv_times = zeros(Float64,length(rvfile_str))
  {rv_times[i] = parse(rvfile_str[i]) for i in 1:length(rvfile_str) }
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



