module TTVFast

function __init__()
  const global LIBTTV = find_library(["libttvfast.so"],["/usr/local/lib",".","../c_version/"])
  const global ttvfast_NA_DEFAULT = convert(Cdouble,-2.0); # /* value for transit information that is not determined by TTVFast*/
  const global TTVFAST_INITIALIZED = true
end

export ttvfast_RV_entry, ttvfast_transit_entry, ttvfast_outputs_type, ttvfast_inputs_type
export ttvfast!, get_event, get_rv, get_rv_time, free_ttvfast_outputs
export ttvfast_NA_DEFAULT, TTVFAST_INITIALIZED

# Immutable Types to match structs from TTVFast
immutable ttvfast_RV_entry
  time::Cdouble
  RV::Cdouble
end

immutable ttvfast_transit_entry
  planet::Cint
  epoch::Cint
  time::Cdouble
  rsky::Cdouble
  vsky::Cdouble
end

# types to package together inputs and outputs to TTVFast
type ttvfast_outputs_type
  transit_workspace::Ptr{ttvfast_transit_entry}
  rv_workspace::Ptr{ttvfast_RV_entry}
  max_num_events::Cint
  num_rv_obs::Cint
end

type ttvfast_inputs_type
  param::Ptr{Cdouble}
  time_step::Cdouble
  t_start::Cdouble
  t_stop::Cdouble
  num_planets::Cint
  input_flag::Cint
end


function ttvfast_outputs_type(nevents::Integer, rv_times::Vector{Float64} = Array(Float64,0))
  @assert(nevents>=1);
  local nrvs = length(rv_times)
  local transit_workspace = c_calloc(nevents,sizeof(ttvfast_transit_entry)) 
  local rv_workspace = c_calloc(nrvs,sizeof(ttvfast_RV_entry))
  @assert(transit_workspace != C_NULL);
  @assert(rv_workspace != C_NULL);
  for i in 1:nrvs
    pointer_to_array(convert(Ptr{ttvfast_RV_entry},rv_workspace),nrvs)[i] = ttvfast_RV_entry(rv_times[i],0.)
  end
  return ttvfast_outputs_type(transit_workspace,rv_workspace,convert(Cint,nevents),convert(Cint,nrvs) );
end

function free_ttvfast_outputs(data::ttvfast_outputs_type)
  try
    c_free(data.rv_workspace);
    c_free(data.transit_workspace);
  catch
    println(STDERR, "# ERROR: Something funky happened while trying to deallocate the space reserved for TTVFast outputs.")
  end
  data.max_num_events = 0;
  data.num_rv_obs = 0;
end

function get_event(data::ttvfast_outputs_type, i::Integer)  
  @assert( (i>=1) && (i<=data.max_num_events) )
  # What not use data?
  pointer_to_array(data.transit_workspace,data.max_num_events)[i]
  #pointer_to_array(ttvfast_output.transit_workspace,data.max_num_events)[i]
end

function get_rv(data::ttvfast_outputs_type, i::Integer, convert_to_mps::Bool=true)  
  @assert( (i>=1) && (i<=data.num_rv_obs) )
  # What not use data?
  rv = pointer_to_array(data.rv_workspace,data.num_rv_obs)[i].RV
  #pointer_to_array(ttvfast_output.rv_workspace,data.num_rv_obs)[i].RV
  if(convert_to_mps)
    rv *= 1731456.84
  end
  return rv
end

function get_rv_time(data::ttvfast_outputs_type, i::Integer)  
  @assert( (i>=1) && (i<=data.num_rv_obs) )
  # What not use data?
  pointer_to_array(data.rv_workspace,data.num_rv_obs)[i].time
  #pointer_to_array(ttvfast_output.rv_workspace,data.num_rv_obs)[i].time
end

function ttvfast_inputs_type(p::Array{Cdouble,1} ; inflag::Integer = 0, t_start::Real=0., t_stop::Real=4*365.2425, dt::Real = 0.02)
  const local max_steps = 300000;
  local np = fld(length(p)-2,7);
  @assert(np>=2);
  @assert( length(p) == np*7+2 );
  @assert( (inflag==0) || (inflag==1) || (inflag==2) );
  @assert (t_stop>t_start);
  @assert( (t_stop-t_start) < dt*max_steps);
  return ttvfast_inputs_type(pointer(p),dt,t_start,t_stop,np,inflag);
end

function ttvfast!(in::ttvfast_inputs_type, out::ttvfast_outputs_type)
  # void TTVFast(double *params,double dt, double Time, double total,int n_planets, CalcTransit *transit, CalcRV *RV_struct,int n,int m, int input);
  ccall( (:TTVFast, LIBTTV ), Void, (Ptr{Void},Cdouble,Cdouble,Cdouble,Cint,Ptr{Void},Ptr{Void},Cint,Cint,Cint,),
                  in.param,in.time_step,in.t_start,in.t_stop,in.num_planets,
                  #out.transit_workspace,out.rv_workspace,out.num_rv_obs,out.max_num_events,in.input_flag )
                  convert(Ptr{Void},out.transit_workspace),convert(Ptr{Void},out.rv_workspace),out.num_rv_obs,out.max_num_events,in.input_flag )
end

end # Module TTVFast

