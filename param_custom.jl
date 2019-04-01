

function add_param_custom( sim_param::SimParam )
   add_param_active(sim_param,"sigma_hk",0.03)
   #add_param_active(sim_param,"sigma_hk_one",0.3)
   #add_param_active(sim_param,"sigma_hk_multi",0.03)
   add_param_active(sim_param,"sigma_incl",1.0)   # degrees; 0 = coplanar w/ generate_kepler_target_simple; ignored by generate_planetary_system_uncorrelated_incl
   return sim_param
end
 
