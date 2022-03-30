from pyswmm import Nodes, Simulation, Links
from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# adapted from : https://github.com/bemason/StormReactor-CaseStudy-Nitrate/blob/main/toolbox_NO.py

def CSTR_tank(t,C,Qin,Cin,Qout,V,k):
    dCdt = (Qin*Cin - Qout*C) / V - k*C
    return dCdt


def analysis(control_type,scenario=None):

    ## lists for result storage ##
    upstream_basin_inflow = []
    upstream_basin_NO = []
    upstream_basin_depth = []
    upstream_basin_outflow = []
    upstream_basin_cuload = []

    downstream_basin_inflow = []
    downstream_basin_NO = []
    downstream_basin_depth = []
    downstream_basin_outflow = []
    downstream_basin_cuload = []

    wetland_inflow = []
    wetland_NO = []
    wetland_depth = []
    wetland_volume = []
    wetland_outflow = []
    wetland_cuload = []
    wetland_DO = []  

    outlet_channel_flow = []
    outlet_channel_NO = []
    outlet_channel_depth = []
    outlet_channel_cuload = []

    wetland_outlet_valve_record = []
    upstream_basin_valve_record = []

    ##

    with Simulation("./No.inp") as sim:

        # create asset objects
        upstream_basin = Nodes(sim)["93-50408"]
        upstream_basin_valve = Links(sim)["95-70951"]
        upstream_control_box = Links(sim)["95-68287"]
        downstream_basin = Nodes(sim)["93-50404"]
        wetland = Nodes(sim)["93-49759"]
        wetland_bypass = Links(sim)["95-70294"]
        wetland_outlet_valve = Links(sim)["95-70293"]
        outlet_channel = Links(sim)["95-70277"]
   
        solver = [] # CSTR solvers
        for i in [0,1,2,3,4,5]: #0-2 are for DO, 3-5 are for NO
            solver.append(ode(CSTR_tank))
            solver[i].set_integrator("dopri5")
            solver[i].set_initial_value(0.0,0.0)

        # for iteration 0 on dt calcs
        start_time = sim.start_time
        last_timestep = start_time

        # tracking time to take control actions every 15 minutes (for a 5 second timestep)
        _tempcount = 180 

        for index,step in enumerate(sim):

            # calculate dt
            dt = (sim.current_time - last_timestep).total_seconds()
            last_timestep = sim.current_time

            # record NO concentrations
            upstream_basin_NO.append(upstream_basin.pollut_quality["NO"])
            downstream_basin_NO.append(downstream_basin.pollut_quality["NO"])
            wetland_NO.append(wetland.pollut_quality["NO"])
            outlet_channel_NO.append(outlet_channel.pollut_quality["NO"])

            # record quantity measures (flows and volumes)
            upstream_basin_inflow.append(upstream_basin.total_inflow)
            upstream_basin_depth.append(upstream_basin.depth)
            upstream_basin_outflow.append(upstream_basin.total_outflow)

            downstream_basin_inflow.append(downstream_basin.total_inflow)
            downstream_basin_depth.append(downstream_basin.depth)
            downstream_basin_outflow.append(downstream_basin.total_outflow)

            wetland_inflow.append(wetland.total_inflow)
            wetland_depth.append(wetland.depth)
            wetland_outflow.append(wetland.total_outflow)
            wetland_volume.append(wetland.volume)

            outlet_channel_flow.append(outlet_channel.flow)
            outlet_channel_depth.append(outlet_channel.depth)

            # record valve positions
            wetland_outlet_valve_record.append(wetland_outlet_valve.current_setting)
            upstream_basin_valve_record.append(upstream_basin_valve.current_setting)

            # the wetland is modeled as three CSTR's in series
            # DO transformation modeling
            # there is a cascade between the downstream_basin and the wetland so assume water enters well oxygenated
            k_DO = 0.000278  # rate/5 sec
            Cin_DO = 9.6     # mg/L

            for i in [0,1,2]: # DO ode instances
                if i ==0: # first section of wetland
                    solver[i].set_f_params(wetland_inflow[-1], Cin_DO, wetland_outflow[-1], wetland_volume[-1], k_DO)
                else: # middle and end
                    solver[i].set_f_params(wetland_inflow[-1], solver[i-1].y, wetland_outflow[-1], wetland_volume[-1], k_DO)
                solver[i].integrate(solver[i].t+dt)

            # this is the effluent DO concentration (conc in last "third" of wetland)
            wetland_DO.append(solver[2].y) 

            # assume because of minimal development that the wetland's subcatchment doesn't contribute nitrate
            # then influent nitrate concentration is that in the downstream_basin pool that flows into the wetland
            Cin_NO = downstream_basin_NO[-1] # most recent value

            for i in [3,4,5]: # NO ode instances
                if solver[i-3].y <= 1.0: # are anoxic conditions present in this third of the wetland?
                    k_ni = 0.000087
                else:
                    k_ni = 0.0

                if i ==3: # first section of wetland
                    solver[i].set_f_params(wetland_inflow[-1], Cin_NO, wetland_outflow[-1], wetland_volume[-1], k_ni)
                else: # middle and end
                    solver[i].set_f_params(wetland_inflow[-1], solver[i-1].y, wetland_outflow[-1], wetland_volume[-1], k_ni)
                solver[i].integrate(solver[i].t+dt)

            # this is the effluent NO concentration (conc in last "third" of wetland)
            wetland.pollut_quality["NO"] = solver[5].y
            # the internal (first and middle third) NO values will be tracked by the ode solvers
            # this is the value that determines loading into the channel so it makes sense to define it as the outlet


            ## control dynamics

            if (control_type == "uncontrolled") and (index == 0):
                print("\n\nno Real Time Control, passive operation")


            ## QUALITY ##
            elif (control_type == "nitrate_removal" and _tempcount==180): # rule-based control for water quality as specified in StormReactor Paper
                if index == 0:
                    print("\n\nControlling for Removal of Nitrate (Water Quality)\n")

                # Wetland & Retention basin Control Actions (every 15 mins - 5 sec timesteps)

                # If DO level is not anoxic
                if wetland_DO[-1] > 1.0:
                    # And if the wetland has capacity
                    if wetland_depth[-1] <= 9.5:
                        # Close the wetland valve and proportionally open retention basin valve C = Qmax/(A*sqrt(2*g*d))
                        wetland_outlet_valve.target_setting = 0.0
                        upstream_basin_valve.target_setting = 1.75*(70.6/(np.sqrt(2*32.2*upstream_basin_depth[-1])*78.5))
                    else:
                        # If not, open the wetland valve and close the RBasin valve
                        wetland_outlet_valve.target_setting = 1.75*(70.6/(np.sqrt(2*32.2*wetland_depth[-1])*12.6))
                        upstream_basin_valve.target_setting = 0.0
                # If DO level is anoxic
                elif wetland_DO[-1] <= 1.0:
                    # And if the (last third) wetland NO concentration is low, open both valves proportionally
                    if solver[5].y[0] <= 5.0:
                        wetland_outlet_valve.target_setting = 1.75*(70.6/(np.sqrt(2*32.2*wetland_depth[-1])*12.6))
                        upstream_basin_valve.target_setting = 1.75*(70.6/(np.sqrt(2*32.2*upstream_basin_depth[-1])*78.5))
                    # Else if the wetland NO concentration is high
                    else:
                        # And if the wetland still has capacity, close both valves
                        if wetland_depth[-1] <= 9.5:
                            wetland_outlet_valve.target_setting = 0.0
                            upstream_basin_valve.target_setting = 0.0
                        # If not, open the wetland valve propotionally and close retention basin
                        else:
                            wetland_outlet_valve.target_setting = 1.75*(70.6/(np.sqrt(2*32.2*wetland_depth[-1])*12.6))
                            upstream_basin_valve.target_setting = 0.0

                    _tempcount= 0 # reset timer


            ## QUANTITY ##
            # rule-based control to cap max discharge and reduce localized flooding (does not consider water quality)
            # this controller tries to create a set point discharge from the asset
            elif (control_type == "quantity" and _tempcount==180): 
                if index == 0:
                    print("\n\nControlling for Flooding and Discharge (Water Quantity)\n")
                    upstream_basin_valve.target_setting = 0 # start closed
                    wetland_outlet_valve.target_setting = 0 # start closed

                # Q_max for water quality was 70.6, take a tighter bound on this, say 50
                Q_desired = 50 # flow (cfs)
                
                # proportional control for flow rate
                K_p = 0.05

                # release if about to flood
                if upstream_basin_depth[-1] > 9.5:
                    upstream_basin_valve.target_setting = 1
                    _tempcount = -180*7 # leave it open for 2 hours to avoid oscillations
                else: # proportional control
                    upstream_error = Q_desired - upstream_control_box.flow # cfs
                    upstream_basin_valve.target_setting = upstream_basin_valve.current_setting + K_p*(upstream_error/Q_desired)
                    _tempcount= 0 # reset timer

                if wetland_depth[-1] > 9.5:
                    wetland_outlet_valve.target_setting = 1
                    _tempcount = -180*7 # leave it open for 2 hours to avoid oscillations
                else:
                    wetland_error = Q_desired - outlet_channel.flow
                    wetland_outlet_valve.target_setting = wetland_outlet_valve.current_setting + K_p*(wetland_error/Q_desired)
                    _tempcount= 0 # reset timer


                
              
            ## ECOLOGY ##
            ## Fish and Wildlife Waterfowl Management Scenarios A, B, and C
            # ref: https://www.nrcs.usda.gov/Internet/FSE_DOCUMENTS/nrcs142p2_016986.pdf page 27
            elif(control_type == "FWS" and _tempcount==180):
                    if index == 0:
                        print("\n\nControlling for Waterfowl Management (USFWS)\n")
                        print("Scenario: ", scenario, "\n")
                        upstream_basin_valve.target_setting = 0 # start closed
                        wetland_outlet_valve.target_setting = 0 # start closed

                    K_p = 0.01 # tune to correspond to 1-2 inches / week flooding 
                    # this will never quite be true, as the response to proportional control will be exponential, not linear

                    permanent_depth = 6 # the portion of the wetland depth we'll never drain

                    # Scenarios (Table 5)
                    if (scenario == "A"):
                        if(sim.current_time.month == 9): # early fall = september
                            desired_depth = 0 # "Dry"
                        elif(sim.current_time.month == 10): #mid fall = october
                            desired_depth = 0 # "Dry"
                        elif(sim.current_time.month == 11): # late fall = november
                            desired_depth = 1.25 # see reference for details
                            K_p = 2*K_p # 2-4 inches / week
                        elif(sim.current_time.month == 12 or sim.current_time.month == 1): # winter (but not late winter)
                            desired_depth = 2.5 # just below full capacity
                        elif(sim.current_time.month >= 2 and sim.current_time.month <= 5): # late winter through spring
                            K_p = K_p/10 # "slow drawdown"
                            desired_depth=0
                        else:
                            desired_depth=0
                            K_p = K_p/2 # allow for natural fluctuations, try to make fluctuations slow
                    elif(scenario == "B"):
                        if(sim.current_time.month == 9): # early fall = september
                            desired_depth = 0 # dry
                        elif(sim.current_time.month == 10): #mid fall = october
                            desired_depth = (1.5*4)/12 # see reference for details
                            # default K_p
                        elif(sim.current_time.month == 11): # late fall = november
                            desired_depth = 2
                        elif(sim.current_time.month == 12 or sim.current_time.month == 1): # winter (but not late winter)
                            desired_depth = 2.5
                        elif(sim.current_time.month >= 2 and sim.current_time.month <= 4): # late winter through spring
                            K_p = K_p/5
                            desired_depth= 0
                        else:
                            desired_depth= 0
                            K_p = K_p/2 # allow for natural fluctuations, try to make fluctuations slow
                    elif(scenario == "C"):
                        if(sim.current_time.month == 9): # early fall = september
                            desired_depth = 4/12
                            K_p = K_p / 2 # "gradual flooding"
                        elif(sim.current_time.month == 10): #mid fall = october
                            desired_depth = 1.5
                            K_p = K_p / 2
                        elif(sim.current_time.month == 11): # late fall = november
                            desired_depth = 2.5 # full functional capacity (9.5 feet / 10 feet)
                        elif(sim.current_time.month == 12 or sim.current_time.month == 1): # winter (but not late winter)
                            desired_depth = 3 # full pool
                        elif(sim.current_time.month >= 2 and sim.current_time.month <= 4): # late winter through spring
                            K_p = K_p/5
                            desired_depth= 0
                        else:
                            desired_depth= 0
                            K_p = K_p/2 # allow for natural fluctuations, try to make fluctuations slow


                    # control logic is just a proportional feedback, but on depth this time
                    error =  wetland_depth[-1] - (permanent_depth+desired_depth)
                    wetland_outlet_valve.target_setting = wetland_outlet_valve.current_setting + K_p*(error/wetland_depth[-1]) # adjust the wetland valve
                    # if error < 0 this will close the wetland valve and retain more water
                    # if error > 0 this will open the wetland valve and release more water
                    _tempcount= 0 # reset timer
                    if (error < 0.25): # more than four inches below our desired setpoint
                        upstream_basin_valve.target_setting = upstream_basin_valve.current_setting - K_p*(error/wetland_depth[-1]) # open the upstream valve
                        # we can get rid of water with just local control, but we need the upstream asset if we want to add water

            





            _tempcount+= 1


        sim._model.swmm_end()
        print("surface runoff error = ",sim.runoff_error, " %")
        print("flow routing error = ",sim.flow_routing_error, " %")
        print("quality routing error = ",sim.quality_error, " %")

    # Convert inflow rate from cfs to m3/s
    conv_cfs_cms = [0.02832]*len(upstream_basin_inflow)
    upstream_basin_inflow_m = [a*b for a,b in zip(upstream_basin_inflow,conv_cfs_cms)]
    downstream_basin_inflow_m = [a*b for a,b in zip(downstream_basin_inflow,conv_cfs_cms)]
    wetland_inflow_m = [a*b for a,b in zip(wetland_inflow,conv_cfs_cms)]
    outlet_channel_flow_m = [a*b for a,b in zip(outlet_channel_flow,conv_cfs_cms)]

    # Convert outflow rate from cfs to m3/s
    conv_cfs_cms = [0.02832]*len(upstream_basin_inflow)
    upstream_basin_outflow_m = [a*b for a,b in zip(upstream_basin_outflow,conv_cfs_cms)]
    downstream_basin_outflow_m = [a*b for a,b in zip(downstream_basin_outflow,conv_cfs_cms)]
    wetland_outflow_m = [a*b for a,b in zip(wetland_outflow,conv_cfs_cms)]

    # Convert depth from ft to m
    conv_ft_m = [0.3048]*len(upstream_basin_inflow)
    upstream_basin_depth_m = [a*b for a,b in zip(upstream_basin_depth,conv_ft_m)]
    downstream_basin_depth_m = [a*b for a,b in zip(downstream_basin_depth,conv_ft_m)]
    wetland_depth_m = [a*b for a,b in zip(wetland_depth,conv_ft_m)]
    outlet_channel_depth_m = [a*b for a,b in zip(outlet_channel_depth,conv_ft_m)]

    # Calculate load each timestep
    conv_mgs_kgs = [0.000001]*len(upstream_basin_inflow)
    timestep = [5]*len(upstream_basin_inflow)
    upstream_basin_load = [a*b*c*d*e for a,b,c,d,e in zip(upstream_basin_NO,upstream_basin_outflow,conv_cfs_cms, conv_mgs_kgs,timestep)]
    downstream_basin_load = [a*b*c*d*e for a,b,c,d,e in zip(downstream_basin_NO,downstream_basin_outflow,conv_cfs_cms,conv_mgs_kgs,timestep)]
    wetland_load = [a*b*c*d*e for a,b,c,d,e in zip(wetland_NO,wetland_outflow,conv_cfs_cms, conv_mgs_kgs,timestep)]
    outlet_channel_load = [a*b*c*d*e for a,b,c,d,e in zip(outlet_channel_NO,outlet_channel_flow,conv_cfs_cms,conv_mgs_kgs,timestep)]

    # Calculate cumulative load (dt = 1)
    upstream_basin_cuload = np.cumsum(upstream_basin_load)
    downstream_basin_cuload = np.cumsum(downstream_basin_load)
    wetland_cuload = np.cumsum(wetland_load)
    outlet_channel_cuload = np.cumsum(outlet_channel_load)
    cols = ["upstream_basin_inflow_m", "downstream_basin_inflow_m","wetland_inflow_m", "outlet_channel_flow_m",
                            "upstream_basin_outflow_m", "downstream_basin_outflow_m", "wetland_outflow_m", 
                            "upstream_basin_depth_m", "downstream_basin_depth_m", "wetland_depth_m", "outlet_channel_depth_m",
                            "upstream_basin_cuload", "downstream_basin_cuload", "wetland_cuload", "outlet_channel_cuload",
                            "wetland_outlet_valve_record","upstream_basin_valve_record"]
    results = pd.DataFrame(np.transpose([upstream_basin_inflow_m, downstream_basin_inflow_m,wetland_inflow_m, outlet_channel_flow_m,
                            upstream_basin_outflow_m, downstream_basin_outflow_m, wetland_outflow_m, 
                            upstream_basin_depth_m, downstream_basin_depth_m, wetland_depth_m, outlet_channel_depth_m,
                            upstream_basin_cuload, downstream_basin_cuload, wetland_cuload, outlet_channel_cuload,
                            wetland_outlet_valve_record, upstream_basin_valve_record]),
                           columns=cols)
    return results



uncontrolled_scenario = analysis("uncontrolled")
uncontrolled_scenario.to_csv("uncontrolled.csv")

FWS_A = analysis("FWS","A")
FWS_A.to_csv("FWS_A.csv")

FWS_B = analysis("FWS","B")
FWS_B.to_csv("FWS_B.csv")

FWS_C = analysis("FWS","C")
FWS_C.to_csv("FWS_C.csv")
    
quantity_control = analysis("quantity")
quantity_control.to_csv("quantity_control.csv")

nitrate_removal_control = analysis("nitrate_removal")
nitrate_removal_control.to_csv("nitrate_removal.csv")



