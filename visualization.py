import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



def viz(scenario_name):

    pyswmm_output = pd.read_csv(str(scenario_name + ".csv"))

    with open(str(scenario_name + "_NO_results.txt"),'w') as f:
        f.write("NO removal efficiency:\n")
        f.write("Cumulative Input = ")
        f.write(str(pyswmm_output["wetland_cuload"].values[-1:]))
        f.write("\nCumulative Output = ")
        f.write( str(pyswmm_output["wetland_effluent_cuload"].values[-1:]))
        f.write("\n1 - (Input / Output) = ")
        f.write(str(1-pyswmm_output["wetland_effluent_cuload"].values[-1:]/pyswmm_output["wetland_cuload"].values[-1:]))


    fig, ax = plt.subplots(4,1)
    fig.set_size_inches(7,4)
    fig.dpi = 500 # high res plot
    plt.tight_layout()
    

    
    ax[0].plot(pyswmm_output["upstream_basin_depth_m"],'r')
    ax[0].plot(pyswmm_output["wetland_depth_m"],'b')
    #ax[0].legend(["upstream depth (m)","wetland depth (m)"])
    ax[0].legend(["upstream", "wetland"],loc="upper left",fontsize="x-small")
    ax[0].set_title("depth (m)",fontsize="x-small")


    ax[1].plot(pyswmm_output["upstream_basin_valve_record"],"r")
    ax[1].plot(pyswmm_output["wetland_outlet_valve_record"],"b")
    #ax[1].legend(["upstream valve setting","wetland valve setting"])
    ax[1].legend(["upstream", "wetland"],loc="upper left",fontsize="x-small")
    ax[1].set_title("valve open %",fontsize="x-small")


    ax[2].plot(pyswmm_output["wetland_cuload"],"b")
    ax[2].plot(pyswmm_output["wetland_effluent_cuload"],'g')
    #ax[2].legend(["cumulative load (NO) into wetland","cumulative load (NO) out of wetland"])
    ax[2].legend(["inflow", "outflow"],loc="upper left",fontsize="x-small")
    ax[2].set_title("cumulative load (kg NO) - wetland",fontsize="x-small")

    ax[3].plot(pyswmm_output["wetland_inflow_m"],"b")
    ax[3].plot(pyswmm_output["wetland_outflow_m"],"g")
    #ax[3].legend(["wetland inflow (cms)","wetland outflow (cms)"])
    ax[3].legend(["inflow", "outflow"],loc="upper left",fontsize="x-small")
    ax[3].set_title("flow (cms) - wetland",fontsize="x-small")
    
    

    plt.savefig(str(scenario_name + "_results.jpg"),dpi=500)

    print("\nfig done\n")


viz("uncontrolled")
viz("quantity_control")
viz("nitrate_removal")
viz("FWS_A")
viz("FWS_B")
viz("FWS_C")