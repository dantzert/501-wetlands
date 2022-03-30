import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


FWS_scen_A = pd.read_csv("FWS_A.csv")

fig, ax = plt.subplots(4,1)


ax[0].plot(FWS_scen_A["upstream_basin_depth_m"],'r')
ax[0].plot(FWS_scen_A["wetland_depth_m"],'b')
ax[0].legend(["upstream","wetland"])
## ax[0].ylabel("depth (m)")
## ax[0].xlabel("time")
#ax[0].title("FWS A depths")

ax[1].plot(FWS_scen_A["upstream_basin_valve_record"],"r")
ax[1].plot(FWS_scen_A["wetland_outlet_valve_record"],"b")
ax[1].legend(["upstream","wetland"])
##ax[1].ylabel("valve position")
## ax[1].xlabel("time")
#ax[1].title("FWS A Valve Actions")

ax[2].plot(FWS_scen_A["downstream_basin_cuload"],'k')
ax[2].plot(FWS_scen_A["wetland_cuload"],"b")
ax[2].plot(FWS_scen_A["outlet_channel_cuload"],'g')
##ax[2].ylabel("cumulative NO load")
##ax[2].title("FWS A Nitrate Removal")




plt.show()


print("yuh")