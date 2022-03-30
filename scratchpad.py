from pyswmm import Nodes, Simulation, Links
from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

sim = Simulation("./NO.inp")



integers = [1,2,3,4,5,6]
odds = [11,13,15,17,19,21]
evens = [20,22,24,26,28,30]
cols = ["integers","odds","evens"]

results = pd.DataFrame(np.transpose([integers,odds,evens]),columns=cols)

print(results["integers"][3])
print(results["evens"][4])
print("Yes")