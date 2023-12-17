from qumpy.simulators import (
    AnalyticDoubleSymmetricPotentialTransmissionSimulator,
    PotentialSimulator,
    TransmissionSimulator,
)


################################### Using the analytic simulator ###################################


transmission_sim = AnalyticDoubleSymmetricPotentialTransmissionSimulator(
    V1 = 0.5,  # eV
    V2 = 1.13,  # eV
    )
transmission_sim.w_1 = 2.0  # nm
transmission_sim.w_2 = 4.0  # nm
transmission_sim.w_3 = 3.6  # nm
transmission_sim.w_4 = 3.0  # nm
transmission_sim.w_5 = 1.4  # nm

E_arr_analytic = transmission_sim.compute_energy_array()
T_arr_analytic = transmission_sim.compute_transmission_array()


################################# Using the propagation simulator ##################################


# Define the potential

def double_asymmetric_barrier_potential(x, V1, V2):
    barrier_1 = {
        "start": 2,  # nm
        "end": 6,  # nm
    }
    barrier_2 = {
        "start": 9.6,  # nm
        "end": 12.6,  # nm
    }
    if barrier_1["start"] <= x < barrier_1["end"]:
        return V1
    elif barrier_2["start"] <= x < barrier_2["end"]:
        return V2
    return 0.0  # eV


V1 = 0.5  # eV
V2 = 1.13  # eV
def V(x): return double_asymmetric_barrier_potential(x, V1, V2)


# Simulate the potential and transmission

potential_sim = PotentialSimulator(V)
potential_sim.x_max = 14 # nm
transmission_sim = TransmissionSimulator(potential_simulator=potential_sim)

E_arr_propagation = transmission_sim.compute_energy_array()
T_arr_propagation = transmission_sim.compute_transmission_array()


############################################### Plot ###############################################


import matplotlib.pyplot as plt

plt.plot(E_arr_propagation, T_arr_propagation, "b-", label="Propagation")
plt.plot(E_arr_analytic, T_arr_analytic, "r--", label="Analytic")
plt.xlabel("E (eV)")
plt.ylabel("T")
plt.title("Transmission coefficient as a function of the electron energy")
plt.legend()

plt.show()