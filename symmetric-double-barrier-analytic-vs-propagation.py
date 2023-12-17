from qumpy.simulators import (
    AnalyticDoubleSymmetricPotentialTransmissionSimulator,
    PotentialSimulator,
    TransmissionSimulator,
)


################################### Using the analytic simulator ###################################


transmission_sim = AnalyticDoubleSymmetricPotentialTransmissionSimulator()
E_arr_analytic = transmission_sim.compute_energy_array()
T_arr_analytic = transmission_sim.compute_transmission_array()


################################# Using the propagation simulator ##################################


# Define the potential

def double_symmetric_barrier_potential(x, V0):
    barrier_1 = {
        "start": 0.2,  # nm
        "end": 0.6,  # nm
    }
    barrier_2 = {
        "start": 1.2,  # nm
        "end": 1.6,  # nm
    }
    if (barrier_1["start"] <= x < barrier_1["end"]) or (
            barrier_2["start"] <= x < barrier_2["end"]):
        return V0
    return 0.0  # eV


V0 = 1.0  # eV
def V(x): return double_symmetric_barrier_potential(x, V0)


# Simulate the potential and transmission

potential_sim = PotentialSimulator(V)
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