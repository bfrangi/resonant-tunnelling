from qumpy.simulators import PotentialSimulator, TransmissionSimulator


# Define the potential


def double_asymmetric_barrier_potential(x, V1, V2):
    barrier_1 = {
        "start": 1,  # nm
        "end": 5,  # nm
    }
    barrier_2 = {
        "start": 9,  # nm
        "end": 12,  # nm
    }
    bias = -0.63/12*x
    if x < barrier_1["start"]:
        return 0.0
    elif barrier_1["start"] <= x < barrier_1["end"]:
        return V1 + bias
    elif barrier_2["start"] <= x < barrier_2["end"]:
        return V2 + bias
    return bias


V1 = 0.5  # eV
V2 = 1.13  # eV
def V(x): return double_asymmetric_barrier_potential(x, V1, V2)


# Simulate the potential and transmission

potential_sim = PotentialSimulator(V)
potential_sim.x_max = 14 # nm
transmission_sim = TransmissionSimulator(potential_simulator=potential_sim)

potential_sim.show_potential_plot()
transmission_sim.show_transmission_plot()