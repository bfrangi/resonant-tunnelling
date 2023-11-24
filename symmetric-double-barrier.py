from qumpy.simulators import PotentialSimulator, TransmissionSimulator


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

potential_sim.show_potential_plot()
transmission_sim.show_transmission_plot()