from qumpy.simulators import PotentialSimulator, TransmissionSimulator


# Define the potential


def single_barrier_potential(x, V0):
    barrier = {
        "start": 1.0,  # nm
        "end": 1.4,  # nm
    }
    if (barrier["start"] <= x < barrier["end"]):
        return V0
    return 0.0  # eV


V0 = 1.0  # eV
def V(x): return single_barrier_potential(x, V0)


# Simulate the potential and transmission

potential_sim = PotentialSimulator(V)
transmission_sim = TransmissionSimulator(potential_simulator=potential_sim)

potential_sim.show_potential_plot()
transmission_sim.show_transmission_plot()