from qumpy.simulators import PotentialSimulator, TransmissionSimulator


# Define the potential


def parabolic_potential(x, V0):
    return max(0.0, V0 - 0.5*V0*(x-1.0)**2)


V0 = 1.0  # eV
def V(x): return parabolic_potential(x, V0)


# Simulate the potential and transmission

potential_sim = PotentialSimulator(V)
transmission_sim = TransmissionSimulator(potential_simulator=potential_sim)

potential_sim.show_potential_plot()
transmission_sim.show_transmission_plot()