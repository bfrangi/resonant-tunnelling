from qumpy.simulators import AnalyticDoubleSymmetricPotentialTransmissionSimulator

# Simulate the transmission

transmission_sim = AnalyticDoubleSymmetricPotentialTransmissionSimulator(
    V1 = 0.5,  # eV
    V2 = 1.13,  # eV
    )
transmission_sim.w_1 = 2.0  # nm
transmission_sim.w_2 = 4.0  # nm
transmission_sim.w_3 = 3.6  # nm
transmission_sim.w_4 = 3.0  # nm
transmission_sim.w_5 = 1.4  # nm

transmission_sim.show_transmission_plot()