from qumpy.simulators import AnalyticDoubleSymmetricPotentialTransmissionSimulator

# Simulate the transmission

transmission_sim = AnalyticDoubleSymmetricPotentialTransmissionSimulator(
    V1 = 0.3427625,  # eV
    V2 = 0.579025,  # eV
    )
transmission_sim.w_1 = 1.0  # nm
transmission_sim.w_2 = 4.0  # nm
transmission_sim.w_3 = 4.0  # nm
transmission_sim.w_4 = 3.0  # nm
transmission_sim.w_5 = 2.0  # nm

transmission_sim.show_transmission_plot()