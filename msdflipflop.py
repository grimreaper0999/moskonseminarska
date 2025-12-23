import numpy as np
import matplotlib.pyplot as plt
import grn, simulator
import networkx as nx
import faulthandler

faulthandler.enable()

mycell = grn.grn()

def registercell(name, network, clkname=None, inputname=None):

    # ------------------------
    # NODES
    # ------------------------

    # Initializes data and clock if input names not provided
    if not inputname:
        inputname = f"{name}_D"
        network.add_input_species(inputname)
    if not clkname:
        clkname = f"{name}_CLK"
        network.add_input_species(clkname)

    # Master latch internal nodes
    network.add_species(f"{name}_MNAND1", 0.1)
    network.add_species(f"{name}_MNAND2", 0.1)
    network.add_species(f"{name}_MNAND3", 0.1)
    network.add_species(f"{name}_MNAND4", 0.1)

    # Slave latch internal nodes
    network.add_species(f"{name}_SNAND1", 0.1)
    network.add_species(f"{name}_SNAND2", 0.1)

    # Outputs
    network.add_species(f"{name}_Q", 0.1)
    network.add_species(f"{name}_QBAR", 0.1)

    # ------------------------
    # MASTER LATCH (CLK = 1)
    # ------------------------

    # MNAND1 = NAND(D, CLK)
    mnand1_reg = [
        {'name': inputname,   'type': -1, 'Kd': 5, 'n': 3},
        {'name': clkname, 'type': -1, 'Kd': 5, 'n': 3}
    ]
    network.add_gene(10, mnand1_reg, [{"name": f"{name}_MNAND1"}], logic_type="or")

    # MNAND2 = NAND(~D, CLK)
    mnand2_reg = [
        {'name': inputname,   'type': 1,  'Kd': 5, 'n': 3},
        {'name': clkname, 'type': -1, 'Kd': 5, 'n': 3}
    ]
    network.add_gene(10, mnand2_reg, [{"name": f"{name}_MNAND2"}], logic_type="or")

    # MNAND3 = NAND(MNAND1, MNAND4)
    mnand3_reg = [
        {'name': f'{name}_MNAND1', 'type': -1, 'Kd': 5, 'n': 3},
        {'name': f'{name}_MNAND4', 'type': -1, 'Kd': 5, 'n': 3}
    ]
    network.add_gene(10, mnand3_reg, [{"name": f"{name}_MNAND3"}], logic_type="or")

    # MNAND4 = NAND(MNAND2, MNAND3)
    mnand4_reg = [
        {'name': f'{name}_MNAND2', 'type': -1, 'Kd': 5, 'n': 3},
        {'name': f'{name}_MNAND3', 'type': -1, 'Kd': 5, 'n': 3}
    ]
    network.add_gene(10, mnand4_reg, [{"name": f"{name}_MNAND4"}], logic_type="or")

    # ------------------------
    # SLAVE LATCH (CLK = 0)
    # ------------------------

    # SNAND1 = NAND(MNAND3, ~CLK)
    snand1_reg = [
        {'name': f'{name}_MNAND3', 'type': -1, 'Kd': 5, 'n': 3},
        {'name': clkname,    'type': 1,  'Kd': 5, 'n': 3}
    ]
    network.add_gene(10, snand1_reg, [{"name": f"{name}_SNAND1"}], logic_type="or")

    # SNAND2 = NAND(MNAND4, ~CLK)
    snand2_reg = [
        {'name': f'{name}_MNAND4', 'type': -1, 'Kd': 5, 'n': 3},
        {'name': clkname,    'type': 1,  'Kd': 5, 'n': 3}
    ]
    network.add_gene(10, snand2_reg, [{"name": f"{name}_SNAND2"}], logic_type="or")

    # ------------------------
    # OUTPUT LATCH
    # ------------------------

    # Q = NAND(SNAND1, QBAR)
    q_reg = [
        {'name': f'{name}_SNAND1', 'type': -1, 'Kd': 5, 'n': 3},
        {'name': f'{name}_QBAR',   'type': -1, 'Kd': 5, 'n': 3}
    ]
    network.add_gene(10, q_reg, [{"name": f"{name}_Q"}], logic_type="or")

    # QBAR = NAND(SNAND2, Q)
    qbar_reg = [
        {'name': f'{name}_SNAND2', 'type': -1, 'Kd': 5, 'n': 3},
        {'name': f'{name}_Q',      'type': -1, 'Kd': 5, 'n': 3}
    ]
    network.add_gene(10, qbar_reg, [{"name": f"{name}_QBAR"}], logic_type="or")

# ------------------------
# DEBUG / VISUALIZATION
# ------------------------

registercell("mycell", mycell)

print(mycell.genes)
mycell.plot_network()
plt.draw()

# ------------------------
# SIMULATION
# ------------------------

#preverjanje: nalaganje, hranjenje za enko in ničlo
#nalaganje - ob spremembi clock
#preveri hranjenje - če ni clock
#(negiran) izhod na vhod - preveri delovanje

clks = [0 if i%2==0 else 100 for i in range(10)]
data = [0 if i<5 else 100 for i in range(10)]
T, Y = simulator.simulate_sequence(mycell, [(data[i], clks[i]) for i in range(len(clks))], t_single = 250)