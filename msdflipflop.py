import numpy as np
import matplotlib.pyplot as plt
import grn, simulator
import networkx as nx
import faulthandler

def registercell(name, network, clkname=None, inputname=None, decays=None, Kds=None, ns=None):

    # ------------------------
    # PARAMETERS
    # ------------------------

    # Initializes base parameters if not provided
    if not decays:
        decays=[0.1]
    if not Kds:
        Kds=[5]
    if not ns:
        ns=[3]

    # Checks for different patterns of initializing
    match len(decays):
        # One value for all
        case 1:
            decays = np.repeat(decays,8)

    match len(Kds):
        # One value for all
        case 1:
            Kds = np.repeat(Kds, 16)
        # One value for activators, one for repressors
        case 2:
            Kds = [Kds[0] if i in [2, 9, 11] else Kds[1]
                   for i in range(16)]

    match len(ns):
        # One value for all
        case 1:
            ns = np.repeat(ns, 16)
        # One value for activators, one for repressors
        case 2:
            ns = [ns[0] if i in [2, 9, 11] else ns[1]
                  for i in range(16)]

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
    network.add_species(f"{name}_MNAND1", decays[0])
    network.add_species(f"{name}_MNAND2", decays[1])
    network.add_species(f"{name}_MNAND3", decays[2])
    network.add_species(f"{name}_MNAND4", decays[3])

    # Slave latch internal nodes
    network.add_species(f"{name}_SNAND1", decays[4])
    network.add_species(f"{name}_SNAND2", decays[5])

    # Outputs
    network.add_species(f"{name}_Q", decays[6])
    network.add_species(f"{name}_QBAR", decays[7])

    # ------------------------
    # MASTER LATCH (CLK = 1)
    # ------------------------

    # MNAND1 = NAND(D, CLK)
    mnand1_reg = [
        {'name': inputname,   'type': -1, 'Kd': Kds[0], 'n': ns[0]},
        {'name': clkname,     'type': -1, 'Kd': Kds[1], 'n': ns[1]}
    ]
    network.add_gene(10, mnand1_reg, [{"name": f"{name}_MNAND1"}], logic_type="or")

    # MNAND2 = NAND(~D, CLK)
    mnand2_reg = [
        {'name': inputname,   'type':  1, 'Kd': Kds[2], 'n': ns[2]},
        {'name': clkname,     'type': -1, 'Kd': Kds[3], 'n': ns[3]}
    ]
    network.add_gene(10, mnand2_reg, [{"name": f"{name}_MNAND2"}], logic_type="or")

    # MNAND3 = NAND(MNAND1, MNAND4)
    mnand3_reg = [
        {'name': f'{name}_MNAND1', 'type': -1, 'Kd': Kds[4], 'n': ns[4]},
        {'name': f'{name}_MNAND4', 'type': -1, 'Kd': Kds[5], 'n': ns[5]}
    ]
    network.add_gene(10, mnand3_reg, [{"name": f"{name}_MNAND3"}], logic_type="or")

    # MNAND4 = NAND(MNAND2, MNAND3)
    mnand4_reg = [
        {'name': f'{name}_MNAND2', 'type': -1, 'Kd': Kds[6], 'n': ns[6]},
        {'name': f'{name}_MNAND3', 'type': -1, 'Kd': Kds[7], 'n': ns[7]}
    ]
    network.add_gene(10, mnand4_reg, [{"name": f"{name}_MNAND4"}], logic_type="or")

    # ------------------------
    # SLAVE LATCH (CLK = 0)
    # ------------------------

    # SNAND1 = NAND(MNAND3, ~CLK)
    snand1_reg = [
        {'name': f'{name}_MNAND3', 'type': -1, 'Kd': Kds[8], 'n': ns[8]},
        {'name': clkname,          'type':  1, 'Kd': Kds[9], 'n': ns[9]}
    ]
    network.add_gene(10, snand1_reg, [{"name": f"{name}_SNAND1"}], logic_type="or")

    # SNAND2 = NAND(MNAND4, ~CLK)
    snand2_reg = [
        {'name': f'{name}_MNAND4', 'type': -1, 'Kd': Kds[10], 'n': ns[10]},
        {'name': clkname,          'type':  1, 'Kd': Kds[11], 'n': ns[11]}
    ]
    network.add_gene(10, snand2_reg, [{"name": f"{name}_SNAND2"}], logic_type="or")

    # ------------------------
    # OUTPUT LATCH
    # ------------------------

    # Q = NAND(SNAND1, QBAR)
    q_reg = [
        {'name': f'{name}_SNAND1', 'type': -1, 'Kd': Kds[12], 'n': ns[12]},
        {'name': f'{name}_QBAR',   'type': -1, 'Kd': Kds[13], 'n': ns[13]}
    ]
    network.add_gene(10, q_reg, [{"name": f"{name}_Q"}], logic_type="or")

    # QBAR = NAND(SNAND2, Q)
    qbar_reg = [
        {'name': f'{name}_SNAND2', 'type': -1, 'Kd': Kds[14], 'n': ns[14]},
        {'name': f'{name}_Q',      'type': -1, 'Kd': Kds[15], 'n': ns[15]}
    ]
    network.add_gene(10, qbar_reg, [{"name": f"{name}_QBAR"}], logic_type="or")


def truthgenerator(clks):
    
    # Extract rising edges
    pts = [0] + [i for i in range(len(clks)-1) if clks[i+1]-clks[i]==100] + [len(clks)-1]
    truth = np.array([])

    # Construct ground truth based on rising edges
    for p in range(len(pts)-1):
        truth = np.concatenate([truth, np.repeat(100-(100 if len(truth)==0 else truth[-1]), pts[p+1]-pts[p])])
    return np.concatenate([clks[0], truth])



if __name__ == "__main__":

    faulthandler.enable()

    # ------------------------
    # DEBUG / VISUALIZATION
    # ------------------------

    mycell = grn.grn()

    registercell("mycell", mycell, inputname="mycell_QBAR")

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
    clks = [0 for i in range(10)]
    clks = [0 if i%2==0 and i < 5 else 100 for i in range(10)]
    clks = np.concat([np.repeat(0, 5),[100 if i%4==0 else 0 for i in range(10)]])
    data = [0 if i<5 else 100 for i in range(10)]
    data = [100 for i in range(10)]
    T, Y = simulator.simulate_sequence(mycell, clks, t_single = 250, plot_on=False) #[(data[i], clks[i]) for i in range(len(clks))]
    t = truthgenerator(Y[:, [0]])
    print(t)
    plt.plot(T, Y[:, [0]], '--')
    plt.plot(T, Y[:, [7, 8]])
    plt.plot(T, t)
    plt.legend(np.array(mycell.species_names)[[0, 7, 8]])
    plt.show()