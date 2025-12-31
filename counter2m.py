import numpy as np
import matplotlib.pyplot as plt
import grn, simulator
import networkx as nx
import faulthandler
from msdflipflop import registercell

def counterregister(network, d, clkname=None, cell_decays=None, cell_Kds=None, cell_ns=None, instr_decays=None, conn_Kds=None, conn_ns=None):

    # ------------------------
    # PARAMETERS
    # ------------------------

    # Initializes base parameters if not provided
    if not instr_decays:
        instr_decays=[0.1]
    if not conn_Kds:
        conn_Kds=[5]
    if not conn_ns:
        conn_ns=[3]

    # Checks for different patterns of initializing
    match len(instr_decays):
        # One value for all
        case 1:
            instr_decays = np.repeat(instr_decays, 2*d)

    match len(conn_Kds):
        # One value for all
        case 1:
            conn_Kds = np.repeat(conn_Kds, 4*d)
        # One value for activators, one for repressors
        case 2:
            conn_Kds = [conn_Kds[0] if i in [0, 1]+[(2*i-1) for i in range(2, d+1)]+[(2*d+2*i+1) for i in range(2, d+1)] else conn_Kds[1]
                   for i in range(4*d)]

    match len(conn_ns):
        # One value for all
        case 1:
            conn_ns = np.repeat(conn_ns, 4*d)
        # One value for activators, one for repressors
        case 2:
            conn_ns = [conn_ns[0] if i in [0, 1]+[(2*i-1) for i in range(2, d+1)]+[(2*d+2*i+1) for i in range(2, d+1)] else conn_ns[1]
                   for i in range(4*d)]


    # ------------------------
    # BUILDING THE COUNTER
    # ------------------------

    # Initializing the clock input if not given
    if not clkname:
        clkname="CLK"
        network.add_input_species(clkname)

    # Initialize the first cell with clk the negated output of the last 
    registercell("CELL_1", network, clkname, f"CELL_{d}_QBAR", cell_decays, cell_Kds, cell_ns)

    # Initialize the rest of the cells with clk and output of the previous cell
    for c in range(2, d+1):
        registercell(f"CELL_{c}", network, clkname, f"CELL_{c-1}_Q", cell_decays, cell_Kds, cell_ns)

    # ------------------------
    # BUILDING THE DECODER
    # ------------------------

    # Initialize the numbered instruction outputs
    for i in range(1, 2*d+1):
        network.add_species(f"INSTRUCTION_{i}", instr_decays[i-1])

    # Instruction 1 is triggered when first and last cell are 0 ie. NOR
    network.add_gene(10, [{'name': f'CELL_1_Q', 'type': -1, 'Kd': conn_Kds[0], 'n': conn_ns[0]},
                            {'name': f'CELL_{d}_Q','type': -1, 'Kd': conn_Kds[1], 'n': conn_ns[1]}],
                        [{"name": f"INSTRUCTION_1"}])

    # Instructions 2 to D are triggered at a falling edge for respective positions
    # ie. for instruction i the i-th cell is 0, but the one before is 1 
    for i in range(2, d+1):
        network.add_gene(10, [{'name': f'CELL_{i-1}_Q', 'type': 1, 'Kd': conn_Kds[2*i-2], 'n': conn_ns[2*i-2]},
                                {'name': f'CELL_{i}_Q', 'type': -1, 'Kd': conn_Kds[2*i-1], 'n': conn_ns[2*i-1]}],
                            [{"name": f"INSTRUCTION_{i}"}])

    # Instruction D+1 is triggered when first and last cell are 1 ie. AND
    network.add_gene(10, [{'name': f'CELL_1_Q', 'type': 1, 'Kd': conn_Kds[2*d], 'n': conn_ns[2*d]},
                            {'name': f'CELL_{d}_Q','type': 1, 'Kd': conn_Kds[2*d+1], 'n': conn_ns[2*d+1]}],
                        [{"name": f"INSTRUCTION_{d+1}"}])

    # Instructions D+2 to 2*D are triggered at a rising edge for respective positions
    # ie. for instruction D+i the i-th cell is 1, but the one before is 0 
    for i in range(2, d+1):
        network.add_gene(10, [{'name': f'CELL_{i-1}_Q', 'type': -1, 'Kd': conn_Kds[2*d+2*i-2], 'n': conn_ns[2*d+2*i-2]},
                                {'name': f'CELL_{i}_Q',  'type': 1, 'Kd': conn_Kds[2*d+2*i-1], 'n': conn_ns[2*d+2*i-1]}],
                            [{"name": f"INSTRUCTION_{d+i}"}])


def truthgenerator(clks, d):

    # Extract every rising edge
    rising = [i for i in range(len(clks)-1) if clks[i+1]-clks[i]==100]

    # Initialize the array for ground truths
    truth = np.array([])

    for i in range(2*d):

        # Keeping only the rising edges that would activate or deactivate an instruction
        cntr = [rising[pt] for pt in range(len(rising)) if (pt%(2*d)==i or pt%(2*d)==((i+1)%(2*d)))]

        # Removing the first rising edge from last instruction as it's unnecessary
        if i == (2*d)-1:
            cntr = cntr[1:]

        cntr = [0] + cntr + [len(clks)-1]

        inst = np.array([])

        # Construct ground truth based on clock inputs
        for p in range(len(cntr)-1):
            inst = np.concatenate([inst, np.repeat(100-(100 if len(inst)==0 else inst[-1]), cntr[p+1]-cntr[p])])

        inst = np.concatenate([[clks[0]], inst])

        if len(truth) == 0:
            truth = np.array([inst])
        else:
            truth = np.append(truth, [inst], axis=0)
    return truth



if __name__ == "__main__":
    # ------------------------
    # DEBUG / VISUALIZATION
    # ------------------------

    # Initializing the network
    myregister = grn.grn()
    counterregister(myregister, 3)

    # Because of the mpl bug
    faulthandler.enable()

    print(myregister.genes)
    myregister.plot_network()
    plt.show()

    # ------------------------
    # SIMULATION
    # ------------------------

    clks = [0 if i%2==0 else 100 for i in range(20)]
    T, Y = simulator.simulate_sequence(myregister, clks, t_single = 250, plot_on=False)
    spe = np.array(myregister.species_names)
    num = len(spe)
    t = truthgenerator(Y[:, 0], 3)
    t = np.array([t[:, i] for i in range(len(T))])
    plt.plot(T, Y[:, [0]], '--')
    plt.plot(T, Y[:, [7, 15, 23]])
    plt.legend(np.concatenate([spe[[0]], spe[[7, 15, 23]]]))
    plt.show()
    plt.plot(T, Y[:, [0]], '--')
    plt.plot(T, Y[:, num-6:])
    plt.legend(np.concatenate([spe[[0]], spe[num-6:]]))
    plt.show()
    plt.plot(T, t)
    plt.show()