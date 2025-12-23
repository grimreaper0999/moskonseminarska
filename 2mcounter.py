import numpy as np
import matplotlib.pyplot as plt
import grn, simulator
import networkx as nx
import faulthandler
from msdflipflop import registercell

# Because of the mpl bug
faulthandler.enable()

# Number of register cells
D = 1

# ------------------------
# BUILDING THE COUNTER
# ------------------------

# Initializing the network
myregister = grn.grn()

# Initializing the clock input (should i build it myself)
myregister.add_input_species("CLK")

# Initialize the first cell with clk the negated output of the last 
registercell("CELL_1", myregister, "CLK", f"CELL_{D}_QBAR")

# Initialize the rest of the cells with clk and output of the previous cell
for c in range(2, D+1):
    registercell(f"CELL_{c}", myregister, "CLK", f"CELL_{c-1}_Q")
    print(f"CELL_{c}")

# ------------------------
# BUILDING THE DECODER
# ------------------------

# Initialize the numbered instruction outputs
for i in range(1, 2*D+1):
    myregister.add_species(f"INSTRUCTION_{i}", 0.1)

# Instruction 1 is triggered when first and last cell are 0 ie. NOR
myregister.add_gene(10, [{'name': f'CELL_1_Q', 'type': -1, 'Kd': 5, 'n': 3},
                         {'name': f'CELL_{D}_Q','type': -1, 'Kd': 5, 'n': 3}],
                    [{"name": f"INSTRUCTION_1"}])

# Instructions 2 to D are triggered at a falling edge for respective positions
# ie. for instruction i the i-th cell is 0, but the one before is 1 
for i in range(2, D+1):
    myregister.add_gene(10, [{'name': f'CELL_{i-1}_Q', 'type': 1, 'Kd': 5, 'n': 3},
                             {'name': f'CELL_{i}_Q',      'type': -1, 'Kd': 5, 'n': 3}],
                        [{"name": f"INSTRUCTION_{i}"}])

# Instruction D+1 is triggered when first and last cell are 1 ie. AND
myregister.add_gene(10, [{'name': f'CELL_1_Q', 'type': 1, 'Kd': 5, 'n': 3},
                         {'name': f'CELL_{D}_Q','type': 1, 'Kd': 5, 'n': 3}],
                    [{"name": f"INSTRUCTION_{D+1}"}])

# Instructions D+2 to 2*D are triggered at a rising edge for respective positions
# ie. for instruction D+i the i-th cell is 1, but the one before is 0 
for i in range(2, D+1):
    myregister.add_gene(10, [{'name': f'CELL_{i-1}_Q', 'type': -1, 'Kd': 5, 'n': 3},
                             {'name': f'CELL_{i}_Q',      'type': 1, 'Kd': 5, 'n': 3}],
                        [{"name": f"INSTRUCTION_{D+i}"}])


# ------------------------
# DEBUG / VISUALIZATION
# ------------------------

print(myregister.genes)
myregister.plot_network()
plt.draw()

# ------------------------
# SIMULATION
# ------------------------

clks = [0 if i%2==0 else 100 for i in range(10)]
T, Y = simulator.simulate_sequence(myregister, clks, t_single = 250)