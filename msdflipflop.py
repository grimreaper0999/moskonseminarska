import numpy as np
import grn, simulator

cell = grn.grn()

cell.add_input_species("X")
cell.add_input_species("CLK")

#cell.add_species("Y", 0.1)
cell.add_species("MNAND1", 0.1)
cell.add_species("MNAND2", 0.1)
cell.add_species("MNAND3", 0.1)
cell.add_species("MNAND4", 0.1)

# prvi master nand gate oz. not x or not clk
mnand1_reg = [{'name': 'X', 'type': -1, 'Kd': 5, 'n': 3},
             {'name': 'CLK', 'type': -1, 'Kd': 5, 'n': 3}]
cell.add_gene(10, nand1_reg, [{"name":"MNAND1"}], logic_type="or")

# drugi master nand gate z neg. X oz. x or not clk
mnand2_reg = [{'name': 'X', 'type': 1, 'Kd': 5, 'n': 2},
             {'name': 'CLK', 'type': -1, 'Kd': 5, 'n': 3}]
cell.add_gene(10, nand2_reg, [{"name":"MNAND2"}], logic_type="or")

# tretji master nand gate
mnand2_reg = [{'name': 'MNAND1', 'type': -1, 'Kd': 5, 'n': 3},
             {'name': 'MNAND4', 'type': -1, 'Kd': 5, 'n': 3}]
cell.add_gene(10, nand2_reg, [{"name":"MNAND3"}], logic_type="or")

# četrti master nand gate
mnand2_reg = [{'name': 'MNAND2', 'type': -1, 'Kd': 5, 'n': 3},
             {'name': 'MNAND3', 'type': -1, 'Kd': 5, 'n': 3}]
cell.add_gene(10, nand2_reg, [{"name":"MNAND4"}], logic_type="or")

print(cell.genes)

cell.plot_network()