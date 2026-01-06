import pygad
import grn
import simulator
import counter2m
import numpy as np
import tqdm

# ------------------------
# HYPERPARAMETERS
# ------------------------

# Counter size
D = 3

# Basic genetic algo hyperparameters
GENERATIONS = 10
PARENTS_MATING = 4
POPULATION_SIZE = 20
PARENT_SELECTION_TYPE = "sss"
KEEP_PARENTS = 2
CROSSOVER_TYPE = "scattered"
MUTATION_TYPE = "random"
MUTATION_PROBABILITY = 0.2

# Gene value pool initialization
CELL_DECAY_VALUES = np.array([0, 0.1, 0.2, 0.5])
CELL_KD_VALUES = np.array([0.1, 0.2, 0.5, 1, 2, 5, 10])
CELL_N_VALUES = np.array([1, 2, 3, 5, 10])
INSTR_DECAY_VALUES = np.array([0, 0.1, 0.2, 0.5])
CONN_KD_VALUES = np.array([0.1, 0.2, 0.5, 1, 2, 5, 10])
CONN_N_VALUES = np.array([1, 2, 3, 5, 10])

# Parameters for converting genes to cell parameters
# See msdflipflop.py and 2mcounter.py for valid values
CELL_DECAY_GENES = 1
CELL_KD_GENES = 2
CELL_N_GENES = 2
INSTR_DECAY_GENES = 1
CONN_KD_GENES = 2
CONN_N_GENES = 2

# Clock inputs and ground truth to test the cell behavior
clks = np.concat([np.repeat(0, 5),[100 if i%4==0 else 0 for i in range(13)]])

# ------------------------
# ALGORITHM PREPARATION
# ------------------------

# Fitness function - determining how good a solution is
def fitness_func(ga_instance, solution, solution_idx):

    sol = solution.astype(np.int32)

    # Preparation of parameters
    temp1 = np.array([CELL_DECAY_GENES, CELL_KD_GENES, CELL_N_GENES, INSTR_DECAY_GENES, CONN_KD_GENES, CONN_N_GENES])
    temp2 = [np.sum(temp1[:i+1]) for i in range(len(temp1))]
    cell_decays = [CELL_DECAY_VALUES[i] for i in sol[0:temp2[0]]]
    cell_Kds = [CELL_KD_VALUES[i] for i in sol[temp2[0]:temp2[1]]]
    cell_ns = [CELL_N_VALUES[i] for i in sol[temp2[1]:temp2[2]]]
    instr_decays=[INSTR_DECAY_VALUES[i] for i in sol[temp2[2]:temp2[3]]]
    conn_Kds=[CONN_KD_VALUES[i] for i in sol[temp2[3]:temp2[4]]]
    conn_ns=[CONN_N_VALUES[i] for i in sol[temp2[4]:]]

    register = grn.grn()

    counter2m.counterregister(register, D, None, cell_decays, cell_Kds, cell_ns, instr_decays, conn_Kds, conn_ns)
    _, Y = simulator.simulate_sequence(register, clks, t_single = 250, plot_on=False)
    ground_truth = counter2m.truthgenerator(Y[:, 0], D)
    
    err = 0
    for i in range(-2*D, 0):
        err += -np.mean((Y[:, i]-ground_truth)**2)
    
    return err 



# ------------------------
# RUNNING THE ALGORITHM
# ------------------------

ga_instance = pygad.GA(num_generations=GENERATIONS,
                    sol_per_pop=POPULATION_SIZE,
                    num_parents_mating=PARENTS_MATING,
                    num_genes=CELL_DECAY_GENES + CELL_KD_GENES + CELL_N_GENES + INSTR_DECAY_GENES + CONN_KD_GENES + CONN_N_GENES,
                    gene_space=np.repeat([range(len(CELL_DECAY_VALUES))], CELL_DECAY_GENES, axis=0).tolist()
                                + np.repeat([range(len(CELL_KD_VALUES))], CELL_KD_GENES, axis=0).tolist()
                                + np.repeat([range(len(CELL_N_VALUES))], CELL_N_GENES, axis=0).tolist()
                                + np.repeat([range(len(INSTR_DECAY_VALUES))], INSTR_DECAY_GENES, axis=0).tolist()
                                + np.repeat([range(len(CONN_KD_VALUES))], CONN_KD_GENES, axis=0).tolist()
                                + np.repeat([range(len(CONN_N_VALUES))], CONN_N_GENES, axis=0).tolist(),
                    fitness_func=fitness_func,
                    parent_selection_type=PARENT_SELECTION_TYPE,
                    keep_parents=KEEP_PARENTS,
                    crossover_type=CROSSOVER_TYPE,
                    mutation_type=MUTATION_TYPE,
                    mutation_probability=MUTATION_PROBABILITY)

ga_instance.run()

temp1 = np.array([CELL_DECAY_GENES, CELL_KD_GENES, CELL_N_GENES, INSTR_DECAY_GENES, CONN_KD_GENES, CONN_N_GENES])
temp2 = [np.sum(temp1[:i+1]) for i in range(len(temp1))]
solution, solution_fitness, solution_idx = ga_instance.best_solution()
print(("Parameters of the best solution :\n" +
       "Cell decay values: {cedecay}\n" +
       "Cell Kd values: {ceKds}\n" +
       "Cell n values: {cens}\n" +
       "Instruction decay values: {idecay}\n" +
       "Connection Kd values: {coKds}\n" +
       "Connection n values: {cons}").format(cedecay=CELL_DECAY_VALUES[(solution[:temp2[0]].astype(np.int32))],
                                             ceKds=CELL_KD_VALUES[solution[temp2[0]:temp2[1]].astype(np.int32)],
                                             cens=CELL_N_VALUES[solution[temp2[1]:temp2[2]].astype(np.int32)],
                                             idecay=INSTR_DECAY_VALUES[solution[temp2[2]:temp2[3]].astype(np.int32)],
                                             coKds=CONN_KD_VALUES[solution[temp2[3]:temp2[4]].astype(np.int32)],
                                             cons=CONN_N_VALUES[solution[temp2[4]:].astype(np.int32)]))
print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))