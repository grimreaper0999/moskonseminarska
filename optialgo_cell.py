import pygad
import grn
import simulator
import msdflipflop
import numpy as np
import tqdm

# ------------------------
# HYPERPARAMETERS
# ------------------------

# Basic genetic algo hyperparameters
GENERATIONS = 20
PARENTS_MATING = 4
POPULATION_SIZE = 20
PARENT_SELECTION_TYPE = "rws"
KEEP_PARENTS = 1
CROSSOVER_TYPE = "scattered"
MUTATION_TYPE = "random"
MUTATION_PROBABILITY = 0.3

# Gene value pool initialization
DECAY_VALUES = np.array([0, 0.1, 0.2, 0.5])
KD_VALUES = np.array([0.1, 0.2, 0.5, 1, 2, 5, 10])
N_VALUES = np.array([1, 2, 3, 5, 10])

# Parameters for converting genes to cell parameters
# See msdflipflop.py for valid values
DECAY_GENES = 1
KD_GENES = 2
N_GENES = 2

# Clock inputs, data inputs and ground truth to test the cell behavior
clks = np.concat([np.repeat(0, 5),[100 if i%4==0 else 0 for i in range(10)]])
data = None

# ------------------------
# ALGORITHM PREPARATION
# ------------------------

# Fitness function - determining how good a solution is
def fitness_func(ga_instance, solution, solution_idx):

    sol = solution.astype(np.int32)

    # Preparation of parameters
    decays = [DECAY_VALUES[i] for i in sol[0:DECAY_GENES]]
    Kds = [KD_VALUES[i] for i in sol[DECAY_GENES:DECAY_GENES+KD_GENES]]
    ns = [N_VALUES[i] for i in sol[DECAY_GENES+KD_GENES:]]

    cell = grn.grn()

    # Simulating the cell using given clock and input data, return MSE of Q vs. ground truth
    if not data:
        msdflipflop.registercell("cell", cell, inputname="cell_QBAR", decays=decays, Kds=Kds, ns=ns)
        _, Y = simulator.simulate_sequence(cell, clks, t_single = 250, plot_on=False)
        return -np.mean((Y[:, 7]-msdflipflop.truthgenerator(Y[:, 0]))**2)

    else:
        msdflipflop.registercell("cell", cell, decays=decays, Kds=Kds, ns=ns)
        _, Y = simulator.simulate_sequence(cell, [(data[i], clks[i]) for i in range(len(clks))], t_single = 250, plot_on=False)
        return -np.mean((Y[:, 8]-msdflipflop.truthgenerator(Y[:, 0]))**2)



# ------------------------
# RUNNING THE ALGORITHM
# ------------------------

ga_instance = pygad.GA(num_generations=GENERATIONS,
                    sol_per_pop=POPULATION_SIZE,
                    num_parents_mating=PARENTS_MATING,
                    num_genes=DECAY_GENES + KD_GENES + N_GENES,
                    gene_space=np.repeat([range(len(DECAY_VALUES))], DECAY_GENES, axis=0).tolist()
                                + np.repeat([range(len(KD_VALUES))], KD_GENES, axis=0).tolist()
                                + np.repeat([range(len(N_VALUES))], N_GENES, axis=0).tolist(),
                    fitness_func=fitness_func,
                    parent_selection_type=PARENT_SELECTION_TYPE,
                    keep_parents=KEEP_PARENTS,
                    crossover_type=CROSSOVER_TYPE,
                    mutation_type=MUTATION_TYPE,
                    mutation_probability=MUTATION_PROBABILITY)

ga_instance.run()

temp1 = np.array([DECAY_GENES, KD_GENES, N_GENES])
temp2 = [np.sum(temp1[:i+1]) for i in range(len(temp1))]
solution, solution_fitness, solution_idx = ga_instance.best_solution()
print(("Parameters of the best solution :\n" +
       "Decay values: {decays}\n" +
       "Kd values: {Kds}\n" +
       "n values: {ns}").format(decays=DECAY_VALUES[solution[:temp2[0]].astype(np.int32)],
                                Kds=KD_VALUES[solution[temp2[0]:temp2[1]].astype(np.int32)],
                                ns=N_VALUES[solution[temp2[1]:].astype(np.int32)]))
print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))