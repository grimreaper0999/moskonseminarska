import pygad
import grn
import simulator
import msdflipflop
import numpy as np

# ------------------------
# HYPERPARAMETERS
# ------------------------

# Basic genetic algo hyperparameters
GENERATIONS = 1000
PARENTS_MATING = 500
POPULATION_SIZE = 2000
PARENT_SELECTION_TYPE = "sss"
KEEP_PARENTS = 200
CROSSOVER_TYPE = "scattered"
MUTATION_TYPE = "random"
MUTATION_PROBABILITY = 0.2

# Gene value pool initialization
DECAY_VALUES = [0, 0.1, 0.2, 0.5]
KD_VALUES = [0.1, 0.2, 0.5, 1, 2, 5, 10]
N_VALUES = [1, 2, 3, 5, 10]

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

    print(solution)

    # Preparation of parameters
    decays = [DECAY_VALUES[i] for i in solution[0:DECAY_GENES]]
    Kds = [KD_VALUES[i] for i in solution[DECAY_GENES:DECAY_GENES+KD_GENES]]
    ns = [N_VALUES[i] for i in solution[DECAY_GENES+KD_GENES:]]

    cell = grn.grn()

    # Simulating the cell using given clock and input data, return MSE of Q vs. ground truth
    if not data:
        msdflipflop.registercell("cell", cell, inputname="cell_QBAR", decays=decays, Kds=Kds, ns=ns)
        _, Y = simulator.simulate_sequence(cell, clks, t_single = 250, plot_on=False)
        return -np.mean((Y[:, 7]-msdflipflop.truthgenerator(Y[0]))**2)

    else:
        msdflipflop.registercell("cell", cell, decays=decays, Kds=Kds, ns=ns)
        _, Y = simulator.simulate_sequence(cell, [(data[i], clks[i]) for i in range(len(clks))], t_single = 250, plot_on=False)
        return -np.mean((Y[:, 8]-msdflipflop.truthgenerator(Y[0]))**2)



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

solution, solution_fitness, solution_idx = ga_instance.best_solution()
print("Parameters of the best solution : {solution}".format(solution=solution))
print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))