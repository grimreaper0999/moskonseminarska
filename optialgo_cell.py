import pygad
import grn

# ------------------------
# INITIALIZATION
# ------------------------

#1. možnost - optimizacija parametrov za obstoječo topologijo npr. degradation rate, Kd, n
#2. možnost, za topologijo:
#mutacija - odvzem ali generiranje nove povezave
#križanje - pol povezav od enega, pol od drugega, nov osebek


# Fitness function - determining how good a solution is
def fitness_func(ga_instance, solution, solution_idx):
    return 0

# ------------------------
# RUNNING THE ALGORITHM
# ------------------------

ga_instance = pygad.GA(num_generations=num_generations,
                       num_parents_mating=num_parents_mating,
                       fitness_func=fitness_func,
                       sol_per_pop=sol_per_pop,
                       num_genes=num_genes,
                       init_range_low=init_range_low,
                       init_range_high=init_range_high,
                       parent_selection_type=parent_selection_type,
                       keep_parents=keep_parents,
                       crossover_type=crossover_type,
                       mutation_type=mutation_type,
                       mutation_percent_genes=mutation_percent_genes)

ga_instance.run()

solution, solution_fitness, solution_idx = ga_instance.best_solution()
print("Parameters of the best solution : {solution}".format(solution=solution))
print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))