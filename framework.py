# initial commit
########### FRAMEWORK BRAINSTORM ###########

# make matrix 15 x 15
arena_pos = 
arena_food = 
arena_pheromone = 
# nr foragers = 3
# colony = 4 blokjes in matrix
# food source = on 4 locations you want a value of 25, randomly allocated

# ant movement: leaves from random point around colony
start = random waarde rond kolonie

def subtract_food():
    new_food = food_source - 1

# if it is a matrix, check for all positions in the matrix
def random_walk():
    position = start
    for ant in arena :
        while pos != food_source and food_source op pos != 0 :
            new_pos = start + np.random tss -1 en 1
        subtract_food()
        return_walk(new_pos)
# nog in rekening brengen dat mier niet uit de arena kan lopen

def return_walk_pheromones(new_pos):
    """ Route back to colony from food source """
# shortest way back?
# ant needs to ignore other pheromones
# pheromone needs to be dropped and decline as ant moves back to colony
# pheromone needs to diffuse to the surrounding cells
# the cells in which the pheromone diffuse need to decline to according to evaporation rate

def evaporation_rate(x,y):
    """ Rate at which the pheromone declines """

def diffusion_pheromones():
    """ Spreading of the pheromones to surrounding cells """

def found_pheromones(new_pos):
    if arena_pheromones[new_pos] > threshold :
        largest_pheromones = arena_pheromone(pos)
        for i = mogelijk posities rond huidige positie
            if arena_pheromone[i] > current arena_pheromone
            dan verplaatsen
            else random walk terug


#### QUESTIONS ####
## 1. Discrete or continuous time? --> upon reaching the food source and starting return, 
## will the evaporation rate decline continuously or in discrete time steps? Pheromone declines exponentially
## either contunuously or the value declines at every time the ant moves to the next grid. In that case, the time in the evaporation rate
## depends on the time it takes for the ant to move to the next box. Can we choose this time?

## 2. How to define the shortest way back? Compute shortest path algorithm? Or do we calculate the difference between the positions
## and manually make the ant return in terms of -1, 0, 1 vector:
    while x,y != hive.x, hive.y:
        pos -1
   # back in colony

### mwamwamwam