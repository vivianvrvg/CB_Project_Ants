#########################################################################
## Title: ANT FORAGING BEHAVIOUR - with pheromones                     ##
## Course: Project Computational Biology                               ##
## Authors: Camille Timmermans & Vivian Van Reybrouck Van Gelder       ##
#########################################################################

## Adaptations compared to previous model:
##  - Added code to store the obtained dataframe in the folder of your computer
##  - Added code to ensure that an ant cannot get stuck against the arena borders (break after if statement in random_walk)
##  - Added code to let pheromone diffuse further at the moment of deposition
##  - I haven't included diagonal movement (yet)
##  - Discuss diffusion issue

## Import packages
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt

## Define parameters
arena_size = 20                      # Size 16x16 => should be even since middle calculation (for colony) is based on even grid sizes 16
colony_size = 4                      # Size 2x2
ants = 10                            # Number of foraging ants in arena 3
food_sources = 4                     # Number of food sources distributed throughout arena
food_value = 25                      # Value of the food source at beginning of run
steps_per_ant = 1000                  # Maximum number of steps for each ant (every ant has done 100 steps = end of run) 100
pheromone_start = 100                # Level of pheromone when just deposited (highest level)
pheromone_min_detectable = 25        # Below this value the ants can no longer detect the pheromone
evaporation_rate = 0.05              # Pheromone evaporation rate per step (declines by 5%) ===> CHANGE WITH EQUATION
runs = 100                           # Number of replicates 50

# Build different classes which represent different agents of the model.
# To each agent we attribute characteristis that describe the agent.

####### CLASS 1: ARENA #######
class Arena:
    def __init__(self, arena_size = arena_size, colony_size = colony_size, food_sources = food_sources, food_value = food_value):
        # Initialise the arena pararmeters
        self.arena_size = arena_size    # Grid size
        self.colony_size =colony_size   # Colony size
        self.pheromone_grid = np.zeros((arena_size, arena_size)) # Stores the pheromone levels for each position in the arena
        self.food_sources = {} # Dictionary to hold food source positions and values --> store as follows: {(3, 5): 25} = food source at position (3,5) with value 25
        self.food_value = food_value # Initial value of food source
        self.init_colony() # Initialise the colony area
        self.init_food_sources(food_sources) # Place food sources randomly in the arena
    
    def init_colony(self): 
        # Calculate middle index of arena
        middle = self.arena_size // 2    # Middle sized index for even sized grid
        col_start = middle - self.colony_size // 2  # Start index (coordinates) of the colony
        col_end = col_start + self.colony_size      # End index (coordinates) of the colony
        self.colony_area = [(x,y) for x in range(col_start, col_end) for y in range(col_start, col_end)] # Create list of coordinates that represent the colony area
    
    def init_food_sources(self, food_sources):
        # Randomly place the food sources in the arena 
        for i in range(food_sources):                       # For every food source (defined earlier), the while loop will be run until it finds a valid position, in which case it will assign the food value
            while True:                                     # Creates an infinite loop that will keep running until it encounters a 'break'
                x = random.randint(0, self.arena_size - 1)  # Generate random x-coordinate ==> since indexing starts at 0, valid indices range from 0 to self.size - 1
                y = random.randint(0, self.arena_size - 1)  # Generate random y-coordinate
                # Ensure that food source is not randomly placed in the colony area or on another food source
                if (x,y) not in self.colony_area and (x,y) not in self.food_sources:
                    self.food_sources[(x,y)] = self.food_value
                    break

    def update_pheromones(self):
        # Evaporate pheromones by reducing the values based on the previously defined evaporation rate
        self.pheromone_grid *= (1 - evaporation_rate)       # Multiply each element in the pheromone grid by the value on the right-hand side and updates the array (1 - evaporation_rate = fraction that remains)
        # Remove pheromones below minimum detectable threshold
        self.pheromone_grid[self.pheromone_grid < pheromone_min_detectable] = 0
    
    def get_pheromone_value(self, position):                # Retrieve the current pheromone level at a specific position in the arena grid (used later on)
        x, y = position                                     # The position refers to the value that will be given when the funsction is called
        return self.pheromone_grid[x,y]
    
    def set_pheromone_value(self, position, value):         # Set or update the pheromone concentration at a specific position in the arena grid (used later on)
        x, y = position                                     # Same as above, when the function is called, the values for 'position' and 'value' will be given
        self.pheromone_grid[x, y] = value                   # Here you change the value of the pheromone level at a specific position

####### CLASS 2: ANT #######
class Ant:
    def __init__(self, arena):
        self.arena = arena                                  # Reference to the arena class
        self.position = self.random_start_position()        # Ant's current position, obtained using function which is defined below
        self.has_food = False                               # Indicator for whether the ant is carrying food, initializes the has_food attribute to False when the ant is created
        self.path = []                                      # Path taken by the ant which is used for leaving the pheromones

    def random_start_position(self):                        # Start from a random point within the colony area
        return random.choice(self.arena.colony_area)        # self.arena.colony_area = reference to colony_area attribute of the Arena class

    def move(self):                                         # Decide action based on whether the ant is carrying food
        if not self.has_food:
            self.search_for_food()                          # Ant searches for food
        else:
            self.return_to_colony()                         # Ant returns to the colony with food

    def search_for_food(self):                              # Movement of ant when looking for food
        # Check presence of pheromones
        pheromone_value = self.arena.get_pheromone_value(self.position)
        if pheromone_value >= pheromone_min_detectable:
            self.follow_pheromones()                        # Follow pheromone trail
        else:
            self.random_walk()                              # If no pheromone (or not detectable), initiate random walk
        # Check if current position has food source
        if self.position in self.arena.food_sources and self.arena.food_sources[self.position] > 0: # If ant is on food source and there is food available then ...
            self.arena.food_sources[self.position] -= 1     # Decrease food value by 1
            self.has_food = True                            # Change has_food attribute to TRUE because ant is now carrying food
            self.path = [self.position]                     # Start recording path to leave the pheromones
    
    def random_walk(self):                                  # Random walk of the ant when no pheromone is detected
        x, y = self.position
        moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]          # Define possible moves
        while True:                                         # Loop will continue as long as a valid movement has not been found
            direction_x, direction_y = random.choice(moves) # Choose random move
            new_x = x +  direction_x
            new_y = y + direction_y
        # Ensure that new position is within the arena boundaries, since ants cannot move outside of the arena
            if 0 <= new_x < self.arena.arena_size and 0 <= new_y < self.arena.arena_size:
                self.position = (new_x, new_y)              # Update position
                break                                       # Exit loop when a valid movement is found (so ant can no longer remain at same position when it bumps into the wall)

    def follow_pheromones(self):                            # Movement of ants based on pheromone level
        # Ant moves to adjacent cell with highest pheromone value
        x, y = self.position                                # Ant's current position
        max_pheromone = 0                                   # Variable to keep track of the highest pheromone value found among adjacent cells
        next_positions = []                                 # Initialise next positions
        moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]          # Define cardinal directions
        for direction_x, direction_y in moves:              
            new_x, new_y = x + direction_x, y + direction_y # Movement of ant
            if 0 <= new_x < self.arena.arena_size and 0 <= new_y < self.arena.arena_size:   # Check if the position is within the arena size
                pheromone_value = self.arena.get_pheromone_value((new_x, new_y))
                if pheromone_value > max_pheromone:         # Update if higher pheromone value is found: moves to the adjacent cell with the highest pheromone value, even if it's lower than the pheromone level at its current position
                    max_pheromone = pheromone_value         # Update max pheromone level
                    next_positions = [(new_x, new_y)]       # Reset list with this position
                elif pheromone_value == max_pheromone:      # If pheromone level equals max
                    next_positions.append((new_x, new_y))   # Add this position to the list of options
        if next_positions:                                  # If there are positions with max pheromone level
            self.position = random.choice(next_positions)   # Move to one randomly
        else:                                               # If no positions have pheromone > current position
            self.random_walk()                              # Perform a random walk

    def return_to_colony(self):                             # Movement of ant after encountering food source
        x, y = self.position                                # Ant's current position
        # Find closes point in colony area to move towards
        colony_positions = self.arena.colony_area           # Retrieve list of positions that make up colony area
        min_distance = float('inf')                         # Initialise the minimum distance to a large number (infinity)
        closest_colony_position = None                      # Variable to start closest colony position found

        # Loop through each position in colony area to find the closest one based on Manhattan distance 
        # = gives the exact number of steps required to reach a destination (when only cardinal movement are allowed)
        for c_pos in colony_positions:
            # Calculate the Manhattan distance (sum of absolute differences in x and y) to the colony position
            distance = abs(c_pos[0] - x) + abs(c_pos[1] - y)    # Take absolute values of x an y positions of the colony
            if distance < min_distance:
                min_distance = distance                     # Update minimum distance found
                closest_colony_position = c_pos             # Update the closest colony position, no else statement because if untrue, it doesn't need to be updated

        # Determine the direction to move towards the closest colony position
        direction_x = np.sign(closest_colony_position[0] - x)
        direction_y = np.sign(closest_colony_position[1] - y)

        # Adjust direction x and y to prevent diagonal movement by choosing to move along one axis only
        if direction_x != 0 and direction_y != 0:            # If both x and y are > 0, the ant would need to move diagonally. In that case pick one axis based on smallest distance to colony.
            # Decide which action to move along based on greater remaining distance
            if abs(closest_colony_position[0] - x) > abs(closest_colony_position[1] - y):
                direction_y = 0
            else:
                direction_x = 0
        
        # Calculate the new position by adding direction_x and y to the current position
        new_x = x + direction_x
        new_y = y + direction_y

        # Ensure the new position is within the bounds of the arena grid
        if 0 <= new_x < self.arena.arena_size and 0 <= new_y < self.arena.arena_size:
            self.position = (new_x, new_y)
        
        # Record the current position in the ant's paths for pheromone deposition later
        self.path.append(self.position)

        # Check if the ant has reached the colony area
        if self.position in self.arena.colony_area:
            self.has_food = False                           # Ant has delivered the food so has_food is no longer True
            self.leave_pheromone()                          # Ant depositis pheromones when going back to colony
            self.path = []                                  # Reset the path since ant will start new foraging trip
            self.position = self.random_start_position()    # Start a new search for food from random position is colony
    
    def leave_pheromone(self):                              
        trail_length = len(self.path)                           # Get the total number of positions in the recorded path.
        if trail_length == 0:                                   # If the path is empty (ant didn't move), no pheromone is deposited.
            return  # Exit the function.
        for i, pos in enumerate(self.path):                     # Loop through each position in the ant's path
            pheromone_value = pheromone_start * np.exp(-i / trail_length)     # Calculate the pheromone value for the current position, decreases exponentially as it moves along the path
            if pheromone_value < pheromone_min_detectable:      # Skip depositing pheromone if the calculated value is below the detectable threshold
                continue
            self.pheromone_diffusion(pos, pheromone_value)       # Call the recursive diffusion function to deposit and diffuse pheromone

    def pheromone_diffusion(self,position, value):
        if value < pheromone_min_detectable:                    # Stop the recursion if the value drops below the minimum detectable threshold
            return
        x, y = position                                         # Update the pheromone value at the current position
        existing_value = self.arena.get_pheromone_value(position)  # Get existing pheromone value
        self.arena.set_pheromone_value(position, max(existing_value, value))  # Set the max of existing and new value
        half_value = value / 2                                  # Calculate the diffused value for neighboring cells --> each neighbor gets 50% of the current value
        if half_value >= pheromone_min_detectable:              # If the diffused value is still detectable, continue diffusion.
            for direction_x in [-1, 0, 1]:                      # Loop through all directions in x-axis
                for direction_y in [-1, 0, 1]:                  # Loop through all directions in y-axis
                    if direction_x == 0 and direction_y == 0:   # Skip the current cell
                        continue
                    new_x, new_y = x + direction_x, y + direction_y     # Calculate the new position
                    if 0 <= new_x < self.arena.arena_size and 0 <= new_y < self.arena.arena_size:       # Ensure the new position is within the arena boundaries
                        self.pheromone_diffusion((new_x, new_y), half_value)     # Recursively diffuse pheromone to the neighboring cell


####### MAIN SIMULATION #######
results = []                                    # Create an empty list to store the results of every simulation run

for run in range(runs) :
    arena = Arena()                             # Initialise the arena for the current run
    ants_list = [Ant(arena) for _ in range(ants)]    # Create a list of ants, each tied to the arena
    for _ in range(steps_per_ant) :             # Perform the specified nr of steps/ant
        for ant in ants_list:                        # Loop for every ant
            ant.move()                          # Move the ant based on its behavior
        arena.update_pheromones()               # Update the pheromone grid to stimulate evaporation

    total_remaining_food = sum(arena.food_sources.values()) # Calculate the food left in the arena 
    results.append({'run' : run + 1, 'remaining_food': total_remaining_food}) # Store the results

# Print the results of all simulation runs
for result in results:
    print(f"Run {result['run']}: Total remaining food = {result['remaining_food']}")

# Store the dataframe as a .csv file in the current working folder
df_pheromones = pd.DataFrame(results)
df_pheromones.to_csv('ant_with_pheromones.csv', index = False) # Index=F indicates that we don't want the index column as a column in the final dataset

####### GENERATE PLOT #######
# Extract the data for plotting
run_nr = [result['run'] for result in results]  # Extract run numbers from the results
remaining_food = [result['remaining_food'] for result in results]  # Extract remaining food values

plt.figure(figsize=(8, 6))  # Set the figure size
plt.plot(run_nr, remaining_food, marker='o', linestyle='-', label='Remaining Food')  # Plot with markers and lines
plt.title('Total Remaining Food per Run', fontsize=14)  # Title for the plot
plt.xlabel('Run Number', fontsize=12)  # Label for the x-axis
plt.ylabel('Total Remaining Food', fontsize=12)  # Label for the y-axis
plt.grid(True)  # Add a grid to the plot
plt.legend()  # Add a legend
plt.tight_layout()  # Adjust layout to prevent clipping

plt.show()

        

