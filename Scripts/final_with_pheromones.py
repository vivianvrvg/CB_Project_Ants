#######################################################################################
## Title: ANT FORAGING BEHAVIOUR - with pheromones                                   ##
## Course: Project Computational Biology                                             ##
## Authors: Camille Timmermans, Vivian Van Reybrouck Van Gelder & Zaza Terpstra      ##
#######################################################################################

## Import packages
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

## Define parameters
arena_size = 20                     # Size of arena => should be even since middle calculation (for colony) is based on even grid sizes
colony_size = 4                     # Size of grid
ants = 10                           # Number of foraging ants in arena 
food_sources = 4                    # Number of food sources distributed throughout arena
food_value = 25                     # Value of the food source at beginning of run
pheromone_start = 100               # Level of pheromone when just deposited (highest level)
pheromone_min_detectable = 25       # Below this value the ants can no longer detect the pheromone
pheromone_max = 1000                # Maximal deposition of the pheromone
evaporation_rate = 0.10             # Pheromone evaporation rate per step (declines by this proportion) 
diffusion_distance = 2              # Cells until where the pheromone can diffuse
runs = 100                          # Number of replicates


# Build different classes which represent different agents of the model.
# To each agent we attribute characteristis that describe the agent.

####### CLASS 1: ARENA #######
class Arena:
    def __init__(self, arena_size = arena_size, colony_size = colony_size, food_sources = food_sources, food_value = food_value):
        # Initialise the arena pararmeters
        self.arena_size = arena_size    # Grid size
        self.colony_size = colony_size   # Colony size
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
                x = random.randint(0, self.arena_size - 2)  # Generate random x-coordinate ==> since indexing starts at 0, valid indices range from 0 to self.size - 2 (avoid placing on border)
                y = random.randint(0, self.arena_size - 2)  # Generate random y-coordinate
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
        current_value = self.pheromone_grid[x,y]            # Select the current pheromone value
        if current_value < pheromone_max:                   # Allow deposition only when the current value is under the max
            self.pheromone_grid[x, y] = min(current_value + value, pheromone_max)     # adapt the current pheromone value and cap at the maximal (so if new deposit is higher than max, it get's set back to max)             # Here you change the value of the pheromone level at a specific position

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
        moves = [(-1, 0), (1, 0), (0, -1), (0, 1),          # Define possible moves
                 (-1, -1), (-1, 1), (1, -1), (1, 1)]          
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
        valid_positions = []                                 # Initialise next positions
        pheromone_values = []                               # Initialise pheromone values of adjacent cells
        moves = [(-1, 0), (1, 0), (0, -1), (0, 1),          # Define possible moves
                 (-1, -1), (-1, 1), (1, -1), (1, 1)]          

        for direction_x, direction_y in moves:              
            new_x, new_y = x + direction_x, y + direction_y # Movement of ant
            if 0 <= new_x < self.arena.arena_size and 0 <= new_y < self.arena.arena_size:   # Check if the position is within the arena size
                valid_positions.append((new_x, new_y))                                      # update the valid positions
                pheromone_values.append(self.arena.get_pheromone_value((new_x, new_y)))     # update the pheromone values of these positions

        #convert pheromone values to probabilities
        total_pheromone = sum(pheromone_values)
        if total_pheromone > 0:                                  # If there are positions with max pheromone level
            probabilities = [value / total_pheromone for value in pheromone_values]
        else:                                               # If no positions have pheromone, they all get the same probability
            probabilities = [1 / len(valid_positions)] * len(valid_positions)

        self.position = random.choices(valid_positions, weights = probabilities, k = 1)[0]      # Choose a random position based on the probabilities


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
        
        # Calculate the new position by adding direction_x and y to the current position
        new_x = x + direction_x
        new_y = y + direction_y

        # Ensure the new position is within the bounds of the arena grid
        if 0 <= new_x < self.arena.arena_size and 0 <= new_y < self.arena.arena_size:
            self.position = (new_x, new_y)
        
        # Record the current position in the ant's paths for pheromone deposition later
        self.path.append(self.position)

        # deposit pheromone when going back to the colony, only when ant has food
        if self.has_food == True:
            self.leave_pheromone()

        # Check if the ant has reached the colony area
        if self.position in self.arena.colony_area:
            self.has_food = False                           # Ant has delivered the food so has_food is no longer True
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
            self.pheromone_diffusion(pos, pheromone_value)      # Call the recursive diffusion function to deposit and diffuse pheromone

    def pheromone_diffusion(self,position, value, depth = 0):
        if value < pheromone_min_detectable or depth > diffusion_distance:    # Stop the recursion if the value drops below the minimum detectable threshold or when it's too far away from the pheromone trail
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
                        self.pheromone_diffusion((new_x, new_y), half_value, depth + 1)     # Recursively diffuse pheromone to the neighboring cell


####### MAIN SIMULATION #######
results = []                                    # Create an empty list to store the results of every simulation run 
food_left_all_runs = []                         # Create an empty list to store the food left after each 100 steps for every run

for run in range(runs) :
    arena = Arena()                             # Initialise the arena for the current run
    ants_list = [Ant(arena) for _ in range(ants)]    # Create a list of ants, each tied to the arena
    total_steps = 0                             # Initialise the total number of steps taken
    food_left = []                                  # Create an empty list to store the food still in the arena after each 100 steps

    while sum(arena.food_sources.values())> 0:         #run until all food is gone
        for ant in ants_list:                        # Loop for every ant
            ant.move()                          # Move the ant based on its behavior
        arena.update_pheromones()               # Update the pheromone grid to stimulate evaporation
        total_steps += 1                        # Update total number of steps taken

        if total_steps % 100 == 0:              # for each 100 steps in a run
            remaining_food = sum(arena.food_sources.values())           # calculate how much food is left in the arena
            food_left.append({'steps':total_steps, 'remaining_food': remaining_food})  # store the remaining food 

    if sum(arena.food_sources.values()) == 0:  #when all the food is gone, add the last value to the food_left list
        food_left.append({'steps': total_steps, 'remaining_food': 0})

    print('run:', run)

    results.append({'run' : run + 1, 'total_steps': total_steps})  # Store the results
    food_left_all_runs.append({'run': run + 1, 'food_left': food_left}) # Store the remaining food for each 100 steps for every run

# Print the results of all simulation runs
for result in results:
    print(f"Run {result['run']}: Total steps = {result['total_steps']}")

# Store the dataframe as a .csv file in the current working folder
df_pheromones = pd.DataFrame(results)           # dataframe for number of steps until all food is gone
df_pheromones.to_csv('ant_with_pheromones.csv', index = False) # Index=F indicates that we don't want the index column as a column in the final dataset

df_food_tracking = pd.DataFrame(food_left_all_runs)   # dataframe for remaing food after each 100 steps
df_food_tracking.to_csv('food_tracking.csv', index=False)

####### GENERATE PLOT #######
# Extract the data for plotting the amount of steps until all food was cleared
run_nr = [result['run'] for result in results]  # Extract run numbers from the results
total_steps = [result['total_steps'] for result in results]  # Extract remaining food values

plt.figure(figsize=(8, 6))  # Set the figure size
plt.plot(run_nr, total_steps, marker='o', linestyle='-', label='Total steps')  # Plot with markers and lines
plt.title('Steps until all food is consumed, with pheromones', fontsize=14)  # Title for the plot
plt.xlabel('Run Number', fontsize=12)  # Label for the x-axis
plt.ylabel('Total Steps', fontsize=12)  # Label for the y-axis
plt.grid(True)  # Add a grid to the plot
plt.legend()  # Add a legend
plt.tight_layout()  # Adjust layout to prevent clipping

plt.show()

#extract data for plotting the leftover food after each 100 steps
# Plot remaining food after each 100 steps for every run
plt.figure(figsize=(10, 8))

for run_data in food_left_all_runs:  # Iterate over food data for all runs
    run_number = run_data['run']
    food_data = run_data['food_left']

    steps = [entry['steps'] for entry in food_data]
    remaining_food = [entry['remaining_food'] for entry in food_data]

    # Plot each run on the same graph
    plt.plot(steps, remaining_food, marker='o', linestyle='-', label=f'Run {run_number}')

# Configure the graph after all runs are plotted
plt.title('Remaining Food After Each 100 Steps (All Runs), with pheromones', fontsize=16)
plt.xlabel('Steps', fontsize=14)
plt.ylabel('Remaining Food', fontsize=14)
plt.grid(True)
plt.legend(title='Runs', fontsize=12, loc='upper right')
plt.tight_layout()

# Show the combined plot
plt.show()

####### ANIMATION #######        
# Visualization for Arena and Ants 
def plot_arena(ax, arena, ants_list):
    ax.clear()
    ax.set_title("Ant Foraging Simulation, with pheromones")
    ax.set_xlim(0, arena.arena_size)
    ax.set_ylim(0, arena.arena_size)
    
    # Plot colony
    colony_x, colony_y = zip(*arena.colony_area)            # I changed this, since I think the animation was off because our colony size was larger than just one element 
    ax.scatter(colony_x, colony_y, c='blue', s=200, label="Colony", marker='s')

    # Plot food sources
    for (y, x), value in arena.food_sources.items():  # y and x changed to make the animation match up
        if value > 0:
            ax.scatter(x, y, c='green', s=100, label="Food" if 'Food' not in ax.get_legend_handles_labels()[1] else "", marker='o')

    # Plot ants
    ant_positions = [ant.position for ant in ants_list]
    ant_y, ant_x = zip(*ant_positions)    # y and x changed to make the animation match up 
    ax.scatter(ant_x, ant_y, c='red', s=50, label="Ants", marker=(5, 1))

    # Plot pheromone trails
    pheromones = arena.pheromone_grid
    ax.imshow(pheromones, cmap="plasma", origin="lower", alpha=0.5, interpolation="nearest")

    # Add legend
    ax.legend(loc="upper right")

# Animation Function
def animate_simulation(arena, ants_list, steps):  # I changed this to total steps (new code in the main function)
    fig, ax = plt.subplots(figsize=(8, 8))

    def update(frame):
        for ant in ants_list:
            ant.move()
        arena.update_pheromones()
        plot_arena(ax, arena, ants_list)

    anim = FuncAnimation(fig, update, frames = steps, interval=200, repeat=False)
    plt.show()

# Heatmap Function
def plot_pheromone_heatmap(arena):
    plt.figure(figsize=(8, 8))
    plt.title("Pheromone Levels in Arena", fontsize=16)
    plt.imshow(arena.pheromone_grid, cmap="plasma", origin="lower", interpolation="nearest")
    plt.colorbar(label="Pheromone Intensity")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.tight_layout()
    plt.show()

# Main Function
if __name__ == "__main__":            # I changed this since the code wouldn't run and chat GPT said I had to implement 2 __ instead of 1 _
    arena = Arena(colony_size=1)  # Adjusted colony size to 1x1
    ants_list = [Ant(arena) for _ in range(ants)]

    # Run animation
    animate_simulation(arena, ants_list, steps = 1000)

    # Plot final pheromone heatmap
    plot_pheromone_heatmap(arena)



