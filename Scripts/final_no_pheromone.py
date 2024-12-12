####################################################################################
## Title: ANT FORAGING BEHAVIOUR - no pheromones                                  ##
## Course: Project Computational Biology                                          ##
## Authors: Camille Timmermans, Vivian Van Reybrouck Van Gelder & Zaza Terpstra   ##
####################################################################################

## Import packages
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

## Define parameters
arena_size = 20                     # Size of arena => should be even since middle calculation (for colony) is based on even grid sizes 
colony_size = 4                     # Size of grid
ants = 10                           # Number of foraging ants in arena 
food_sources = 4                    # Number of food sources distributed throughout arena
food_value = 25                     # Value of the food source at beginning of run
runs = 100                          # Number of replicates 

# Build different classes which represent different agents of the model.
# To each agent we attribute characteristics that describe the agent.

####### CLASS 1: ARENA #######
class Arena:
    def __init__(self, arena_size=arena_size, colony_size=colony_size, food_sources=food_sources, food_value=food_value):
        # Initialise the arena parameters
        self.arena_size = arena_size    # Grid size
        self.colony_size = colony_size   # Colony size
        self.food_sources = {} # Dictionary to hold food source positions and values --> store as follows: {(3, 5): 25} = food source at position (3,5) with value 25
        self.food_value = food_value # Initial value of food source
        self.init_colony() # Initialise the colony area
        self.init_food_sources(food_sources) # Place food sources randomly in the arena
    
    def init_colony(self): 
        # Calculate middle index of arena
        middle = self.arena_size // 2    # Middle-sized index for even-sized grid
        col_start = middle - self.colony_size // 2  # Start index (coordinates) of the colony
        col_end = col_start + self.colony_size      # End index (coordinates) of the colony
        self.colony_area = [(x, y) for x in range(col_start, col_end) for y in range(col_start, col_end)] # Create list of coordinates that represent the colony area
    
    def init_food_sources(self, food_sources):
        # Randomly place the food sources in the arena 
        for i in range(food_sources):                       # For every food source (defined earlier), the while loop will be run until it finds a valid position, in which case it will assign the food value
            while True:                                     # Creates an infinite loop that will keep running until it encounters a 'break'
                x = random.randint(0, self.arena_size - 2)  # Generate random x-coordinate ==> since indexing starts at 0, valid indices range from 0 to self.size - 1
                y = random.randint(0, self.arena_size - 2)  # Generate random y-coordinate
                # Ensure that food source is not randomly placed in the colony area or on another food source
                if (x, y) not in self.colony_area and (x, y) not in self.food_sources:
                    self.food_sources[(x, y)] = self.food_value
                    break

####### CLASS 2: ANT #######
class Ant:
    def __init__(self, arena):
        self.arena = arena                                  # Reference to the arena class
        self.position = self.random_start_position()        # Ant's current position, obtained using function which is defined below
        self.has_food = False                               # Indicator for whether the ant is carrying food, initializes the has_food attribute to False when the ant is created

    def random_start_position(self):                        # Start from a random point within the colony area
        return random.choice(self.arena.colony_area)        # self.arena.colony_area = reference to colony_area attribute of the Arena class

    def move(self):                                         # Decide action based on whether the ant is carrying food
        if not self.has_food:
            self.search_for_food()                          # Ant searches for food
        else:
            self.return_to_colony()                         # Ant returns to the colony with food

    def search_for_food(self):                              # Movement of ant when looking for food
        self.random_walk()                                  # Random walk when looking for food
        # Check if current position has food source
        if self.position in self.arena.food_sources and self.arena.food_sources[self.position] > 0: # If ant is on food source and there is food available then ...
            self.arena.food_sources[self.position] -= 1     # Decrease food value by 1
            self.has_food = True                            # Change has_food attribute to TRUE because ant is now carrying food

    def random_walk(self):                                  # Random walk of the ant
        x, y = self.position
        moves = [(-1, 0), (1, 0), (0, -1), (0, 1),          # Define possible moves
                 (-1, -1), (-1, 1), (1, -1), (1, 1)]          
        while True:                                         # Loop will continue as long as a valid movement has not been found
            direction_x, direction_y = random.choice(moves) # Choose random move
            new_x = x + direction_x
            new_y = y + direction_y
            # Ensure that new position is within the arena boundaries
            if 0 <= new_x < self.arena.arena_size and 0 <= new_y < self.arena.arena_size:
                self.position = (new_x, new_y)              # Update position
                break                                       # Exit loop when a valid movement is found

    def return_to_colony(self):                             # Movement of ant after encountering food source
        x, y = self.position                                # Ant's current position
        # Find closest point in colony area to move towards
        colony_positions = self.arena.colony_area           # Retrieve list of positions that make up colony area
        min_distance = float('inf')                         # Initialise the minimum distance to a large number (infinity)
        closest_colony_position = None                      # Variable to start closest colony position found

        # Loop through each position in colony area to find the closest one based on Manhattan distance 
        for c_pos in colony_positions:
            # Calculate the Manhattan distance
            distance = abs(c_pos[0] - x) + abs(c_pos[1] - y)
            if distance < min_distance:
                min_distance = distance                     # Update minimum distance found
                closest_colony_position = c_pos             # Update the closest colony position

        # Determine the direction to move towards the closest colony position
        direction_x = np.sign(closest_colony_position[0] - x)
        direction_y = np.sign(closest_colony_position[1] - y)

        # Calculate the new position by adding direction_x and y to the current position
        new_x = x + direction_x
        new_y = y + direction_y

        # Ensure the new position is within the bounds of the arena grid
        if 0 <= new_x < self.arena.arena_size and 0 <= new_y < self.arena.arena_size:
            self.position = (new_x, new_y)

        # Check if the ant has reached the colony area
        if self.position in self.arena.colony_area:
            self.has_food = False                           # Ant has delivered the food so has_food is no longer True
            self.position = self.random_start_position()    # Start a new search for food from random position in colony

####### MAIN SIMULATION #######
results = []                                    # Create an empty list to store the results of every simulation run 
food_left_all_runs = []                         # Create an empty list to store the food left after each 100 steps for every run

for run in range(runs):
    arena = Arena()                             # Initialise the arena for the current run
    ants_list = [Ant(arena) for _ in range(ants)]    # Create a list of ants, each tied to the arena
    total_steps = 0                             # Initialise the total number of steps taken
    food_left = []                                  # Create an empty list to store the food still in the arena after each 100 steps

    while sum(arena.food_sources.values()) > 0:         # Run until all food is gone
        for ant in ants_list:                        # Loop for every ant
            ant.move()                          # Move the ant based on its behavior
        total_steps += 1                        # Update total number of steps taken

        if total_steps % 100 == 0:              # For each 100 steps in a run
            remaining_food = sum(arena.food_sources.values())           # Calculate how much food is left in the arena
            food_left.append({'steps': total_steps, 'remaining_food': remaining_food})  # Store the remaining food 

    if sum(arena.food_sources.values()) == 0:  # When all the food is gone, add the last value to the food_left list
        food_left.append({'steps': total_steps, 'remaining_food': 0})

    print('run:', run)

    results.append({'run': run + 1, 'total_steps': total_steps})  # Store the results
    food_left_all_runs.append({'run': run + 1, 'food_left': food_left}) # Store the remaining food for each 100 steps for every run

# Print the results of all simulation runs
for result in results:
    print(f"Run {result['run']}: Total steps = {result['total_steps']}")

# Store the dataframe as a .csv file in the current working folder
df_no_pheromones = pd.DataFrame(results)           # Dataframe for number of steps until all food is gone
df_no_pheromones.to_csv('ant_no_pheromones.csv', index=False) # Index=F indicates that we don't want the index column as a column in the final dataset

# Store the food tracking dataframe (currently stored as a dictionary, not compatible with R)
flattened_data = []

# Loop through each run and expand the nested structure
for run_data in food_left_all_runs:
    run_number = run_data['run']
    for record in run_data['food_left']:
        flattened_data.append({
            'run': run_number,
            'after_nr_steps': record['steps'],
            'remaining_food': record['remaining_food']
        })

# Convert the flattened data into a DataFrame
df_flattened = pd.DataFrame(flattened_data)
# Save the resulting DataFrame to a CSV file
df_flattened.to_csv('food_tracking_no_pheromones.csv', index=False)
# Print the first few rows of the DataFrame for verification
print(df_flattened.head())

####### GENERATE PLOT #######
# Extract the data for plotting the amount of steps until all food was cleared
run_nr = [result['run'] for result in results]  # Extract run numbers from the results
total_steps = [result['total_steps'] for result in results]  # Extract total steps values

# Plot for total steps
plt.figure(figsize=(8, 6))  # Set the figure size
plt.plot(run_nr, total_steps, marker='o', linestyle='-', label='Total steps')  # Plot with markers and lines
plt.title('Steps until all food is consumed, no pheromones', fontsize=14)  # Title for the plot
plt.xlabel('Run Number', fontsize=12)  # Label for the x-axis
plt.ylabel('Total Steps', fontsize=12)  # Label for the y-axis
plt.grid(True)  # Add a grid to the plot
plt.legend()  # Add a legend
plt.tight_layout()  # Adjust layout to prevent clipping
plt.show()

# Extract data for plotting the leftover food after each 100 steps
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
plt.title('Remaining Food After Each 100 Steps (All Runs), no pheromones', fontsize=16)
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
    ax.set_title("Ant Foraging Simulation, no pheromones")
    ax.set_xlim(0, arena.arena_size)
    ax.set_ylim(0, arena.arena_size)
    
    # Plot colony
    colony_x, colony_y = zip(*arena.colony_area)  # Adjust for colony size
    ax.scatter(colony_x, colony_y, c='blue', s=200, label="Colony", marker='s')

    # Plot food sources
    for (y, x), value in arena.food_sources.items():  # Adjusted for correct coordinates
        if value > 0:
            ax.scatter(x, y, c='green', s=100, label="Food" if 'Food' not in ax.get_legend_handles_labels()[1] else "", marker='o')

    # Plot ants
    ant_positions = [ant.position for ant in ants_list]
    ant_y, ant_x = zip(*ant_positions)  # Adjusted for correct coordinates
    ax.scatter(ant_x, ant_y, c='red', s=50, label="Ants", marker=(5, 1))

    # Add legend
    ax.legend(loc="upper right")

# Animation Function
def animate_simulation(arena, ants_list, steps):  # Adjusted for the model without pheromones
    fig, ax = plt.subplots(figsize=(8, 8))

    def update(frame):
        for ant in ants_list:
            ant.move()  # Move ants based on current model behavior
        plot_arena(ax, arena, ants_list)  # Plot arena without pheromone trails

    anim = FuncAnimation(fig, update, frames=steps, interval=200, repeat=False)
    plt.show()

# Main Function
if __name__ == "__main__":  # Adjusted to ensure the code runs correctly
    arena = Arena(colony_size=1)  # Adjusted colony size to 1x1
    ants_list = [Ant(arena) for _ in range(ants)]

    # Run animation
    animate_simulation(arena, ants_list, steps=1000)  # Animate for 1000 steps
