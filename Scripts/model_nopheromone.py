#########################################################################
## Title: ANT FORAGING BEHAVIOUR - no pheromones                       ##
## Course: Project Computational Biology                               ##
## Authors: Camille Timmermans & Vivian Van Reybrouck Van Gelder       ##
#########################################################################

## Import packages
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt

## Define parameters
arena_size = 20                      # Size of the arena
colony_size = 4                      # Size of the colony
ants = 10                            # Number of foraging ants in the arena
food_sources = 4                     # Number of food sources distributed throughout the arena
food_value = 25                      # Initial value of each food source
steps_per_ant = 1000                 # Maximum number of steps for each ant
runs = 100                           # Number of simulation runs

####### CLASS 1: ARENA #######
class Arena:
    def __init__(self, arena_size=arena_size, colony_size=colony_size, food_sources=food_sources, food_value=food_value):
        # Initialize arena parameters
        self.arena_size = arena_size    # Grid size
        self.colony_size = colony_size  # Colony size
        self.food_sources = {}          # Dictionary to hold food source positions and values
        self.food_value = food_value    # Initial value of food source
        self.init_colony()              # Initialize the colony area
        self.init_food_sources(food_sources)  # Place food sources randomly in the arena

    def init_colony(self):
        # Calculate middle index of the arena
        middle = self.arena_size // 2               # Middle sized index for even sized grid
        col_start = middle - self.colony_size // 2  # Start index (coordinates) of the colony
        col_end = col_start + self.colony_size      # End index (coordinates) of the colony
        self.colony_area = [(x, y) for x in range(col_start, col_end) for y in range(col_start, col_end)]  # Create list of coordinates that represent the colony area

    def init_food_sources(self, food_sources):
        # Randomly place the food sources in the arena
        for _ in range(food_sources):               # For every food source (defined earlier), the while loop will be run until it finds a valid position, in which case it will assign the food value
            while True:                             # Creates an infinite loop that will keep running until it encounters a 'break'
                x = random.randint(0, self.arena_size - 1)  # Generate random x-coordinate
                y = random.randint(0, self.arena_size - 1)  # Generate random y-coordinate
                # Ensure that the food source is not placed in the colony area or on another food source
                if (x, y) not in self.colony_area and (x, y) not in self.food_sources:
                    self.food_sources[(x, y)] = self.food_value
                    break

####### CLASS 2: ANT #######
class Ant:
    def __init__(self, arena):
        self.arena = arena                                  # Reference to the arena instance
        self.position = self.random_start_position()        # Ant's current position
        self.has_food = False                               # Indicates whether the ant is carrying food

    def random_start_position(self):                        # Start from a random point within the colony area
        return random.choice(self.arena.colony_area)

    def move(self):                                         # Decide action based on whether the ant is carrying food
        if not self.has_food:
            self.search_for_food()                          # Ant searches for food
        else:
            self.return_to_colony()                         # Ant returns to the colony with food

    def search_for_food(self):                              # Ant performs a random walk when searching for food
        self.random_walk()
        if self.position in self.arena.food_sources and self.arena.food_sources[self.position] > 0:     # Check if current position has a food source
            self.arena.food_sources[self.position] -= 1     # Decrease food value by 1
            self.has_food = True                            # Ant is now carrying food

    def random_walk(self):
        x, y = self.position
        moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]          # Define possible moves
        while True:
            direction_x, direction_y = random.choice(moves)  # Choose a random move
            new_x = x + direction_x
            new_y = y + direction_y
            # Ensure that the new position is within the arena boundaries
            if 0 <= new_x < self.arena.arena_size and 0 <= new_y < self.arena.arena_size:
                self.position = (new_x, new_y)  # Update position
                break  # Exit loop when a valid movement is found

    def return_to_colony(self):
        x, y = self.position
        colony_positions = self.arena.colony_area  # List of colony positions
        # Find the closest colony position based on Manhattan distance
        closest_colony_position = min(colony_positions, key=lambda c_pos: abs(c_pos[0] - x) + abs(c_pos[1] - y))

        # Determine the direction to move towards the closest colony position
        direction_x = np.sign(closest_colony_position[0] - x)
        direction_y = np.sign(closest_colony_position[1] - y)

        # Adjust direction to prevent diagonal movement
        if direction_x != 0 and direction_y != 0:
            if abs(closest_colony_position[0] - x) > abs(closest_colony_position[1] - y):
                direction_y = 0
            else:
                direction_x = 0

        # Calculate the new position
        new_x = x + direction_x
        new_y = y + direction_y

        # Ensure the new position is within the arena boundaries
        if 0 <= new_x < self.arena.arena_size and 0 <= new_y < self.arena.arena_size:
            self.position = (new_x, new_y)

        # Check if the ant has reached the colony area
        if self.position in self.arena.colony_area:
            self.has_food = False  # Ant has delivered the food
            self.position = self.random_start_position()  # Start a new search from the colony

####### MAIN SIMULATION #######
results = []  # List to store the results of each simulation run

for run in range(runs):
    arena = Arena()  # Initialize the arena for the current run
    ants_list = [Ant(arena) for _ in range(ants)]  # Create a list of ants
    for _ in range(steps_per_ant):
        for ant in ants_list:
            ant.move()  # Move each ant based on its behavior
        # No pheromone update needed since pheromones are not used

    total_remaining_food = sum(arena.food_sources.values())  # Calculate the remaining food
    results.append({'run': run + 1, 'remaining_food': total_remaining_food})  # Store the results

# Print the results of all simulation runs
for result in results:
    print(f"Run {result['run']}: Total remaining food = {result['remaining_food']}")

# Save the results to a CSV file
df_no_pheromones = pd.DataFrame(results)
df_no_pheromones.to_csv('ant_without_pheromones.csv', index=False)

####### GENERATE PLOT #######
# Extract the data for plotting
run_nr = [result['run'] for result in results]  # Run numbers
remaining_food = [result['remaining_food'] for result in results]  # Remaining food values

plt.figure(figsize=(8, 6))  # Set the figure size
plt.plot(run_nr, remaining_food, marker='o', linestyle='-', label='Remaining Food')  # Plot the data
plt.title('Total Remaining Food per Run (No Pheromones)', fontsize=14)
plt.xlabel('Run Number', fontsize=12)
plt.ylabel('Total Remaining Food', fontsize=12)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
