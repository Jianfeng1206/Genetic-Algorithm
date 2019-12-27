# Author Faizan Ali
# Improvements yet needed to be done

import math
from cv2 import cv2
import numpy as np
class Position():
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def calculate_distance(self, pos_B):
        return math.sqrt((self.x - pos_B.x)**2 + (self.y - pos_B.y)**2)


class Node():
    def __init__(self, position, node_type='Obstacle', size=8, color=(0,0,255), thickness=-1):
        self.position = position
        self.type = node_type
        self.size = size
        self.color = color
        self.thickness = thickness
        self.neighbours = [self]
    
    def draw(self, image):
        cv2.circle(image, (self.position.x, self.position.y), self.size, self.color, self.thickness)

    def draw_connections(self, image, line_thickness=1, line_color=(0,0,0)):
        for neighbour in self.neighbours:
            cv2.line(image, (self.position.x,self.position.y), 
                    (neighbour.position.x,neighbour.position.y), line_color, line_thickness)
    
    def add_neighbour(self, node):
        self.neighbours.append(node)


def generate_random_nodes(number, area, padding=20, min_space=10):
    nodes = []
    while len(nodes) < number:
        pos = Position(
                    np.random.randint(padding, area[0]-padding), 
                    np.random.randint(padding, area[1]-padding)
                    )
        
        min_dist = math.inf
        # Finding the minimum distance for this position
        for node in nodes:
            distance = Position.calculate_distance(pos, node.position)
            if distance < min_dist:
                min_dist = distance

        if min_dist >= min_space:
            nodes.append(Node(pos))
    
    return nodes


def connect_nodes(nodes, max_distance=100):
    for node in nodes:
        for n_node in nodes:
            if node is not n_node:
                distance = Position.calculate_distance(node.position, n_node.position)
                if distance < max_distance and n_node not in node.neighbours:
                    node.add_neighbour(n_node)
                    n_node.add_neighbour(node)

            
def select_random_start_goal(nodes, min_dist=50):
    # Randomly pick Start and Goal
    Start = nodes[0]
    Start.type = 'Start'
    Start.color = (255,0,0)
    Start.size = 20
    
    Goal = None
    for node in nodes:
        if Position.calculate_distance(Start.position, node.position) > min_dist:
            Goal = node
            Goal.type = 'Goal'
            Goal.color = (0,255,0)
            Goal.size = 20
            break
    
    return Start, Goal


refPt = []
def click_select_points(image):
    clone = image.copy()

    def mouse_callback(event,x,y,flags,param):
        global refPt
        if event == cv2.EVENT_LBUTTONUP:
            refPt.append(Position(x,y))

    cv2.destroyAllWindows()
    cv2.namedWindow('Select start and goal')
    cv2.setMouseCallback('Select start and goal', mouse_callback)

    prev_len = len(refPt)
    while True:
        cv2.imshow('Select start and goal',clone)
        cv2.waitKey(20)
        
        if len(refPt) > prev_len:
            break
    
    return refPt[-1]


def click_select_start_goal(nodes, image):
    # Selecting start node
    Start_coord = click_select_points(image)

    min_dist = math.inf
    for node in nodes:
        distance = Position.calculate_distance(node.position, Start_coord)
        if distance < min_dist:
            min_dist = distance
            Start = node
    
    Start.type = 'Start'
    Start.color = (255,0,0)
    Start.size = 20

    for node in nodes:
        node.draw(image)
    
    cv2.imshow('Select start and goal', image)
    cv2.waitKey(20)

    # Selecting goal node
    Goal_coord = click_select_points(image)
    
    min_dist = math.inf
    for node in nodes:
        distance = Position.calculate_distance(node.position, Goal_coord)
        if distance < min_dist and node is not Start:
            min_dist = distance
            Goal = node

    Goal.type = 'Goal'
    Goal.color = (0,255,0)
    Goal.size = 20

    for node in nodes:
        node.draw(image)
    
    cv2.imshow('Select start and goal', image)
    cv2.waitKey(20)

    return Start, Goal


def draw_line(pos_A, pos_B, image, thickness=1, color=(0,0,0)):
    cv2.line(image, (pos_A.x, pos_A.y), (pos_B.x, pos_B.y), color, thickness)


def draw_path(path, image, thickness=3, color=(50,100,50)):
    for i in range(len(path)):
        if i > 0:
            draw_line(path[i].position, path[i-1].position, image, thickness, color)


def generate_random_path(Start, Goal, image):
    current_node = Start
    current_node.visited = True
    path_nodes = [current_node]
    while current_node is not Goal:
        neighbours = current_node.neighbours

        # Calculating fitnesses
        fitnesses = []

        legit_node = False
        for i in range(len(neighbours)):
            if neighbours[i] not in path_nodes:
                legit_node = True
                distance = Position.calculate_distance(neighbours[i].position, Goal.position) + 0.00001
                fitnesses.append(1/(distance))
            else:
                fitnesses.append(0.0000000000001)

        if not legit_node:
            current_node = Start
            path_nodes = [current_node]

        else:
            # PDF
            pdf = []
            total_sum = sum(fitnesses)
            for i in range(len(neighbours)):
                pdf.append(fitnesses[i]/total_sum)
            
            # CDF
            cdf = []
            c_sum = 0
            for i in range(len(neighbours)):
                c_sum += pdf[i]
                cdf.append(c_sum)

            # Selecting based on probability
            for i in range(len(neighbours)):
                if np.random.rand() < cdf[i]:
                    current_node = neighbours[i]
                    path_nodes.append(current_node)
                    break

    return path_nodes


def calculate_fitness(path):
    distance = 0.00001
    for i in range(len(path)):
        if i > 0:
            distance += Position.calculate_distance(path[i].position, path[i-1].position)
    
    return 1/(distance**2)


def get_common_nodes(node_A, node_B):
    common_nodes = []
    for node in node_A.neighbours:
        if node in node_B.neighbours:
            common_nodes.append(node)
    
    return common_nodes


if __name__ == "__main__":
    world_size = (640, 480)
    display = np.ones((world_size[1], world_size[0], 3), dtype=np.uint8) * 255 # Display init as white background

    # Creating random nodes list
    nodes = generate_random_nodes(60, world_size, min_space=50)

    # Connecting nodes to neighbours
    connect_nodes(nodes, max_distance=100)

    # Draw nodes
    for node in nodes:
        node.draw(display)

    # Draw connections
    line_thickness = 1
    line_color = (150,0,0)
    for node in nodes:
        node.draw_connections(display, line_thickness, line_color)

    cv2.imshow('display', display)
    cv2.waitKey(15)

    # Start, Goal = select_random_start_goal(nodes, min_dist=240)
    Start, Goal = click_select_start_goal(nodes, display)

    cv2.destroyAllWindows()
    # Draw nodes
    for node in nodes:
        node.draw(display)

    cv2.imshow('display', display)
    cv2.waitKey(15)


    # Generating inital population
    population_size = 100

    gen = []
    max_path_length = 30
    for i in range(population_size):
        debug_image = display.copy()
        path = generate_random_path(Start, Goal, debug_image)
        while len(path) > max_path_length or len(path) == 0:
            path = generate_random_path(Start, Goal, debug_image)

        while len(path) < max_path_length:
            path.append(Goal)

        gen.append(path)

    cv2.destroyAllWindows()
    iterations = 1000
    for it in range(iterations):
        # Check population fitness
        fitnesses = []
        for i in range(population_size):
            fitnesses.append(calculate_fitness(gen[i]))

        # Display chromosomes above fitness threshold
        d = display.copy()
        for i in range(population_size):
            if fitnesses[i] == max(fitnesses):
                color = (np.random.randint(0,256),np.random.randint(0,256),np.random.randint(0,256))
                draw_path(gen[i], d, color=color)

        cv2.imshow('itterations', d)
        cv2.waitKey(20)

        # Generating PDF
        pdf = []
        total_sum = sum(fitnesses)
        for i in range(population_size):
            pdf.append(fitnesses[i]/total_sum)

        # Generating CDF
        cdf = []
        c_sum = 0
        for i in range(population_size):
            c_sum += pdf[i]
            cdf.append(c_sum)

        # Selecting parents based on cdf
        parents = []
        for i in range(population_size):
            for k in range(population_size):
                if np.random.rand() <= cdf[k]:
                    parents.append(gen[k])
                    break
        
        # Generating new population by cross-over and mutation
        fitnesses_copy = fitnesses.copy()
        top_ten = [max(fitnesses_copy)]
        for i in range(10):
            fitnesses_copy.remove(max(fitnesses_copy))
            top_ten.append(max(fitnesses_copy))

        for i in range(population_size):
            if not fitnesses[i] in top_ten:
                parent_A = parents[np.random.randint(0,population_size)]
                parent_B = parents[np.random.randint(0,population_size)]

                gen[i] = parent_A
                
                # Constructive Cross-over
                for c in range(len(gen[i])):
                    if c > 0 and (c+1) < len(gen[i]):
                        if gen[i][c+1] in gen[i][c-1].neighbours:
                            gen[i][c] = gen[i][c+1]

                        if parent_B[c] in gen[i][c+1].neighbours and parent_B[c] in gen[i][c-1].neighbours:
                            distance_A = Position.calculate_distance(gen[i][c-1].position, gen[i][c].position) + Position.calculate_distance(gen[i][c].position, gen[i][c+1].position) 
                            distance_B = Position.calculate_distance(gen[i][c-1].position, parent_B[c].position) + Position.calculate_distance(parent_B[c].position, gen[i][c+1].position) 
                            if distance_B < distance_A:
                                gen[i][c] = parent_B[c]
                
                # Constructive Mutation
                for c in range(len(gen[i])):
                    if c > 0 and (c+1) < len(gen[i]):
                        common_nodes = get_common_nodes(gen[i][c+1], gen[i][c-1])
                        distance_A = Position.calculate_distance(gen[i][c-1].position, gen[i][c].position) + Position.calculate_distance(gen[i][c].position, gen[i][c+1].position) 
                        min_dist = math.inf     
                        if len(common_nodes) > 0:
                            for node in common_nodes:
                                distance_B = Position.calculate_distance(gen[i][c-1].position, node.position) + Position.calculate_distance(node.position, gen[i][c+1].position)
                                if distance_B < min_dist and distance_B < distance_A:
                                    min_dist = distance_B
                                    gen[i][c] = common_nodes[np.random.randint(0,len(common_nodes))]