#!/usr/bin/env python
# -*- coding:utf-8 -*-
#@Time  : 2019/11/4 17:41
#@Author: Jianfeng.Sun
#@File  : Genetic Algorithm basic.py

# learning Genetic algorithm

# 求解 0-5 范围内的一段函数的最大数值和最小的数值
import numpy as np
import matplotlib.pyplot as plt
DNA_SIZE = 10;
POP_SIZE = 100;
CROSS_RATE = 0.8;
MUTATTION_RATE = 0.003;
N_GENERATIONS = 200;
X_BOUND = [0,5];
m=X_BOUND [1];   # x upper and lower bounds

# to define a function to get the fitness

def F(x):return np.sin(10*x)*x + np.cos(2*x)*x     # to find the maximum of this function


#m=2 ** np.arange(DNA_SIZE);
# 从后向前取所有的数 arrange 得到的是0 1  2  到9
#n=2 ** np.arange(DNA_SIZE)[::-1];

# find non-zero fitness for selection
#  去让函数都大于0， 这样才比较好计算其概率的数值。。。。。。。
def get_fitness(pred): return pred + 1e-3 - np.min(pred)

# convert binary DNA to decimal and normalize it to a range(0, 5)

def translateDNA(pop): return pop.dot(2 ** np.arange(DNA_SIZE)[::-1]) / float(2**DNA_SIZE-1) * X_BOUND[1]

# 根据 fitness 取选择 pop
# random 去选择，根据累计的概率去选择

def select(pop, fitness):    # nature selection wrt pop's fitness
    idx = np.random.choice(np.arange(POP_SIZE), size=POP_SIZE, replace=True,
                           p=fitness/fitness.sum())
    return pop[idx]
# 基因的交叉
def crossover(parent, pop):     # mating process (genes crossover)
    if np.random.rand() < CROSS_RATE:
        # 表示 从0到 pop_size-1 随机出一个，种群索引去交叉
        i_ = np.random.randint(0, POP_SIZE, size=1)                             # select another individual from pop
        # 随机出DNA 序列对 的交叉点
        cross_points = np.random.randint(0, 2, size=DNA_SIZE).astype(np.bool)   # choose crossover points
        # 选择 第 i 个 进行一个交叉然后 给parent

        parent[cross_points] = pop[i_, cross_points]                            # mating and produce one child
    return parent
# dna 开始变异的环节怎么去写了
def mutate(child):
    for point in range(DNA_SIZE):
        if np.random.rand() < MUTATTION_RATE:
            if child[point] == 0: child[point]=1;
            else: child[point] =0
        return child
# to start the algorithms
# 产生 100行 10 列的随机数 相当于矩阵的形式  array。。。
pop = np.random.randint(2, size=(POP_SIZE, DNA_SIZE))   # initialize the pop DNA

#print(pop[1]) 0-99
# 动态图的使用
plt.ion()       # something about plotting
# to plot this figure
# 取 500 个点取画图
x = np.linspace(*X_BOUND, 200)
plt.plot(x, F(x))

# 使用迭代的方法取计算基因遗传算法
#for _ in range(N_GENERATIONS)
# how to translate a DNA to a bound funtion 因为函数是有具体的自变量的范围的
for _ in range(N_GENERATIONS):
    F_values = F(translateDNA(pop))    # compute function value by extracting DNA

    # something about plotting
    if 'sca' in globals(): sca.remove()
    sca = plt.scatter(translateDNA(pop), F_values, s=200, lw=0, c='red', alpha=0.5); plt.pause(0.05)

    # GA part (evolution)
    fitness = get_fitness(F_values)
    print("Most fitted DNA: ", pop[np.argmax(fitness), :])
    pop = select(pop, fitness)
    pop_copy = pop.copy()
    for parent in pop:

        child = crossover(parent, pop_copy)
        # parent 相当于一个引用，改变了里面的数值
        child = mutate(child)
         # 变异后的直接给第一个父亲了
        parent[:] = child       # parent is replaced by its child


plt.ioff(); plt.show()
