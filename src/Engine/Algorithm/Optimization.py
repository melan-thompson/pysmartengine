

# 容易陷入局部最优，所以采用多种优化算法
class multiopt:
    def __init__(self, fun, lb, ub):
        self.fun = fun
        self.lb = lb
        self.ub = ub
        self.dim = len(self.lb)
        from math import inf
        self.best_x = inf
        self.best_y = inf

    def DE(self, size_pop=20, max_iter=10):
        from sko.DE import DE
        de = DE(self.fun, self.dim, size_pop=size_pop, max_iter=max_iter, lb=self.lb, ub=self.ub)
        best_x, best_y = de.run()
        if best_y < self.best_y:
            self.best_y = best_y
            self.best_x = best_x

    def GA(self, size_pop=20, max_iter=10):
        from sko.GA import GA
        ga = GA(self.fun, n_dim=self.dim, size_pop=size_pop, max_iter=max_iter, lb=self.lb, ub=self.ub)
        best_x, best_y = ga.run()
        if best_y < self.best_y:
            self.best_y = best_y
            self.best_x = best_x

    def PSO(self, size_pop=20, max_iter=10):
        from sko.PSO import PSO
        pso = PSO(self.fun, n_dim=self.dim, pop=size_pop, max_iter=max_iter, lb=self.lb, ub=self.ub)
        pso.run()
        if pso.best_y < self.best_y:
            self.best_y = pso.best_y
            self.best_x = pso.best_x

    def SA(self):
        from sko.SA import SA, SAFast,SABoltzmann
        print(self.best_x)
        sa = SABoltzmann(self.fun, x0=self.best_x, L=200, max_stay_counter=500, lb=self.lb, up=self.ub)
        best_x, best_y = sa.run()
        if best_y < self.best_y:
            self.best_y = best_y
            self.best_x = best_x