from math import sqrt
from os import makedirs, path
import numpy as np
import matplotlib.pyplot as plt


class Sceme:
    def __init__(self, name, alpha_1_m1, alpha_0_m1, alpha_0_0, alpha_m1_0):
        self.name = name
        self.alpha_1_m1 = alpha_1_m1
        self.alpha_0_m1 = alpha_0_m1
        self.alpha_0_0 = alpha_0_0
        self.alpha_m1_0 = alpha_m1_0

    def return_weights(self):
        return (self.alpha_1_m1, self.alpha_0_m1, self.alpha_0_0, self.alpha_m1_0)

    def get_name(self):
        return self.name


# Начальные условия
def phi(x):
    if x < 0.4 or x > 0.6:
        return 0
    return sqrt(1-100*(x-0.5)*(x-0.5))


def plot_graph(data_to_graph, graph_name):
    N = list()
    X = list()
    Y = list()

    for i in range(data_to_graph.shape[1]):
        X.append(data_to_graph[0][i])
        Y.append(data_to_graph[1][i])
        N.append("begin")

    for i in range(data_to_graph.shape[1]):
        X.append(data_to_graph[0][i])
        Y.append(data_to_graph[2][i])
        N.append("end")

    fig, ax = plt.subplots()

    for n in set(N):
        x1 = [i for i, j in zip(X, N) if j == n]
        y1 = [i for i, j in zip(Y, N) if j == n]
        ax.plot(x1, y1, label=n)

    plt.xlabel('x')
    plt.ylabel('f(x)')

    plt.grid(visible=True)
    plt.legend()
    plt.savefig("graphs/"+graph_name+".png", dpi=300)
    # plt.show()


def calculate_scheme(scheme, N, h, tau, time_steps):

    u = np.zeros((3, N), dtype=np.double)
    for i in range(N):
        # задаём начальные условия, предполагая, что -1 и 0 слой одинаковые
        u[0][i] = phi(h*i)
        u[1][i] = u[0][i]
    weight = scheme.return_weights()
    for n in range(1, time_steps):
        for m in range(1, N):
            u[2][m] = weight[0]*u[2][m-1] + weight[1] * \
                u[1][m-1] + weight[2]*u[1][m] + weight[3]*u[0][m]

        u[[0, 1]] = u[[1, 0]]
        u[[1, 2]] = u[[2, 1]]
        u[2] = np.zeros(N)

    # np.savetxt(scheme.name()+".csv", u, delimiter=",")

    data_to_graph = np.zeros((3, N), dtype=np.double)
    print(scheme.get_name(), weight)
    for i in range(N):
        data_to_graph[0][i] = i*h
        data_to_graph[1][i] = phi(h*i)
        data_to_graph[2][i] = u[1][i]

    np.savetxt("graphs_data/"+scheme.get_name()+"(graph).csv",
               data_to_graph, delimiter=",")

    plot_graph(data_to_graph, scheme.get_name())

    # print(u)


def main():
    sigma = 0.25
    N = 201
    X = 2
    h = X/(N-1)
    lam = 1
    tau = sigma*h/lam
    time_steps = 200
    print("distance = ", tau*time_steps)

    newpath = "graphs"
    if not path.exists(newpath):
        makedirs(newpath)

    newpath = "graphs_data"
    if not path.exists(newpath):
        makedirs(newpath)

    alpha_0_m1 = 0.25
    alpha_0_0 = 0.75
    alpha_m1_0 = 1/(1+2*sigma) - alpha_0_m1*sigma / \
        (1+2*sigma) - alpha_0_0*(1+sigma)/(1+2*sigma)
    alpha_1_m1 = (2*sigma)/(1+2*sigma) - alpha_0_m1*(1+sigma) / \
        (1+2*sigma) - alpha_0_0*(sigma)/(1+2*sigma)
    first_ord_right_up = Sceme("First_order_approximation_right_up",
                               alpha_0_m1=alpha_0_m1, alpha_0_0=alpha_0_0,
                               alpha_1_m1=alpha_1_m1, alpha_m1_0=alpha_m1_0)

    calculate_scheme(first_ord_right_up, N, h, tau, time_steps)

    alpha_0_m1 = 0.4
    alpha_0_0 = 0
    alpha_m1_0 = 1/(1+2*sigma) - alpha_0_m1*sigma / \
        (1+2*sigma) - alpha_0_0*(1+sigma)/(1+2*sigma)
    alpha_1_m1 = (2*sigma)/(1+2*sigma) - alpha_0_m1*(1+sigma) / \
        (1+2*sigma) - alpha_0_0*(sigma)/(1+2*sigma)
    first_ord_right_down = Sceme("First_order_approximation_right_down",
                                 alpha_0_m1=alpha_0_m1, alpha_0_0=alpha_0_0,
                                 alpha_1_m1=alpha_1_m1, alpha_m1_0=alpha_m1_0)

    calculate_scheme(first_ord_right_down, N, h, tau, time_steps)

    alpha_0_m1 = 0
    alpha_0_0 = 0.8
    alpha_m1_0 = 1/(1+2*sigma) - alpha_0_m1*sigma / \
        (1+2*sigma) - alpha_0_0*(1+sigma)/(1+2*sigma)
    alpha_1_m1 = (2*sigma)/(1+2*sigma) - alpha_0_m1*(1+sigma) / \
        (1+2*sigma) - alpha_0_0*(sigma)/(1+2*sigma)
    first_ord_left_up = Sceme("First_order_approximation_left_up",
                              alpha_0_m1=alpha_0_m1, alpha_0_0=alpha_0_0,
                              alpha_1_m1=alpha_1_m1, alpha_m1_0=alpha_m1_0)

    calculate_scheme(first_ord_left_up, N, h, tau, time_steps)

    alpha_0_m1 = 0
    alpha_0_0 = 0
    alpha_m1_0 = 1/(1+2*sigma) - alpha_0_m1*sigma / \
        (1+2*sigma) - alpha_0_0*(1+sigma)/(1+2*sigma)
    alpha_1_m1 = (2*sigma)/(1+2*sigma) - alpha_0_m1*(1+sigma) / \
        (1+2*sigma) - alpha_0_0*(sigma)/(1+2*sigma)
    first_ord_left_down = Sceme("First_order_approximation_left_down",
                                alpha_0_m1=alpha_0_m1, alpha_0_0=alpha_0_0,
                                alpha_1_m1=alpha_1_m1, alpha_m1_0=alpha_m1_0)

    calculate_scheme(first_ord_left_down, N, h, tau, time_steps)

    # ----------2 порядок по фридрихсу------------------------------------
    alpha_0_0 = 1.05
    alpha_0_m1 = alpha_0_0*(-sigma-1)/(sigma+1) + 2/(sigma+1)
    alpha_m1_0 = alpha_0_0*(-sigma-1)/((sigma+1)*(2*sigma+1)) + \
        (-sigma+1)/((sigma+1)*(2*sigma+1))
    alpha_1_m1 = 1/(2*sigma + 1)*alpha_0_0 + (2*sigma - 2)/(2*sigma+1)
    second_ord_friedr_closest = Sceme("Second_order_approximation_friedr_closest",
                                      alpha_0_m1=alpha_0_m1, alpha_0_0=alpha_0_0,
                                      alpha_1_m1=alpha_1_m1, alpha_m1_0=alpha_m1_0)

    calculate_scheme(second_ord_friedr_closest, N, h, tau, time_steps)

    # ----------2 порядок ещё 2 схемы------------------------------------
    alpha_0_0 = 1.5
    alpha_0_m1 = alpha_0_0*(-sigma-1)/(sigma+1) + 2/(sigma+1)
    alpha_m1_0 = alpha_0_0*(-sigma-1)/((sigma+1)*(2*sigma+1)) + \
        (-sigma+1)/((sigma+1)*(2*sigma+1))
    alpha_1_m1 = 1/(2*sigma + 1)*alpha_0_0 + (2*sigma - 2)/(2*sigma+1)
    second_ord_1 = Sceme("Second_order_approximation_1",
                         alpha_0_m1=alpha_0_m1, alpha_0_0=alpha_0_0,
                         alpha_1_m1=alpha_1_m1, alpha_m1_0=alpha_m1_0)

    calculate_scheme(second_ord_1, N, h, tau, time_steps)

    alpha_0_0 = 0.6
    alpha_0_m1 = alpha_0_0*(-sigma-1)/(sigma+1) + 2/(sigma+1)
    alpha_m1_0 = alpha_0_0*(-sigma-1)/((sigma+1)*(2*sigma+1)) + \
        (-sigma+1)/((sigma+1)*(2*sigma+1))
    alpha_1_m1 = 1/(2*sigma + 1)*alpha_0_0 + (2*sigma - 2)/(2*sigma+1)
    second_ord_2 = Sceme("Second_order_approximation_2",
                         alpha_0_m1=alpha_0_m1, alpha_0_0=alpha_0_0,
                         alpha_1_m1=alpha_1_m1, alpha_m1_0=alpha_m1_0)

    calculate_scheme(second_ord_2, N, h, tau, time_steps)

    # ----------3 порядок ------------------------------------------------
    alpha_0_m1 = (2*sigma)/(sigma + 1)
    alpha_0_0 = -2*(sigma-1)/(sigma+1)
    alpha_m1_0 = (sigma-1)/(2*sigma*sigma + 3*sigma + 1)
    alpha_1_m1 = (2*sigma*(sigma-1))/((sigma+1)*(2*sigma+1))
    third_ord = Sceme("Third_order_approximation",
                      alpha_0_m1=alpha_0_m1, alpha_0_0=alpha_0_0,
                      alpha_1_m1=alpha_1_m1, alpha_m1_0=alpha_m1_0)

    calculate_scheme(third_ord, N, h, tau, time_steps)


if __name__ == '__main__':
    main()
