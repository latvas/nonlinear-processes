from math import sqrt
from os import makedirs, path
import numpy as np
import matplotlib.pyplot as plt

rho_0 = 0.25
c_0 = 2.0
u_L = 1.0
u_R = 0.0
p_L = 5.0
p_R = 2.0

y_L = u_L + (p_L)/(rho_0*c_0)
y_R = u_R + (p_R)/(rho_0*c_0)

z_L = u_L - (p_L)/(rho_0*c_0)
z_R = u_R - (p_R)/(rho_0*c_0)


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
def y_0(x):
    if x <= 0:
        return y_L
    else:
        return y_R


def z_0(x):
    if x <= 0:
        return z_L
    else:
        return z_R


def u_0(x):
    if x <= 0:
        return u_L
    else:
        return u_R


def p_0(x):
    if x <= 0:
        return p_L
    else:
        return p_R


def translate(y, z):
    shape = y.shape
    u = np.zeros(shape, dtype=np.double)
    p = np.zeros(shape, dtype=np.double)
    for i in range(shape[0]):
        u[i] = (y[i]+z[i])/2
        p[i] = rho_0*c_0*(y[i]-z[i])/2
    return u, p


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


def plot_difference_graph(data_to_graph, graph_name):
    N = list()
    X = list()
    Y = list()

    for i in range(data_to_graph.shape[1]):
        X.append(data_to_graph[0][i])
        Y.append(data_to_graph[3][i] - data_to_graph[2][i])
        N.append("difference")

    fig, ax = plt.subplots()

    for n in set(N):
        x1 = [i for i, j in zip(X, N) if j == n]
        y1 = [i for i, j in zip(Y, N) if j == n]
        ax.plot(x1, y1, label=n)

    plt.xlabel('x')
    plt.ylabel('f(x)')

    plt.grid(visible=True)
    plt.legend()
    plt.savefig("graphs/difference_"+graph_name+".png", dpi=300)


def calculate_scheme(scheme, N, h, tau, time_steps, lam):
    y = np.zeros((3, N), dtype=np.double)
    z = np.zeros((3, N), dtype=np.double)

    for i in range(N):
        # задаём начальные условия, предполагая, что -1 и 0 слой одинаковые
        y[0][i] = y_0(h*i-1)
        y[1][i] = y_0((h*i-1) - lam*tau)
        z[0][i] = z_0(h*i-1)
        z[1][i] = z_0((h*i-1) - lam*tau)

    y[2][0] = y_L
    y[2][N-1] = y_R

    z[2][0] = z_L
    z[2][N-1] = z_R

    weight = scheme.return_weights()
    for n in range(1, time_steps):
        for m in range(1, N-1):
            y[2][m] = weight[0]*y[2][m-1] + weight[1] * \
                y[1][m-1] + weight[2]*y[1][m] + weight[3]*y[0][m]

        y[[0, 1]] = y[[1, 0]]
        y[[1, 2]] = y[[2, 1]]
        y[2] = np.zeros(N)
        y[2][0] = y_L
        y[2][N-1] = y_R

    for n in range(1, time_steps):
        for m in range(N-2, 0, -1):
            z[2][m] = weight[0]*z[2][m+1] + weight[1] * \
                z[1][m+1] + weight[2]*z[1][m] + weight[3]*z[0][m]

        z[[0, 1]] = z[[1, 0]]
        z[[1, 2]] = z[[2, 1]]
        z[2] = np.zeros(N)
        z[2][0] = z_L
        z[2][N-1] = z_R

    # np.savetxt(scheme.name()+".csv", u, delimiter=",")
    '''
    data_to_graph = np.zeros((4, N), dtype=np.double)
    print(scheme.get_name(), weight)
    for i in range(N):
        data_to_graph[0][i] = i*h - 1
        data_to_graph[1][i] = y_0(h*i-1)
        data_to_graph[2][i] = y[1][i]
        data_to_graph[3][i] = y_0((h*i-1) - tau*time_steps)

    np.savetxt("graphs_data/y_"+scheme.get_name()+"(graph).csv",
               data_to_graph, delimiter=",")

    plot_graph(data_to_graph, "y_" + scheme.get_name())

    data_to_graph = np.zeros((4, N), dtype=np.double)
    print(scheme.get_name(), weight)
    for i in range(N):
        data_to_graph[0][i] = i*h - 1
        data_to_graph[1][i] = z_0(h*i-1)
        data_to_graph[2][i] = z[1][i]
        data_to_graph[3][i] = z_0((h*i-1) - tau*time_steps)

    np.savetxt("graphs_data/z_"+scheme.get_name()+"(graph).csv",
               data_to_graph, delimiter=",")

    plot_graph(data_to_graph, "z_" + scheme.get_name())
    '''

    u, p = translate(y[1], z[1])

    data_to_graph = np.zeros((4, N), dtype=np.double)
    print(scheme.get_name(), weight)
    for i in range(N):
        data_to_graph[0][i] = i*h - 1
        data_to_graph[1][i] = u_0(h*i-1)
        data_to_graph[2][i] = u[i]
        data_to_graph[3][i] = u_0((h*i-1) - tau*time_steps)

    np.savetxt("graphs_data/u_"+scheme.get_name()+"(graph).csv",
               data_to_graph, delimiter=",")

    plot_graph(data_to_graph, "u_" + scheme.get_name())

    data_to_graph = np.zeros((4, N), dtype=np.double)
    print(scheme.get_name(), weight)
    for i in range(N):
        data_to_graph[0][i] = i*h - 1
        data_to_graph[1][i] = p_0(h*i-1)
        data_to_graph[2][i] = p[i]
        data_to_graph[3][i] = p_0((h*i-1) - tau*time_steps)

    np.savetxt("graphs_data/p_"+scheme.get_name()+"(graph).csv",
               data_to_graph, delimiter=",")

    plot_graph(data_to_graph, "p_" + scheme.get_name())

    # print(u)


def calculate_hybrid_scheme(schemes, N, h, tau, time_steps, lam):
    y = np.zeros((3, N), dtype=np.double)
    z = np.zeros((3, N), dtype=np.double)

    for i in range(N):
        # задаём начальные условия, предполагая, что -1 и 0 слой одинаковые
        y[0][i] = y_0(h*i-1)
        y[1][i] = y_0((h*i-1) - lam*tau)
        z[0][i] = z_0(h*i-1)
        z[1][i] = z_0((h*i-1) - lam*tau)

    y[2][0] = y_L
    y[2][N-1] = y_R

    z[2][0] = z_L
    z[2][N-1] = z_R

    for n in range(1, time_steps):
        for m in range(1, N):
            for scheme in schemes:
                weight = scheme.return_weights()
                y[2][m] = weight[0]*y[2][m-1] + weight[1] * \
                    y[1][m-1] + weight[2]*y[1][m] + weight[3]*y[0][m]
                delta = y[1][m-1] - y[1][m]
                if (delta < 0 and y[1][m-1] < y[2][m] < y[1][m]) or (delta > 0 and y[1][m-1] > y[2][m] > y[1][m]):
                    break

        y[[0, 1]] = y[[1, 0]]
        y[[1, 2]] = y[[2, 1]]
        y[2] = np.zeros(N)
        y[2][0] = y_L
        y[2][N-1] = y_R

    for n in range(1, time_steps):
        for m in range(N-2, 0, -1):
            z[2][m] = weight[0]*z[2][m+1] + weight[1] * \
                z[1][m+1] + weight[2]*z[1][m] + weight[3]*z[0][m]

    for n in range(1, time_steps):
        for m in range(N-2, 0, -1):
            for scheme in schemes:
                weight = scheme.return_weights()
                z[2][m] = weight[0]*z[2][m+1] + weight[1] * \
                    z[1][m+1] + weight[2]*z[1][m] + weight[3]*z[0][m]
                delta = z[1][m+1] - z[1][m]
                if (delta < 0 and z[1][m+1] < z[2][m] < z[1][m]) or (delta > 0 and z[1][m+1] > z[2][m] > z[1][m]):
                    break

        z[[0, 1]] = z[[1, 0]]
        z[[1, 2]] = z[[2, 1]]
        z[2] = np.zeros(N)
        z[2][0] = z_L
        z[2][N-1] = z_R

    # np.savetxt(scheme.name()+".csv", u, delimiter=",")
    '''
    data_to_graph = np.zeros((4, N), dtype=np.double)
    print(scheme.get_name(), weight)
    for i in range(N):
        data_to_graph[0][i] = i*h - 1
        data_to_graph[1][i] = y_0(h*i-1)
        data_to_graph[2][i] = y[1][i]
        data_to_graph[3][i] = y_0((h*i-1) - tau*time_steps)

    np.savetxt("graphs_data/y_"+scheme.get_name()+"(graph).csv",
               data_to_graph, delimiter=",")

    plot_graph(data_to_graph, "y_" + scheme.get_name())

    data_to_graph = np.zeros((4, N), dtype=np.double)
    print(scheme.get_name(), weight)
    for i in range(N):
        data_to_graph[0][i] = i*h - 1
        data_to_graph[1][i] = z_0(h*i-1)
        data_to_graph[2][i] = z[1][i]
        data_to_graph[3][i] = z_0((h*i-1) - tau*time_steps)

    np.savetxt("graphs_data/z_"+scheme.get_name()+"(graph).csv",
               data_to_graph, delimiter=",")

    plot_graph(data_to_graph, "z_" + scheme.get_name())
    '''

    y, p = translate(y[1], z[1])

    data_to_graph = np.zeros((4, N), dtype=np.double)

    print("--------------Hybrid scheme----------------")
    for scheme in schemes:
        print(scheme.get_name(), weight)
    print("-------------------------------------------")

    name = str()
    for scheme in schemes:
        name += scheme.get_name()

    for i in range(N):
        data_to_graph[0][i] = i*h - 1
        data_to_graph[1][i] = u_0(h*i-1)
        data_to_graph[2][i] = y[i]
        data_to_graph[3][i] = u_0((h*i-1) - tau*time_steps)

    np.savetxt("graphs_data/u_"+name+"(graph).csv",
               data_to_graph, delimiter=",")

    plot_graph(data_to_graph, "u_" + name)

    data_to_graph = np.zeros((4, N), dtype=np.double)
    print(scheme.get_name(), weight)
    for i in range(N):
        data_to_graph[0][i] = i*h - 1
        data_to_graph[1][i] = p_0(h*i-1)
        data_to_graph[2][i] = p[i]
        data_to_graph[3][i] = p_0((h*i-1) - tau*time_steps)

    np.savetxt("graphs_data/p_"+name+"(graph).csv",
               data_to_graph, delimiter=",")

    plot_graph(data_to_graph, "p_" + name)


def main():
    sigma = 0.25
    lam = c_0
    h = 0.01
    tau = sigma*h/lam
    time_steps = 100
    N = 201

    print("distance = ", lam*tau*time_steps)

    newpath = "graphs"
    if not path.exists(newpath):
        makedirs(newpath)

    newpath = "graphs_data"
    if not path.exists(newpath):
        makedirs(newpath)

    # ----------1 порядок с мин вязкостью----------------------------------
    alpha_0_m1 = 0.25
    alpha_0_0 = 0.75
    alpha_m1_0 = 1/(1+2*sigma) - alpha_0_m1*sigma / \
        (1+2*sigma) - alpha_0_0*(1+sigma)/(1+2*sigma)
    alpha_1_m1 = (2*sigma)/(1+2*sigma) - alpha_0_m1*(1+sigma) / \
        (1+2*sigma) - alpha_0_0*(sigma)/(1+2*sigma)
    first_ord_right_up = Sceme("First_order_approximation_right_up",
                               alpha_0_m1=alpha_0_m1, alpha_0_0=alpha_0_0,
                               alpha_1_m1=alpha_1_m1, alpha_m1_0=alpha_m1_0)

    calculate_scheme(first_ord_right_up, N, h, tau, time_steps, lam)

    # ----------2 порядок по фридрихсу------------------------------------
    alpha_0_0 = 1.05
    alpha_0_m1 = alpha_0_0*(-sigma-1)/(sigma+1) + 2/(sigma+1)
    alpha_m1_0 = alpha_0_0*(-sigma-1)/((sigma+1)*(2*sigma+1)) + \
        (-sigma+1)/((sigma+1)*(2*sigma+1))
    alpha_1_m1 = 1/(2*sigma + 1)*alpha_0_0 + (2*sigma - 2)/(2*sigma+1)
    second_ord_friedr_closest = Sceme("Second_order_approximation_friedr_closest",
                                      alpha_0_m1=alpha_0_m1, alpha_0_0=alpha_0_0,
                                      alpha_1_m1=alpha_1_m1, alpha_m1_0=alpha_m1_0)

    calculate_scheme(second_ord_friedr_closest, N, h, tau, time_steps, lam)

    # ----------3 порядок ------------------------------------------------
    alpha_0_m1 = (2*sigma)/(sigma + 1)
    alpha_0_0 = -2*(sigma-1)/(sigma+1)
    alpha_m1_0 = (sigma-1)/(2*sigma*sigma + 3*sigma + 1)
    alpha_1_m1 = (2*sigma*(sigma-1))/((sigma+1)*(2*sigma+1))
    third_ord = Sceme("Third_order_approximation",
                      alpha_0_m1=alpha_0_m1, alpha_0_0=alpha_0_0,
                      alpha_1_m1=alpha_1_m1, alpha_m1_0=alpha_m1_0)

    calculate_scheme(third_ord, N, h, tau, time_steps, lam)

    scemes_list = (third_ord, second_ord_friedr_closest, first_ord_right_up)
    calculate_hybrid_scheme(scemes_list, N, h, tau, time_steps, lam)
    '''
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

    alpha_0_0 = 0.8
    alpha_0_m1 = alpha_0_0*(-sigma-1)/(sigma+1) + 2/(sigma+1)
    alpha_m1_0 = alpha_0_0*(-sigma-1)/((sigma+1)*(2*sigma+1)) + \
        (-sigma+1)/((sigma+1)*(2*sigma+1))
    alpha_1_m1 = 1/(2*sigma + 1)*alpha_0_0 + (2*sigma - 2)/(2*sigma+1)
    second_ord_2 = Sceme("Second_order_approximation_2",
                         alpha_0_m1=alpha_0_m1, alpha_0_0=alpha_0_0,
                         alpha_1_m1=alpha_1_m1, alpha_m1_0=alpha_m1_0)

    calculate_scheme(second_ord_2, N, h, tau, time_steps)

    

    # ----------Считаем по гибридным схемам ------------------------------
    # 5п
    scemes_list = (second_ord_1, second_ord_2, first_ord_right_up)
    calculate_hybrid_scheme(scemes_list, N, h, tau, time_steps)
    # 6п
    scemes_list = (third_ord, second_ord_1, first_ord_right_up)
    calculate_hybrid_scheme(scemes_list, N, h, tau, time_steps)

    scemes_list = (third_ord, second_ord_2, first_ord_right_up)
    calculate_hybrid_scheme(scemes_list, N, h, tau, time_steps)
    # 7п
    scemes_list = (third_ord, second_ord_1, second_ord_2, first_ord_right_up)
    calculate_hybrid_scheme(scemes_list, N, h, tau, time_steps)
    '''


if __name__ == '__main__':
    main()
