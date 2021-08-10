import numpy as np
import matplotlib.pyplot as plt
import copy


def hermite(arr):
    return np.conjugate(arr.T)


def plot_transmissions(agf_list, pic_name="agf_test.png"):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlabel("Frequency $\omega$")
    plt.ylabel("Transmission")
    omega_max = np.max([agf_list[i].omega_max for i in range(len(agf_list))])

    plt.xlim(-0.01*omega_max, omega_max)
    plt.ylim(0.0, 1.01)
    for i in range(len(agf_list)):
        plt.plot(agf_list[i].omega, agf_list[i].trans_list,
                 label=agf_list[i].legend)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.savefig(pic_name, dpi=100, bbox_inches='tight')


class AGF_one_d:
    def __init__(self, mc=4.6*1e-26, md=4.6*1e-26, fc=32, fd=32, fcd=32, omega_max=6.0*1e13, omega_len=1000,
                 iterative_num=100, epsilon=1e-6, delta_0=1e-3, legend="Homogeneous", sgf_mode="decimation"):
        self.mc = mc
        self.md = md
        self.fc = fc
        self.fd = fd
        self.fcd = fcd
        self.omega_max = omega_max
        self.omega_len = omega_len
        self.iterative_num = iterative_num
        self.epsilon = epsilon
        self.delta_0 = delta_0
        self.legend = legend
        self.sgf_mode = sgf_mode

        self.h_d = np.array([[(self.fd+self.fcd)/self.md, -self.fd/self.md, 0],
                             [-self.fd/self.md, 2*self.fd /
                                 self.md, -self.fd/self.md],
                             [0, -self.fd/self.md, (self.fd+self.fcd)/self.md]])

    def delta(self, omega):
        return omega**2*(1-omega/self.omega_max)*self.delta_0

    def sgf_iterative(self, omega):
        g_00 = 1/(omega**2+1j*self.delta(omega)-2*self.fc/self.mc)
        g_s = g_00
        g_b = 0
        i = 0
        while np.linalg.norm(g_b - g_s)/np.linalg.norm(g_s) > self.epsilon and i < self.iterative_num:
            g_b = copy.deepcopy(g_s)
            g_s = 1/(omega**2+1j*self.delta(omega)-2 *
                     self.fc/self.mc-g_s*(self.fc/self.mc)**2)
            i += 1
        return g_s

    def sgf_decimation(self, omega):
        h_s0 = omega**2+1j*self.delta(omega)-(self.fc+self.fcd)/self.mc
        h_b0 = omega**2+1j*self.delta(omega)-2*self.fc/self.mc
        tau_0 = self.fc/self.mc
        h_sn = h_s0
        h_sn1 = 0
        h_bn = h_b0
        h_bn1 = 0
        tau_n = tau_0
        tau_n1 = 0
        i = 0
        while np.linalg.norm(tau_n/tau_0) > self.epsilon and i < self.iterative_num:
            h_sn1 = copy.deepcopy(h_sn)
            h_bn1 = copy.deepcopy(h_bn)
            tau_n1 = copy.deepcopy(tau_n)
            h_sn -= tau_n1**2/h_bn1
            h_bn -= 2*tau_n1**2/h_bn1
            tau_n = -tau_n1**2/h_bn1
            i += 1
        return 1/h_sn

    def sgf_analytic(self, omega):
        return (-(2*self.fc/self.mc - (omega**2 + 1j*self.delta(omega)))-np.sqrt((2*self.fc/self.mc
                                                                                  - (omega**2 + 1j*self.delta(omega)))**2 - 4*(self.fc/self.mc)**2))/2/(self.fc/self.mc)**2

    def sigma(self, omega):
        if self.sgf_mode == "decimation":
            g_1s = self.sgf_decimation(omega)
        elif self.sgf_mode == "iterative":
            g_1s = self.sgf_iterative(omega)
        else:
            g_1s = self.sgf_analytic(omega)
        sigma_1 = np.array(
            [[self.fcd*self.fcd*g_1s/self.mc/self.md, 0, 0], [0, 0, 0], [0, 0, 0]])
        sigma_2 = np.array(
            [[0, 0, 0], [0, 0, 0], [0, 0, self.fcd*self.fcd*g_1s/self.mc/self.md]])

        return sigma_1, sigma_2

    def gamma(self, omega):
        sigma_1, sigma_2 = self.sigma(omega)
        return 1j*(sigma_1 - hermite(sigma_1)), 1j*(sigma_2 - hermite(sigma_2))

    def trans_sub(self, omega):
        gamma_1, gamma_2 = self.gamma(omega)
        sigma_1, sigma_2 = self.sigma(omega)
        g_d = np.linalg.pinv((omega**2)*np.eye(3)-self.h_d-sigma_1-sigma_2)
        #g_d = np.linalg.inv((omega**2)*np.eye(3)-h_d-sigma_1-sigma_2)
        return np.trace(np.dot(gamma_1, np.dot(g_d, np.dot(gamma_2, hermite(g_d)))))

    def transmission(self):
        self.omega = np.linspace(0, self.omega_max, self.omega_len)
        self.trans_list = np.zeros(self.omega_len)
        for i in range(self.omega_len):
            self.trans_list[i] = self.trans_sub(self.omega[i])
            if i % 100 == 0:
                print(i)

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.xlabel("Frequency $\omega$")
        plt.ylabel("Transmission")
        plt.xlim(-0.01*self.omega_max, self.omega_max)
        plt.ylim(0.0, 1.01)
        plt.plot(self.omega, self.trans_list, label=self.legend)
        ax.legend()
