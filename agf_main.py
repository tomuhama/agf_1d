import numpy as np
import matplotlib.pyplot as plt


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
    ax.legend()

    plt.savefig(pic_name, dpi=100)


class AGF_one_d:
    def __init__(self, mc=4.6*1e-26, md=4.6*1e-26, fc=32, fd=32, fcd=32, omega_max=6.0*1e13, omega_len=1000, legend="Homogeneous"):
        self.mc = mc
        self.md = md
        self.fc = fc
        self.fd = fd
        self.fcd = fcd
        self.omega_max = omega_max
        self.omega_len = omega_len
        self.legend = legend

        self.h_d = np.array([[(self.fd+self.fcd)/self.md, -self.fd/self.md, 0],
                             [-self.fd/self.md, 2*self.fd /
                                 self.md, -self.fd/self.md],
                             [0, -self.fd/self.md, (self.fd+self.fcd)/self.md]])

    def delta(self, omega):
        return 1e-7*omega**2

    def sigma(self, omega):
        g_1s = (-(2*self.fc/self.mc - (omega**2 + 1j*self.delta(omega)))-np.sqrt((2*self.fc/self.mc
                                                                                  - (omega**2 + 1j*self.delta(omega)))**2 - 4*(self.fc/self.mc)**2))/2/(self.fc/self.mc)**2
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

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.xlabel("Frequency $\omega$")
        plt.ylabel("Transmission")
        plt.xlim(-0.01*self.omega_max, self.omega_max)
        plt.ylim(0.0, 1.01)
        plt.plot(self.omega, self.trans_list, label=self.legend)
        ax.legend()
