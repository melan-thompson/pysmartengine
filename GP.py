import numpy as np


def GaussianDistribution(X, mean=None, covariance=None):
    import numpy as np
    n = len(X)
    if mean is None:
        mean = np.zeros(n)
    if covariance is None:
        covariance = np.eye(n)
    # ev,evec=np.linalg.eig(covariance)
    # print(evec)
    temp = np.dot(np.dot(np.transpose(X - mean), np.linalg.inv(covariance)), X - mean)
    detC = np.linalg.det(covariance)
    return 1 / pow(2 * np.pi, n / 2.) / pow(detC, 1. / 2.) * np.exp(-1. / 2. * temp)


def polt2DGaussianDistribution(plotrange, mean=None, cov=None):
    x = y = np.arange(plotrange[0], plotrange[1], 0.1)
    xx, yy = np.meshgrid(x, y)
    # cov = np.array([[0.9, 0.2], [0.2, 0.3]])
    cov = np.eye(2)
    shape = xx.shape
    zz = np.zeros(shape)
    for i in range(shape[0]):
        for j in range(shape[1]):
            temp = [xx[i][j], yy[i][j]]
            zz[i][j] = GaussianDistribution(temp, covariance=cov)
    import matplotlib.pyplot as plt
    # from matplotlib import cm
    fig = plt.figure(figsize=(8, 8))
    # ax = fig.gca(projection='3d')
    # surf=ax.plot_surface(xx,yy,zz,cmap="Blues",
    #                    linewidth=0, antialiased=False)
    # plt.show()
    ev, evec = np.linalg.eig(cov)
    evecT = evec.transpose()
    unitevec1 = evecT[0] / np.linalg.norm(evecT[0], ord=2)
    unitevec2 = evecT[1] / np.linalg.norm(evecT[1], ord=2)
    vector1 = unitevec1 * ev[0]
    vector2 = unitevec2 * ev[1]
    print(unitevec1, unitevec2)

    plt.contourf(xx, yy, zz, 10, cmap="Blues")
    plt.plot([0, vector1[0]], [0, vector1[1]], 'ro-')
    plt.plot([0, vector2[0]], [0, vector2[1]], 'ro-')
    plt.tight_layout()
    plt.show()


def gaussKenel(delta_x, thita=1, pk=1):
    from math import exp
    return exp(-thita * pow(abs(delta_x), pk))


def threeTimes(delta_x, thita=0.1):
    temp = thita * abs(delta_x)
    if 0 <= temp <= 0.2:
        return 1 - 15. * temp ** 2 + 30. * temp ** 3
    elif 0.2 < temp < 1:
        return 1.25 * (1 - temp) ** 3
    else:
        return 0

def R(theta,p,x1,x2):
    pass

