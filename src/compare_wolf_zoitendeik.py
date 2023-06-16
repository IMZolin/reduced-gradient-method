from zoitendeik_lib.Zoitendeik import *
from Wolf import Wolf
import matplotlib.pyplot as plt

if __name__ == '__main__':
    phi0 = Target_function(lambda x: (x[0] - 8) ** 2 + (x[1] + 2) ** 2,
                           lambda x: [2 * (x[0] - 8), 2 * (x[1] + 2)])

    phi1 = Constraint('ineq',
                      lambda x: -0.5*x[0] - x[1] + 2,
                      lambda x: [-0.5, -1])

    phi2 = Constraint('ineq',
                      lambda x: -1.5*x[0] - x[1] + 3,
                      lambda x: [-1.5, -1])

    phi3 = Constraint('ineq',
                      lambda x: -x[0],
                      lambda x: [-1, 0])

    phi4 = Constraint('ineq',
                      lambda x: -x[1],
                      lambda x: [0, -1])

    phi5 = Constraint('ineq',
                      lambda x: x[0] - x[1] - 4,
                      lambda x: [1, -1])

    z = Zoitendeik_step(phi0, [phi1, phi2, phi3, phi4], [0.0, 0.0], 0.25, 0.5)

    # K, R изменены в Zoitendeik.py

    print('ZOITENDEIK\n')
    ztd_arr = z.minimize(eps=0.01)

    print('\nWOLF\n')
    A2 = [[-0.5, -1, 1, 0],
          [-1.5, -1, 0, 1]]
    b2 = [-2, -3]
    f2 = Target_function(lambda x: (x[0] - 8) ** 2 + (x[1] + 2) ** 2,
                         lambda x: [2 * (x[0] - 8), 2 * (x[1] + 2), 0, 0])
    w2 = Wolf(f2, A2, b2)
    wlf_arr = w2.minimize(eps=0.01)

    x_star = [8, 0]
    plt.plot(range(len(ztd_arr)),
             [norma_calculate(np.array(ztd_arr[i]) - np.array(x_star)) for i in range(len(ztd_arr))],
             '--',
             label='Zoitendeik')
    plt.plot(range(len(wlf_arr)),
             [norma_calculate(np.array(wlf_arr[i])[:2] - np.array(x_star)) for i in range(len(wlf_arr))],
             '-o',
             label='Wolf')
    plt.xlabel('iteration')
    plt.ylabel('||x_k - x*||')
    plt.legend()
    plt.semilogy()

    plt.show()



