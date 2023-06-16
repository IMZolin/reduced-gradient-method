from itertools import combinations
import numpy as np
from one_d_min_lib.One_D_Problem_file import *


class Target_function:
    def __init__(self, f_, grad_):
        self.f = f_
        self.grad = grad_


class Wolf:
    def __init__(self, f_: Target_function, A_, b_):
        self.f0 = f_
        self.A = A_
        self.b = b_
        self.x = []
        self.I = []
        self.r = []
        self.d = [0.0 for _ in self.A[0]]
        self.lmd = None
        self.n = len(self.A[0])
        self.m = len(self.A)
        self.B_ = []  # B^-1

    def find_x0(self):
        x0 = [0.0 for _ in range(self.n)]

        comb = list(combinations(range(len(x0)), self.m))  # варианты (индексы столбцов) квадратных матриц
        for item in comb:  # пытаемся решить квадратную матрицу. Как только решили => stop
            square_matrix = np.matrix(self.A)[:, list(item)]
            try:
                ans = np.linalg.solve(square_matrix, self.b)
                if all([ans[j] >= 0 for j in range(len(ans))]):
                    for i in range(len(item)):
                        x0[item[i]] = ans[i]
                    break
            except Exception as ex:
                print(ex)
                pass
        self.x = x0

    def I_upd(self):
        if self.I:
            I_prev = self.I
        self.I = [0]
        min_element_index = 0
        for i in range(1, self.n):  # проходимся по всем компонентам, но выбираем только m наибольших
            if len(self.I) < self.m:
                self.I.append(i)
                if self.x[i] < self.x[min_element_index]:
                    min_element_index = i
            elif self.x[i] > self.x[min_element_index]:
                self.I[min_element_index] = i
            # elif i not in self.I and any([self.x[i] == self.x[j]] for j in self.I):


    def r_upd(self):
        B = np.matrix(self.A)[:, self.I]
        B_ = np.linalg.inv(B)
        self.B_ = B_

        tmp = [self.f0.grad(self.x)[i] for i in self.I]
        B_A = np.matmul(B_, self.A)

        grad_ = self.f0.grad(self.x)
        grad_b = np.dot(tmp, B_A).tolist()[0]

        self.r = [grad_[i] - grad_b[i] for i in range(self.n)]

    def d_upd(self):
        d_n = []
        I_n = [i for i in range(self.n) if i not in self.I]
        for i in I_n:
            if self.r[i] <= 0:
                d_n.append(-self.r[i])
            else:
                d_n.append(-self.r[i] * self.x[i])

        N = np.matrix(self.A)[:, I_n]
        B_N = np.matmul(self.B_, N)

        d_b = -np.dot(B_N, d_n)
        d_b = d_b.tolist()[0]

        for i in range(len(self.I)):
            self.d[self.I[i]] = d_b[i]
        for i in range(len(I_n)):
            self.d[I_n[i]] = d_n[i]

    def lmd_upd(self, eps=0.001):
        lmd_max_list = [-self.x[i] / self.d[i] for i in range(self.n) if self.d[i] < 0]
        if len(lmd_max_list) != 0:
            lmd_max = min(lmd_max_list)
        else:
            lmd_max = 100  # большое число

        pr = One_D_Problem(0, lmd_max, lambda t: self.f0.f([self.x[i] + t * self.d[i] for i in range(self.n)]))
        self.lmd = pr.golden_search(eps)[0]


    def x_upd(self):
        self.x = [self.x[i] + self.lmd * self.d[i] for i in range(self.n)]

    def minimize(self, eps=0.01):
        self.find_x0()
        x_array = [self.x]
        for i in range(30):
            self.I_upd()
            self.r_upd()
            self.d_upd()
            self.lmd_upd(eps)
            print('N', i)
            self.print_params()
            self.x_upd()

            x_array.append(self.x)

            if all([abs(self.d[i]) < eps for i in range(self.n)]):
                print('STOP: d = 0')
                break

        return x_array

    def print_params(self):
        print('I', self.I)
        print('x', [round(self.x[i], 5) for i in range(self.n)])
        print('f', round(self.f0.f(self.x), 5))
        print('d', [round(self.d[i], 5) for i in range(self.n)])
        print('lmd', self.lmd)
        print()


if __name__ == '__main__':

    # A1 = [[1, 1, 1, 0],
    #       [1, -1, 0, 1]]
    # b1 = [5, -2]
    # f1 = Target_function(lambda x: (x[0] - 3)**2 + (x[1] - 1)**2,
    #                      lambda x: [2*(x[0] - 3), 2*(x[1] - 1), 0, 0])
    # w1 = Wolf(f1, A1, b1)
    # w1.minimize()

    A2 = [[-0.5, -1, 1, 0, 0],
          [-1.5, -1, 0, 1, 0],
          [1, -1, 0, 0, 1]]
    b2 = [-2, -3, 4]
    f2 = Target_function(lambda x: (x[0] - 8)**2 + (x[1] + 2)**2,
                         lambda x: [2*(x[0] - 8), 2*(x[1] + 2), 0, 0, 0])
    w2 = Wolf(f2, A2, b2)
    w2.minimize()

