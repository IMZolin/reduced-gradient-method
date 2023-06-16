from one_d_min_lib.uniform_search import uniform_search_method
from one_d_min_lib.Trial_Point_Method_file import trial_point_method
from one_d_min_lib.golden_egg import golden_search


def target_function(x_):
    """целевая функция"""
    return 10 * (((x_ - 1) ** 2) ** (1 / 3)) / (x_ ** 2 + 9)


def fun(self, x_):
    return x_ + 2


class One_D_Problem:
    def __init__(self, t0=-3, t1=3, targ_fun=lambda a: target_function(a)):
        self.left_border = t0  # интервал
        self.right_border = t1
        self.target_function = targ_fun

    uniform_search_method = uniform_search_method
    trial_point_method = trial_point_method
    golden_search = golden_search

    """
    Здесь добавляем методы решения данной задачи
    Методы возвращают два числа: 1) ответ 2) количество обращений к вычислению функции
    Метод равномерного поиска: (точность accuracy, число разбиений n) -> (ответ, то есть любая точка последнего интервала)
    Метод золотого сечения: (точность accuracy) -> (ответ, то есть любая точка последнего интервала)
    Метод пробных точек: (точность accuracy, ...) -> (ответ, то есть любая точка последнего интервала)
    """
    """
    Для каждого метода надо написать метод красивой печати
    """
