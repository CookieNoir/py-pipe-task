import pipe_task
from decimal import Decimal

# это та функция, которая у каждого должна быть своя
def func01(t, a):
    # при t∈[0, 1] результат на отрезке [0, 1]
    return (t * t * (Decimal(3) - Decimal(2) * t)) ** a

# число шагов в методе Рунге-Кутты
steps = 5000
# начальные параметры метода деления отрезка пополам
left_border = Decimal(0.01)
right_border = Decimal(100)
threshold = Decimal(0.025)

pipe_task.solve_task(steps, left_border, right_border, threshold, func01)