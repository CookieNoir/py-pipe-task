from runge_kutta_4 import rungeKutta
import matplotlib.pyplot as plt
from decimal import Decimal

u = [Decimal(15.19),        #1
    Decimal(8.18),          #2
    Decimal(13.198),        #3
    Decimal(3.543),         #4
    Decimal(4723.7),        #5
    Decimal(423.7),         #6
    Decimal(204.41),        #7
    Decimal(0.001466),      #8
    Decimal(0.013),         #9
    Decimal(0.09),          #10
    Decimal(0.000005428),   #11
    Decimal(0.024),         #12
    Decimal(0.00000592)]    #13
E = [Decimal(25000),        #1
    Decimal(25000),         #2
    Decimal(25000),         #3
    Decimal(25000),         #4
    Decimal(40000),         #5
    Decimal(40000),         #6
    Decimal(40000),         #7
    Decimal(20000),         #8
    Decimal(20000),         #9
    Decimal(20000),         #10
    Decimal(20000),         #11
    Decimal(20000),         #12
    Decimal(20000)]         #13
m = [Decimal(18),           #0
    Decimal(84),            #1
    Decimal(56),            #2
    Decimal(42),            #3
    Decimal(28),            #4
    Decimal(92),            #5
    Decimal(16)]            #6
q = [Decimal(78),           #2
     Decimal(140),          #3
     Decimal(140)]          #4
def p(l):
    return Decimal(5) - l / Decimal(60)

def R(i, x7):
    deg = Decimal(23) - E[i - 1] / x7
    result = u[i - 1] * deg.exp()
    #print("R{0} = {1}".format(i, result))
    return result

func01 = lambda t, a: t
alpha = Decimal(1)
def x7(l):
    return Decimal(373) + Decimal(1127) * (func01(l / Decimal(180), alpha))

def v(l, x):
    G = Decimal(1750)
    Gs = Decimal(3500)
    sum = Decimal(0)
    for i in range(0, 6):
        sum += m[i + 1] * x[i]
    result = Decimal(509.209) * p(l) * (sum / (x7(l) * (G + Gs * sum / m[0])))
    return result

def task(l, x):
    x7calc = x7(l)
    vl = v(l, x)
    R1 = R(1, x7calc)
    R2 = R(2, x7calc)
    R3 = R(3, x7calc)
    R4 = R(4, x7calc)
    R5 = R(5, x7calc)
    R6 = R(6, x7calc)
    R7 = R(7, x7calc)
    R8 = R(8, x7calc)
    R9 = R(9, x7calc)
    R10 = R(10, x7calc)
    R11 = R(11, x7calc)
    R12 = R(12, x7calc)
    R13 = R(13, x7calc)
    dx1 = Decimal(-1) * (R1 + R2 + R3 + R4) * x[0] * vl
    dx2 = (R3 * x[0] - (R6 + R7 + R10 + R13) * x[1]) * vl
    dx3 = (R2 * x[0] + R6 * x[1] - (R5 + R9 + R12) * x[2]) * vl
    dx4 = (R1 * x[0] + R7 * x[1] + R5 * x[2] - (R8 + R11) * x[3]) * vl
    dx5 = (R10 * x[1] + R9 * x[2] + R8 * x[3]) * vl
    dx6 = (R4 * x[0] + R13 * x[1] + R12 * x[2] + R11 * x[3]) * vl
    return [dx1, dx2, dx3, dx4, dx5, dx6]

def quality(x):
    size = len(x)
    result = [Decimal(0)] * size
    for i in range(0, size):
        uppersum = Decimal(0)
        for up in range(0, 3):
            uppersum += q[up] * m[2 + up] * x[i][1 + up]
        lowersum = Decimal(0)
        for down in range(0, 6):
            lowersum += m[1 + down] * x[i][down]
        result[i] = uppersum / lowersum
    return result

l0 = Decimal(0)
l1 = Decimal(180)
x0 = [Decimal(1),
    Decimal(0),
    Decimal(0),
    Decimal(0),
    Decimal(0),
    Decimal(0)]
n = 5000
def iteration(value):
    global alpha
    alpha = value
    rungeResult = rungeKutta(l0, x0, l1, n, task)
    x_quality = quality(rungeResult[1])
    result = max(x_quality)
    return result

# returns best alpha for pipe task
def bisection(left_border, right_border, threshold):
    a = left_border
    b = right_border
    mid = (a + b) / Decimal(2)
    fmid = iteration(mid)
    delta = abs(b - a)
    while delta > threshold:
        y = a + delta / Decimal(4)
        fy = iteration(y)
        z = b - delta / Decimal(4)
        fz = iteration(z)
        if fy > fmid:
            b = mid
            mid = y
            fmid = fy
        else:
            if fz > fmid:
                a = mid
                mid = z
                fmid = fz
            else:
                a = y
                b = z
        delta = abs(b - a)
        print(f"new segment is [{a}, {b}], f(y) is {fy}, f(z) is {fz}")
    return (a + b) / Decimal(2)

def solve_task(steps, left_border, right_border, threshold, custom_func01):
    global n, func01, alpha
    n = steps
    func01 = custom_func01
    alpha = bisection(left_border, right_border, threshold)
    rungeResult = rungeKutta(l0, x0, l1, n, task)
    x_quality = quality(rungeResult[1])
    print(f"max quality is {max(x_quality)}")

    # sum = 1 ?
    res = Decimal(0)
    last = len(rungeResult[1]) - 1
    for i in range(0, len(rungeResult[1][last])):
        res += rungeResult[1][last][i]
    print(f"sum of x[i] is {res}")

    fiqure,axis = plt.subplots(2)
    axis[0].plot(rungeResult[0], rungeResult[1])
    axis[0].legend(["x1","x2","x3","x4","x5","x6"])
    axis[0].set_xlim([l0, l1])
    axis[0].set_title("x(l)")

    axis[1].plot(rungeResult[0], x_quality)
    axis[1].set_xlim([l0, l1])
    axis[1].set_title("quality")

    plt.show()
