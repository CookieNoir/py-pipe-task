from decimal import Decimal

def scalarVectorMultiplication(scalar, vector):
    size = len(vector)
    result = [Decimal(0)] * size
    for i in range(size):
        result[i] = vector[i] * scalar
    return result

def vectorSum(vector1, vector2):
    size = len(vector1)
    result = [Decimal(0)] * size
    for i in range(size):
        result[i] = vector1[i] + vector2[i]
    return result

def rungeKutta(x0, y0, x1, n, dydx):
    h = (x1 - x0) / n
    ycount = len(y0)

    x = [Decimal(0)] * (n + 1)
    x[0] = x0

    y = [[Decimal(0)] * ycount for i in range(n + 1)]
    for i in range(ycount):
        y[0] = y0

    k1 = [Decimal(0)] * ycount
    k2 = [Decimal(0)] * ycount
    k3 = [Decimal(0)] * ycount
    k4 = [Decimal(0)] * ycount

    for i in range(n):
        k1 = scalarVectorMultiplication(h, dydx(x[i],                     y[i]))
        k2 = scalarVectorMultiplication(h, dydx(x[i] + Decimal(0.5) * h,  vectorSum(y[i], scalarVectorMultiplication(Decimal(0.5), k1))))
        k3 = scalarVectorMultiplication(h, dydx(x[i] + Decimal(0.5) * h,  vectorSum(y[i], scalarVectorMultiplication(Decimal(0.5), k2))))
        k4 = scalarVectorMultiplication(h, dydx(x[i] + h,                 vectorSum(y[i], k3)))

        y[i + 1] = vectorSum(
            y[i],
            scalarVectorMultiplication(
                Decimal(1.0 / 6.0),
                vectorSum(
                    k1,
                    vectorSum(
                        scalarVectorMultiplication(Decimal(2),k2), 
                        vectorSum(
                            scalarVectorMultiplication(Decimal(2), k3),
                            k4)
                        )
                    )
                )
            )
        x[i + 1] = x[i] + h
    return x, y
