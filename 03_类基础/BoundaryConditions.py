def Dirichlet(left, right, A, b):
    n = len(b)
    A[0] = 0
    A[0][0] = 1
    b[0] = left
    A[n-1] = 0
    A[n-1][n-1] = 1
    b[n-1] = right
    return A, b


def Left_Neumann(left, right, A, b, c, val):
    n = len(b)
    b[0] -= left*c(val)
    A[n - 1] = 0
    A[n - 1][n - 1] = 1
    b[n - 1] = right
    return A, b


def Right_Neumann(left, right, A, b, c, val):
    n = len(b)
    A[0] = 0
    A[0][0] = 1
    b[0] = left
    b[n - 1] += right*c(val)
    return A, b


def Neumann(left, right, A, b, c, val1, val2):
    pass


def Right_Robbin(left, pb, qb, A, b, c, val):
    n = len(b)
    A[0] = 0
    A[0][0] = 1
    b[0] = left
    b[n - 1] += pb*c(val)
    A[n - 1][n - 1] += qb*c(val)
    return A, b