def dirichlet1d(left, right, A, b):
    n = len(b)
    A[0] = 0
    A[0][0] = 1
    b[0] = left
    A[n-1] = 0
    A[n-1][n-1] = 1
    b[n-1] = right
    return A, b