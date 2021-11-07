from fractions import Fraction

from matrix import Matrix


def solve_with_Haus_method(A: Matrix, b: Matrix) -> Matrix:
    assert A.shape[0] == b.shape[0], 'The number of rows of matrix A must be equal to number of rows b'
    assert b.shape[1] == 1, 'The number of columns b must be 1'

    _x = []
    A_tr = A.copy().append(b, axis=1).to_triangle()
    assert A_tr.rank == A.rank, 'Rank matrix don`t match, no solutions'
    # A_tr = A_tr[:, :A_tr.shape[1]]
    # print(A_tr)
    for i in range(A_tr.shape[0], 0, -1):
        y = list(A_tr[i, A_tr.shape[1]])[0][0]
        for n, x in enumerate(_x):
            y -= x * list(A_tr[i, A_tr.shape[1] - 1 - n])[0][0]
        xi = list(A_tr[i, A_tr.shape[1] - 1 - len(_x)])[0][0]
        _x.append(Fraction(y, xi))

    return Matrix([[x] for x in _x[::-1]])


def solve_with_Kramer(A: Matrix, b: Matrix):
    assert A.shape[0] == A.shape[1], 'Matrix must be square'
    assert A.shape[0] == b.shape[0], 'The number of rows of matrix A must be equal to number of rows b'
    assert b.shape[1] == 1, 'The number of columns b must be 1'

    determinants_A = []
    for j in range(1, A.shape[1] + 1):
        _A = A.copy()
        _A = _A[:, :j].append(b, axis=1).append(_A[:, j + 1:], axis=1)
        # print(_A)
        determinants_A.append(_A.determinant)
    return Matrix([[det / A.determinant] for det in determinants_A])


def solve_with_matrix_transformation(A, b):
    return (A ** -1) * b


if __name__ == '__main__':
    A = Matrix([[2, -1, 0], [1, 2, -1], [0, 1, 1]])
    print(A)

    b = Matrix([[1], [2], [2]])
    print(b)

    x = solve_with_Haus_method(A, b)
    print(x)

    x = solve_with_Kramer(A, b)
    print(x)

    x = solve_with_matrix_transformation(A, b)
    print(x)
