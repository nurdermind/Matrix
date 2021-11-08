from collections import OrderedDict
from fractions import Fraction

from matrix import Matrix


def solve_with_Haus_method(A: Matrix, b: Matrix) -> tuple[Matrix, str]:
    assert A.shape[0] == b.shape[0], 'The number of rows of matrix A must be equal to number of rows b'
    assert b.shape[1] == 1, 'The number of columns b must be 1'

    _x = OrderedDict()
    A_tr = A.copy().append(b, axis=1).to_triangle()
    assert A_tr.rank == A.rank, 'Rank matrix don`t match, no solutions'
    # A_tr = A_tr[:, :A_tr.shape[1]]
    remark = None
    mode_t = True if any(list(A_tr[A_tr.shape[0], :A_tr.shape[1] - 1])[0]) else False

    # print(mode_t)
    # print(A_tr)
    for i in range(A_tr.shape[0], 0, -1):
        if not any(list(A_tr[i])[0]):
            continue
        if mode_t:
            # for x_j in (list(A_tr[i, :A_tr.shape[1]])[0][::-1]):
            for j in range(1, A_tr.shape[1]):
                x_j = list(A_tr[i, j])[0][0]
                if x_j != 0:
                    if not remark:
                        remark = f"* x{j}"
                    _x[j] = 1
            mode_t = False

        if i not in _x:
            y = list(A_tr[i, A_tr.shape[1]])[0][0]
            for j, x_j in _x.items():
                y -= x_j * list(A_tr[i, j])[0][0]
            xi = list(A_tr[i, A_tr.shape[1] - 1 - len(_x)])[0][0]
            _x[i] = (Fraction(y, xi))
    return Matrix([[_x[key]] for key in sorted(_x)]), remark


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
    A = Matrix([
        [2, 5, -8, 1],
        [4, 3, -9, 1],
        [2, 3, -5, 2],
        [1, 8, -7, 1],
        [2, 3, 5, 2],
    ])
    print(A)
    # print(A.to_triangle())

    b = Matrix([[8],
                [9],
                [7],
                [12],
                [17],
                ])
    print(b)

    x = solve_with_Haus_method(A, b)
    print(*x)

    # x = solve_with_Kramer(A, b)
    # print(x)
    #
    # x = solve_with_matrix_transformation(A, b)
    # print(x)
