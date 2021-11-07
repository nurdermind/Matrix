from copy import deepcopy
from fractions import Fraction
from typing import Union, Tuple

from tabulate import tabulate


class Matrix(list):
    """
    The matrix
    """

    def __init__(self, body: Union[int, list], limit_denominator=20):
        super(Matrix, self).__init__()
        if type(body) is int:
            self._body = [[body]]
        elif body and type(body[0]) is int:
            self._body = [deepcopy(body)]
        else:
            self._body = deepcopy(body)
        self._determinant = None
        self._limit_denominator = limit_denominator

    @property
    def shape(self) -> Tuple[int, int]:
        i = len(self._body)
        j = len(self._body[0]) if self._body else 0
        return i, j

    @property
    def determinant(self):
        if self._determinant:
            return self._determinant
        assert self.shape[0] - self.shape[1] == 0, 'Matrix must be square'
        if self.shape == (1, 1):
            return self._body[0][0]
        _determinant = 0
        for j, a in enumerate(self[1]._body[0]):
            j += 1
            a_dop = self.get_algebraic_completion((1, j))
            _determinant += a * a_dop
        return _determinant

    def get_minor(self, position: Tuple[int, int]) -> float:
        _new_matrix = Matrix([[0 for _ in range(self.shape[1] - 1)] for _ in range(self.shape[0] - 1)])
        i_exclude = position[0]
        j_exclude = position[1]

        for i in range(1, self.shape[0] + 1):
            if i == i_exclude:
                continue
            if i > i_exclude:
                _new_matrix[i - 1] = self[i, :j_exclude]._body[0] + self[i, j_exclude + 1:]._body[0]
            else:
                _new_matrix[i] = self[i, :j_exclude]._body[0] + self[i, j_exclude + 1:]._body[0]

        return _new_matrix.determinant

    def get_algebraic_completion(self, position: Tuple[int, int]):
        i, j = position
        return (-1) ** (i + j) * self.get_minor(position)

    @property
    def T(self):
        _new_matrix = Matrix(deepcopy(self._body))
        for i in range(_new_matrix.shape[0]):
            for j in range(_new_matrix.shape[1]):
                if i == j:
                    break
                _new_matrix._body[i][j], _new_matrix._body[j][i] = \
                    _new_matrix._body[j][i], _new_matrix._body[i][j]
        return _new_matrix

    def apply(self, foo):
        _new_matrix = Matrix(deepcopy(self._body))
        for i in range(_new_matrix.shape[0]):
            for j in range(_new_matrix.shape[1]):
                _new_matrix._body[i][j] = foo(self._body[i][j])
        return _new_matrix

    def append(self, other: 'Matrix', axis=0) -> 'Matrix':
        assert axis in (0, 1), 'axis must be 0 or 1'
        if axis == 0:
            assert self.shape[0] == other.shape[0], 'The number of lines must match'
            _other_matrix = Matrix(deepcopy(other._body))
            for i in range(self.shape[0]):
                self._body[i] += other._body[i]
        elif axis == 1:
            assert self.shape[1] == other.shape[1], 'The number of column must match'
            self._body += deepcopy(other._body)
        return self

    def __iter__(self):
        self.__iter = iter(self._body)
        return self

    def __next__(self):
        return next(self.__iter)

    def __add__(self, other: Union['Matrix', int, float]):
        _new_matrix = Matrix(deepcopy(self._body))
        if type(other) in (int, float):
            for i in range(self.shape[0]):
                for j in range(self.shape[1]):
                    _new_matrix._body[i][j] = self._body[i][j] + other
        else:
            assert self.shape == other.shape
            for i in range(self.shape[0]):
                for j in range(self.shape[1]):
                    _new_matrix._body[i][j] = self._body[i][j] + other._body[i][j]
        return _new_matrix

    def __sub__(self, other):
        return self + (other * -1)

    def __mul__(self, other: Union['Matrix', int]):
        if type(other) in (int, float, Fraction):
            _new_matrix = Matrix(deepcopy(self._body))
            for i in range(self.shape[0]):
                for j in range(self.shape[1]):
                    _new_matrix._body[i][j] = self._body[i][j] * Fraction(other)
            return _new_matrix
        elif self.shape[1] == other.shape[0]:
            _new_matrix = Matrix([[0 for _ in range(other.shape[1])] for _ in range(self.shape[0])])
            for i in range(_new_matrix.shape[0]):
                for j in range(_new_matrix.shape[1]):
                    for k in range(self.shape[1]):
                        _new_matrix._body[i][j] += self._body[i][k] * Fraction(other._body[k][j])

            return _new_matrix

    def __pow__(self, power, modulo=None):
        _new_matrix = Matrix([[1 if i == j else 0 for j in range(self.shape[1])] for i in range(self.shape[0])])
        if power >= 0:
            for _ in range(power):
                _new_matrix = _new_matrix * self
            return _new_matrix
        if power == -1:
            assert self.determinant != 0
            _new_matrix = Matrix(deepcopy(self._body))
            for i in range(self.shape[0]):
                for j in range(self.shape[1]):
                    _new_matrix._body[i][j] = self.get_algebraic_completion((i + 1, j + 1))
            return _new_matrix.T * Fraction(1, self.determinant)

    def __len__(self):
        return len(self._body)

    def __getitem__(self, val: Union[int, slice, tuple]):
        # print(val)
        _res = None
        if type(val) is tuple and len(val) == 2:
            x, y = val
            if type(x) is slice:
                pass
            elif type(x) is int:
                x = slice(x, x + 1)
            else:
                raise Exception("Slice is wrong!")
            if type(y) is slice:
                pass
            elif type(y) is int:
                y = slice(y, y + 1)
            else:
                raise Exception("Slice is wrong!")
        elif type(val) is slice:
            x, y = val, None
        elif type(val) is int:
            x, y = slice(val, val + 1), None
        else:
            raise Exception("Slice is wrong!")
        # print(x, y)
        if type(x) is slice:
            i, j, k = x.start, x.stop, x.step
            i = i - 1 if i is not None else None
            j = j - 1 if j is not None else None
            _res = self._body[i:j:k]

        if type(y) is int:
            y = slice(y, y + 1)
        if type(y) is slice:
            i, j, k = y.start, y.stop, y.step
            i = i - 1 if i is not None else None
            j = j - 1 if j is not None else None
            for _index, _value in enumerate(_res):
                _res[_index] = _value[i:j:k]
        result = Matrix(_res)
        # return result._body[0][0] if result.shape == (1, 1) else \
        #         result._body[0] if result.shape[0] == 1 else result
        return result

    def __setitem__(self, key, value):
        self._body[key - 1] = value

    def __str__(self):
        # s = ''
        # for row in self._body:
        #     s += '[ '
        #     for item in row:
        #         s += f"{item} "
        #     s += ']\n'
        return tabulate([[str(item) for item in row] for row in self._body], tablefmt='fancy_grid', )
        # return str(self._body)
        return s

    def __repr__(self):
        return str(self._body)


if __name__ == '__main__':
    m = Matrix([
        [1, 9, 5, 4],
        [4, 2, 3, 7],
        [7, 2, 4, 9],
        [8, 1, 4, 2]
    ])
    n = Matrix([
        [1, 2],
        [4, 0],
        [7, 0],
        [8, 0]
    ])

    print('m.shape')  # (4, 4)
    print(m.shape)  # (4, 4)

    print('-' * 100)
    print('m[1]')  # [7, 2, 4, 9]
    print(m[1])  # [7, 2, 4, 9]

    print('-' * 100)
    print('Срез от 1 строки до конца, от 2 столбца до 4 строки (не включительно)')
    print(m[:3, 2:4])  # [[9, 5] , [2, 3]]

    print('-' * 100)
    print('m.determinant')  # -142
    print(m.determinant)  # -142

    print('-' * 100)
    print('m + 2')
    print(m + 2)

    print('-' * 100)
    print('m + m')
    print(m + m)

    print('-' * 100)
    print('m - 2')
    print(m - 2)

    print('-' * 100)
    print('m - m')
    print(m - m)

    print('-' * 100)
    print('Обратная матрица')
    print((m ** -1))

    print('-' * 100)
    print('m^-1 * m')
    print((m * (m ** -1)))

    print('-' * 100)
    print('m ^ 3')
    print((m ** 3))

    print('-' * 100)
    print('m * 2')
    print(m * 2)

    print('-' * 100)
    print('m * n')
    print(m * n)

    print('-' * 100)
    print('Транспонирование')
    print(m.T)

    print('-' * 100)
    print('Минор 2 строка, 3 столбец')
    print(m.get_minor((2, 3)))

    print('-' * 100)
    print('Алгебраическое дополнение 2 строка, 3 столбец')
    print(m.get_algebraic_completion((2, 3)))
