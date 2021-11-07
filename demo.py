from matrix import Matrix, solve_with_Kramer, solve_with_matrix_transformation, solve_with_Haus_method

m = Matrix([
    [1, 9, 5, 4],
    [4, 2, 3, 7],
    [7, 2, 4, 9],
    [8, 1, 4, 2]
])
k = Matrix([
    [1, 2, 3],
    [4, 5, 6],
    [2, 1, 0]
])
n = Matrix([
    [1, 2],
    [4, 0],
    [7, 0],
    [8, 0]
])

print(m)
print(n)
print(k)

print('m.shape')  # (4, 4)
print(m.shape)  # (4, 4)

print('-' * 100)
print('m[1]')  # [7, 2, 4, 9]
print(m[1])  # [7, 2, 4, 9]

print('-' * 100)
print('Срез от 1 строки до 3, от 2 столбца до 4 строки (не включительно)')
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

print('-' * 100)
print('Приведение к треугольному виду')
print(m.to_triangle())

print('-' * 100)
print('Ранг матрицы')
print(k.rank)

print('-' * 100)
print('##Решение Линейных Алгебраических уровнений\n')


A = Matrix([[2, -1, 0], [1, 2, -1], [0, 1, 1]])
print('Матрица коэфицентов')
print(A)

b = Matrix([[1], [2], [2]])
print('Матрица правой части')
print(b)

print('###Методом Гауса')
x = solve_with_Haus_method(A, b)
print(x)

print('###Методом Крамера')
x = solve_with_Kramer(A, b)
print(x)

print('###Методом матричных преобразований')
x = solve_with_matrix_transformation(A, b)
print(x)
