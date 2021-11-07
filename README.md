# Matrix

## m = Matrix([
        [1, 9, 5, 4],
        [4, 2, 3, 7],
        [7, 2, 4, 9],
        [8, 1, 4, 2]
    ])
## n = Matrix([
        [1, 2],
        [4, 0],
        [7, 0],
        [8, 0]
    ])

## m.shape
(4, 4)
----------------------------------------------------------------------------------------------------
## m[1]
[1, 9, 5, 4]
----------------------------------------------------------------------------------------------------
## Срез от 1 строки до конца, от 2 столбца до 4 строки (не включительно)
[9, 5]
[2, 3]
----------------------------------------------------------------------------------------------------
## m.determinant
-142
----------------------------------------------------------------------------------------------------
## m + 2
[3, 11, 7, 6]
[6, 4, 5, 9]
[9, 4, 6, 11]
[10, 3, 6, 4]
----------------------------------------------------------------------------------------------------
## m + m
[2, 18, 10, 8]
[8, 4, 6, 14]
[14, 4, 8, 18]
[16, 2, 8, 4]
----------------------------------------------------------------------------------------------------
## m - 2
[-1, 7, 3, 2]
[2, 0, 1, 5]
[5, 0, 2, 7]
[6, -1, 2, 0]
----------------------------------------------------------------------------------------------------
## m - m
[0, 0, 0, 0]
[0, 0, 0, 0]
[0, 0, 0, 0]
[0, 0, 0, 0]
----------------------------------------------------------------------------------------------------
## Обратная матрица округленная до сотых
[0.09, -1.49, 1.15, -0.15]
[0.24, -1.73, 1.31, -0.31]
[-0.23, 3.39, -2.68, 0.68]
[-0.02, 0.04, 0.12, -0.12]
----------------------------------------------------------------------------------------------------
## m^-1 * m
[1, 0, 0, 0]
[0, 1, 0, 0]
[0, 0, 1, 0]
[0, 0, 0, 1]
----------------------------------------------------------------------------------------------------
## m ^ 3
[1704, 1274, 1395, 1555]
[1331, 1110, 1152, 1463]
[1870, 1485, 1583, 2077]
[1493, 885, 1116, 1573]
----------------------------------------------------------------------------------------------------
## m * 2
[2, 18, 10, 8]
[8, 4, 6, 14]
[14, 4, 8, 18]
[16, 2, 8, 4]
----------------------------------------------------------------------------------------------------
## m * n
[104, 2]
[89, 8]
[115, 14]
[56, 16]
----------------------------------------------------------------------------------------------------
## Транспонирование
[1, 4, 7, 8]
[9, 2, 2, 1]
[5, 3, 4, 4]
[4, 7, 9, 2]
----------------------------------------------------------------------------------------------------
## Минор 2 строка, 3 столбец
481
----------------------------------------------------------------------------------------------------
## Алгебраическое дополнение 2 строка, 3 столбец
-481
