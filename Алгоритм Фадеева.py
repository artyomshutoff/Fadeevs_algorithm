#artyomshutoff

import time
import numpy as np

A = np.array([[2, -1, 2],
			  [1, 2, -1]], dtype=int)

def tr(matrix):
	"""
	След матрицы
	"""
	out = 0
	for i in range(len(matrix)):
		if len(matrix[i]) == i:
			return out
		out += matrix[i][i]
	return out

def T(matrix):
	"""
	Транспанирование матриц
	"""
	out_matrix = [[] for i in range(len(matrix[0]))]
	for i in range(len(matrix[0])):
		for j in range(len(matrix)):
			out_matrix[i].append(matrix[j][i])
	return np.array(out_matrix, dtype=float)

def matrix_multiplication(A, B):
	"""
	Умножение матриц
	"""
	if len(A[0]) != len(B):
		raise "Кол-во столбцов A не совпадате с кол-вом строк B"
	matrix = [[0 for i in range(len(B[0]))][:] for i in range(len(A))]
	for i in range(len(matrix)):
		for j in range(len(matrix[0])):
			for row in range(len(A[i])):
				matrix[i][j] += A[i][row] * B[row][j]
	return np.array(matrix, dtype=float)

def E_matrix(matrix):
	"""
	Единичная матрица из обычной
	"""
	out_matrix = []
	for i in range(len(matrix)):
		if i < len(matrix[0]):
			out_matrix.append([0 for j in range(len(matrix[0]))])
			out_matrix[i][i] = 1		
	return np.array(out_matrix, dtype=float)

def alpha_multiplication(a, matrix):
	"""
	Умножение на коэфф. 𝞪
	"""
	for i in range(len(matrix)):
		for j in range(len(matrix[i])):
			matrix[i][j] *= a
	return np.array(matrix, dtype=float)

def plus_minus_matrix(A, B, operation = '+'):
	"""
	Вычитание и прибавление матриц
	"""
	matrix = A[:]
	if len(A[0]) != len(B[0]) or len(A) != len(B):
		raise "Матрицы не равны!"

	for i in range(len(A)):
		for j in range(len(A[i])):
			matrix[i][j] = A[i][j] - B[i][j] if operation == '-' else A[i][j] + B[i][j]
	
	return np.array(matrix, dtype=float)

def Fadeevs_algorithm(A, F=1, F_1=1, fi=1, E=1, k = 1):

	def check_last(fi, E, A, F, k):
		A_T_A_F = matrix_multiplication(matrix_multiplication(T(A), A), F)
		E = E_matrix(A_T_A_F)
		F = plus_minus_matrix(alpha_multiplication(fi, E), A_T_A_F, '-')
		fi = tr(matrix_multiplication(matrix_multiplication(T(A), A), F)) / k
		return fi

	if k == 1:
		F = E_matrix(matrix_multiplication(T(A), A))
		fi = tr(matrix_multiplication(T(A), A))
		F_1 = F[:]
		E = F_1[:]
	else:
		A_T_A_F = matrix_multiplication(matrix_multiplication(T(A), A), F)
		E = E_matrix(A_T_A_F)
		F = plus_minus_matrix(alpha_multiplication(fi, E), A_T_A_F, '-')
		fi = tr(matrix_multiplication(matrix_multiplication(T(A), A), F)) / k
		if check_last(fi, E, A, F, k+1) == 0:
			return matrix_multiplication(alpha_multiplication(1/fi, F), T(A))

	return Fadeevs_algorithm(A, F, F_1, fi, E, k + 1)

print('Псевдообратная матрица:')
print('')
print(Fadeevs_algorithm(A))