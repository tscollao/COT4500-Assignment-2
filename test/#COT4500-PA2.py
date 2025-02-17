#Neville's Method
def neville(x_values, f_values, x):
    n = len(x_values)
    P = [[0 for _ in range(n)] for _ in range(n)]
    
    for i in range(n):
        P[i][i] = f_values[i]

    for i in range(1, n):
        for j in range(n-i):
            P[j][j+i] = ((x - x_values[j+i]) * P[j][j+i-1] - (x - x_values[j]) * P[j+1][j+i]) / (x_values[j+i] - x_values[j])
    
    return P[0][n-1]

# Given data
x_values = [3.6, 3.8, 3.9]
f_values = [1.675, 1.436, 1.318]
x = 3.7

# Calculate f(3.7)
result = neville(x_values, f_values, x)
print(f"Using Neville's Method, f(3.7) = {result}")

#Newton's foward method
def forward_difference_table(x_values, f_values):
    n = len(x_values)
    table = [f_values.copy()]  # First row is just f(x) values
    
    # Calculate forward differences
    for i in range(1, n):
        diff = []
        for j in range(n - i):
            diff.append(table[i - 1][j + 1] - table[i - 1][j])
        table.append(diff)
    
    return table

def newton_forward_polynomial(x_values, f_values, degree, x):
    # Create the forward difference table
    table = forward_difference_table(x_values, f_values)
    
    # Start with the first term, f(x_0)
    result = f_values[0]
    product_term = 1
    
    # Compute the terms for polynomial for the given degree
    for i in range(1, degree + 1):
        product_term *= (x - x_values[i - 1])  # Multiply with (x - x_i-1)
        result += table[i][0] * product_term  # Add the corresponding forward difference term

    return result

# Given data
x_values = [7.2, 7.4, 7.5, 7.6]
f_values = [23.5492, 25.3913, 26.8224, 27.4589]

# Print polynomial approximations for degrees 1, 2, and 3 at x = 7.3
for degree in range(1, 4):
    result = newton_forward_polynomial(x_values, f_values, degree, 7.3)
    print(f"Polynomial approximation of degree {degree} at f(7.3) = {result}")

#Divided difference method
def hermite_divided_difference(x_values, f_values, f_prime_values):
    n = len(x_values)
    table = [[0 for _ in range(2*n)] for _ in range(2*n)]
    
    # Fill in initial values for f(x) and f'(x)
    for i in range(n):
        table[2*i][0] = f_values[i]
        table[2*i + 1][0] = f_values[i]
        table[2*i][1] = f_prime_values[i]
    
    # Construct the divided difference table
    for j in range(2, 2*n):
        for i in range(2*n - j):
            table[i][j] = (table[i+1][j-1] - table[i][j-1]) / (x_values[i//2 + j//2] - x_values[i//2])
    
    return table

# Given data
x_values = [3.6, 3.8, 3.9]
f_values = [1.675, 1.436, 1.318]
f_prime_values = [-1.195, -1.188, -1.182]

# Calculate the Hermite divided difference table
table = hermite_divided_difference(x_values, f_values, f_prime_values)

# Print the Hermite polynomial approximation table
for row in table:
    print(row)

#cubic spline interpolation
import numpy as np

def cubic_spline_interpolation(x_values, f_values):
    n = len(x_values)
    A = np.zeros((n, n))
    b = np.zeros(n)
    
    # Fill the matrix A and vector b based on the system of equations for cubic splines
    h = np.diff(x_values)
    
    # Matrix A
    A[0, 0] = 1
    A[n-1, n-1] = 1
    for i in range(1, n-1):
        A[i, i-1] = h[i-1]
        A[i, i] = 2 * (h[i-1] + h[i])
        A[i, i+1] = h[i]
    
    # Vector b
    for i in range(1, n-1):
        b[i] = 3 * ((f_values[i+1] - f_values[i]) / h[i] - (f_values[i] - f_values[i-1]) / h[i-1])
    
    # Solve the system of equations
    x = np.linalg.solve(A, b)
    
    return A, b, x

# Given data
x_values = [2, 5, 8, 10]
f_values = [3, 5, 7, 9]

# Solve for matrix A, vector b, and vector x
A, b, x = cubic_spline_interpolation(x_values, f_values)

# Print the results
print("Matrix A:")
print(A)
print("\nVector b:")
print(b)
print("\nVector x:")
print(x)
