def degree4Vandermond(vectorX):
  """ Construct the degree 4 vandermond matrix using vectorX. 

  We do this by first iterating from 0 to 4 since we are constructing the degree 4 vandermond matrix. At the beginning of each iteration we take a temporary empty set and start to iterate over the elements in the inputed vector. we take each element in this vector and raise it to the current 0 to 4 iteration that we are in and append it to the temporary vector. Whe we have done this for every element, we then append the temporary vector to the output matrix as a list of colums. When complete, we return the matrix as the output.

  Args:
    vectorX: a list of numbers representing a vector. 
  Returns:
    matrixA : degree 4 vandermond matrix generated from the inputed vectorX. The matrix is a list of column vectors. 
  """

  matrixA = []
  for exponent in range(5):
    temp = []
    for element in range(len(vectorX)):
      temp.append(vectorX[element]**exponent)
    matrixA.append(temp)
  return matrixA

def conjugateTranspose(matrixA):
  """ Compute the conjugate transpose of the inputed matrix. 

  We first build the 0 n x m matrixB according to the transpose measurements of the inputed m x n matrixA. We then iterate over the colums of the matrix as index1. for each column that we are in, we then iterate over the elements in the columns as index2. We store the element in matrixA[index2][index1] into the element in matrixB[index1][index2]. We then take the column vectors in matrixB and replace it with the conjugate of the elements by calling the conjugate function.

  Args:
    matrixA: A matrix of size m x n with column vectors.
  Returns:
    matrixB: A matrix of size n x m with column vectors. 

  """

  matrixB = [0] * len(matrixA[0])
  for rows in range(len(matrixA[0])):
    matrixB[rows] = [0] * len(matrixA)
  for index1 in range(len(matrixB)):
    for index2 in range(len(matrixB[0])):
      matrixB[index1][index2] = matrixA[index2][index1]   
  print(matrixB) 
  for element in range(len(matrixB)):
    matrixB[element] = conjugate(matrixB[element])
  return matrixB

def conjugate(vector):
  """ compute the conjugate of all the elements in the vector. 

  Take each element in the vector and copy the real part of that number and save it to a temporary placement holder. We then take that temporary element and subtract the imaginary part of the element in the vector. We then take that number append it to the output vector. 

  Args:
    vector: A list of numbers
  Returns:
    conjugate: A vector as a list of numbers that is the conjugate of the inputed vector.
  """
  conjugate = []
  temp=[]
  for element in range(len(vector)):
    temp = vector[element].real
    temp = temp - ((vector[element].imag)*1j)
    conjugate.append(temp)
  return conjugate

def vectorAddition(vector1, vector2):
  """ Add two vectors together.

  Take the associated elements of each vector and add them 
  together. Take that sum and save it to the associated element location of the new vector. 

  Args:
    vector1: vector with numbers as the elements
    vector2: vector with numbers as the elements
  Returns:
    newVector: vector with numbers as the elements.
  """
  newVector = vector1

  for element in range(len(vector1)):
    newVector[element] = vector1[element] + vector2[element]
  return newVector

def scalarVectorMultiply(vector, scalar):
  """ multiply a vetor by a scalar

  Take each element of the vector and multiply it by the 
  scalar then save it to the new vector

  Args:
    vector: A vector with elements in the form of numbers.
    scalar: A number.
  Returns:
    newVector: a vector with elements as numbers.

  """

  newVector = vector
  for element in range(len(vector)):
    newVector[element] = scalar * vector[element]
  return newVector

def matrixVectorMul(matrix,vector1):
  """Multiply the inputed matrix by the inputed vector. 

  We iterate over the columns of the matrix and multiply that column with the associated element in the vector. We call other functions to complete the process. 

  Args:
    matrix: a m x n matrix with column vectors
    vector1: a vector with numbers as the elements
  Returns
    vector2: a vector of length m
  """
  vector2=[0]*len(matrix[0])
  for element in range(len(vector1)):
    vector2 = vectorAddition(vector2,scalarVectorMultiply(matrix[element],vector1[element]))
  return vector2

def twoNorm(vector):
  """Calculates the 2 norm of a vector.

  Sums the squares of the corresponding elements of a given vector and returns the square root of the sum.

  Args:
    vector: a list of numbers
  Returns:
    A scalar which is the 2-norm of the given vector
  """

  result = 0
  for element in range(len(vector)):
    result = result + (vector[element].real)**2 + (vector[element].imag)**2
  result = result**(1/2)
  return result

def normalize(vector):
  """ Normalize a vector with respect to its 2 normalize

  Takes the inputed vector and its is going to normalize it by first checking the 2 norm of the vector. If the 2 norm is 0, then we cant compute the normalization and print out an error message. If the 2 norm is 1, then it is just going to return the original vector since it is already normalized. If the 2-norm is none of the above cases, then it normalizes the vector by calling the scalarVectorMultiply function  with 1/norm and the original vector as the inputs. 

  Args:
    vector: with elements as a numbers
  
  Returns:
    normalized vector represented as a list of numbers.

  """

  norm = twoNorm(vector)
  if (norm == 0 ):
    print("invalid input")
  elif (norm == 1):
    return vector
  else:
    return scalarVectorMultiply(vector, 1/norm)

def backsub(matrixR,vector_b):
  """ use back substitution to get back the vector x that solves Ax = b. 

  We use the information from the following rows of the current vector that we are in to solve the element X[iterator] corresponding to the current vector. 

  Args: 
    matrixR : a upper triangular matrix with column vectors. 
    Vector_b : a vector with numbers as elements

  Returns:
    result: a vector x that solves Ax=b with numbers as elements. 

    """


  result = vector_b
  for iterator in range(len(matrixR[0])):
    alpha = (len(matrixR[0]-1))
    result[alpha-iterator] = vector_b[alpha-iterator]- summation(matrixR[alpha-iterator][iterator],result[iterator])*(1/matrixR[alpha-iterator][alpha-iterator]

  return result

  
def modifiedGS(matrixA):
  """ Compute matrixQ witch is a unitary matrix and matrixR witch is a upper triangular matrix.

  Args:
    matrixA : a m x n matrix as a list of colun vectors. 
  Returns:
    matrixQ : a unitary matrix with column vectors.
    matrixR : a upper triangular matrix. 

  """
  matrixV = matrixA
  matrixQ =[0] * len(matrixA)
  matrixR = [0]*len(matrixA)
  for elements in range(len(matrixA[0])):
    matrixR[elements] = [0] * len(matrixA)
  for j in range(len(matrixA)):
    matrixR[j][j]= twoNorm(matrixV[j])
    matrixQ[j] = scalarVectorMultiply(matrixV[j],1/matrixR[j][j])
    for k in range(j+1,len(matrixA)):
      matrixR[k][j] = dotProduct(conjugate(matrixQ[j]),matrixV[k])
      matrixV[k] = matrixV[k] + scalarVectorMultiply(matrixQ[j],-matrixR[k][j])
  return matrixR , matrixQ

def dotProduct(vector1,vector2):
  """ Compute the dot product of 2 given vectors of compatible dimmensions.

  Take the corresponding element of each vector and multiply them together, and then compute the sum of the products. 

  Args:
    vector1:a list with numbers as the elements
    vector2:a list with numbers as the elements  
  
  Returns:
    scalar:the computed dot product of the 2 vectors

  """

  scalar = 0
  for element in range(len(vector1)):
    scalar = scalar + vector1[element] * vector2[element]
  return scalar

def d4interpolant(beta,x):
  """ Print the degree 4 polynomial.

  It takes the elements in the iputed vector and uses them as the coefficients of the polynomial. The x is the variable that it is going to use. 

  Args:
    beta: vector with numbers as its elements
    x : is a scalar. 

  returns:
    the polynomial created with the vector and the variable. 

  """


  return (beta[0]+beta[1]*x +beta[2]* (x**2) + beta[3]*(x**3) + beta[4]*(x**4))
