def leastSquares(vectorX, vectorY):
  """Compute the least squares of the 2 vectors.

  Import the basicFunctions file for the subalgorithms that are being used in this process. We first need to construct the degree 4 vandermond matrix with the inputed VectorX. After we have that matrix, we are going to construct MatrixQ and MatrixR using the modified Gram Shmidt process on matrixA. We then use backsubstitution to calculate our beta vector and input that intot the degree4Interpolant function.

  Args: 
    vectorX : A vector with numbers as its elements
    vectorY : A vector with numbers as its elements
  
  Returns:
    result : a scalar
  """
  import basicFunctions.py as bF

  matrixA = bF.degree4Vandermond(vectorX)

  [matrixQ,matrixR] = bF.modifiedGS(matrixA)

  beta = bF.backSub(matrixR, bF.matrixVectorMul(bF.conjugateTranspose(matrixQ),vectorY))
  
  result = bF.degree4Interpolant(beta,d)

  return result


