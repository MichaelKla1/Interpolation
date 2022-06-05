import math

points = [[0,0],[1,0.8415],[2,0.9093],[3,0.1411],[4,-0.7568],[5,-0.9589],[6,-0.2794]] #left number in pair - x, right number - y
XToFind = 2.5

def linearInterpolation(points,XToFind):
    """
    Linear interpolation. Returns f(XToFind) based on the points table
    :param points: table of points
    :param XToFind: X to find f(X) based on the points table
    :return: f(XToFind) based on the points table
    """
    i = 0
    found = False
    while(i < len(points) - 1):
        if (points[i][0] - XToFind)*(points[i+1][0] - XToFind) < 0: #if this is true, X coordinate to find is between 2 these points
            found = True
            result = ((points[i][1] - points[i+1][1])*XToFind)/(points[i][0] - points[i+1][0])+(points[i+1][1]*points[i][0]-points[i][1]*points[i+1][0])/(points[i][0] - points[i+1][0])
            break
        i = i + 1
    if not found:
        return None
    else:
        return result

def polynomialInterpolation(points,XToFind):
    """
    Polynomial interpolation. Returns f(XToFind) based on the points table
    :param points: table of points
    :param XToFind: X to find f(X) based on the points table
    :return: f(XToFind) based on the points table
    """
    i = 0
    matrix = []
    resvec = []
    while (i < len(points)):
        lineToAppend = []
        j=0
        while (j<len(points)):
            lineToAppend.append(math.pow(points[i][0],j))
            j = j + 1
        matrix.append(lineToAppend)
        resvec.append(points[i][1])
        i = i + 1
    prefixValues = solutionOfMatrixByGaussElimination(matrix,resvec)
    i = 0
    result = 0
    while (i < len(prefixValues)):
        result = result + math.pow(XToFind,i)*prefixValues[i]
        i = i + 1
    return result

def lagrangeInterpolation(points,XToFind):
    """
    Interpolation by lagrange method. Returns f(XToFind) based on the points table
    :param points: table of points
    :param XToFind: X to find f(X) based on the points table
    :return: f(XToFind) based on the points table
    """
    i = 0
    LArray = []
    while (i < len(points)):
        j = 0
        resultToAppendToLArray = 1
        while(j < len(points)):
            if not i == j:
                resultToAppendToLArray = resultToAppendToLArray * ((XToFind - points[j][0])/(points[i][0] - points[j][0]))
            j = j + 1
        LArray.append(resultToAppendToLArray)
        i = i + 1
    result = 0
    i = 0
    while i < len(LArray):
        result = result + points[i][1] * LArray[i]
        i = i + 1
    return result

def nivulInterpolation(points,XToFind):
    """
    Interpolation by nivul method. Returns f(XToFind) based on the points table
    :param points: table of points
    :param XToFind: X to find f(X) based on the points table
    :return: f(XToFind) based on the points table
    """
    PArray = []
    tempPArray = []
    step = 0
    while step < len(points):
        i = 0
        while (step + i < len(points)):
            if step == 0:
                tempPArray.append(points[i][1])
            else:
                valueToAppend = ((XToFind - points[i][0])*PArray[i+1] - (XToFind - points[i+step][0]) * PArray[i]) / (points[i+step][0] - points[i][0])
                tempPArray.append(valueToAppend)
            i = i + 1
        PArray = tempPArray.copy()
        tempPArray.clear()
        step = step + 1
    return PArray[0]

def splineInterpolation(points,XToFind):
    """
    Interpolation by spline method. Returns f(XToFind) based on the points table
    :param points: table of points
    :param XToFind: X to find f(X) based on the points table
    :return: f(XToFind) based on the points table
    """
    lambdaVec = []
    muVec = []
    dVec = []
    hVec = []
    #build hVec
    i = 0
    while i < len(points) - 1:
        hVec.append(points[i+1][0] - points[i][0])
        i = i + 1
    #build lambdaVec
    i = 1
    while i < len(points) - 1:
        lambdaVec.append(hVec[i]/(hVec[i-1] + hVec[i]))
        i = i + 1
    # build muVec
    i = 1
    while i < len(points) - 1:
        muVec.append(1-lambdaVec[i-1])
        i = i + 1
    #build dVec
    i = 1
    while i < len(points) - 1:
        dVec.append((6/(hVec[i-1]+hVec[i]))*(((points[i+1][1]-points[i][1])/hVec[i])-(points[i][1]-points[i-1][1])/hVec[i-1]))
        i = i + 1
    #build matrix
    i = 0
    matrix = []
    while i < len(hVec) - 1:
        rowToAppend = []
        j = 0
        while j < len(hVec) - 1:
            if i==j:
                rowToAppend.append(2)
            elif i > j and abs(i-j) == 1:
                rowToAppend.append(muVec[i])
            elif i < j and abs(i-j) == 1:
                rowToAppend.append(lambdaVec[i])
            else:
                rowToAppend.append(0)
            j = j + 1
        matrix.append(rowToAppend)
        i = i + 1
    MVec = []
    MVec.append(0) #natural spline, 0 element equals 0
    temp = solutionOfMatrixByGaussElimination(matrix,dVec)
    MVec.extend(temp)
    MVec.append(0) #natural spline, last element equals 0
    i = 0
    found = False
    while (i < len(points) - 1):
        if (points[i][0] - XToFind) * (points[i + 1][0] - XToFind) < 0:  # if this is true, X coordinate to find is between 2 these points
            found = True
            result = (math.pow(points[i + 1][0] - XToFind,3)*MVec[i]+math.pow(XToFind - points[i][0],3)*MVec[i+1])/(6*hVec[i])+(((points[i + 1][0] - XToFind)*points[i][1])+(XToFind - points[i][0])*points[i+1][1])/hVec[i]-((points[i + 1][0] - XToFind)*MVec[i]+(XToFind - points[i][0])*MVec[i+1])*hVec[i]/6
            break
        i = i + 1
    if not found:
        return None
    else:
        return result


def buildIdentityMatrix(n):
    """
    Builds and returns Identity matrix nxn. if n is not integer returns None
    :param n: size of matrix (nxn)
    :return: Identity matrix nxn
    """
    if not isinstance(n, int):
        return None
    resmat = []
    i = 0

    while i < n:
        resrow = []
        j = 0
        while j < n:
            if i == j:
                resrow.append(1)
            else:
                resrow.append(0)
            j = j + 1
        resmat.append(resrow)
        i = i + 1
    return resmat


def isMatrix(m):
    """
    Returns True if m is matrix and False otherwise
    :param m: parameter to be checked if it is matrix
    :return: True if m is matrix and False otherwise
    """
    try:
        # check if input are matrixes
        iter(m)
        if len(m) < 1:
            return False  # matrix is empty
        matcolcount = len(m[0])
        for i in m:
            iter(i)
            if len(i) != matcolcount:
                return False  # columns of matrix have different number of elements
    except TypeError as e:
        return False  # input is not matrixes
    return True


def multiplyMatrix(matrix1, matrix2):
    """
    Return the result of multiplication matrix1 and matrix2
    :param matrix1: first matrix
    :param matrix2: second matrix
    :return: the result of multiplication matrix1 and matrix2
    """
    # check if input is correct
    try:
        # check if input are matrixes
        if not isMatrix(matrix1) or not isMatrix(matrix2):
            return None  # one or both inputs are not matrixes
        iter(matrix1)
        iter(matrix2)
        mat1colcount = len(matrix1[0])
        mat1rowcount = len(matrix1)
        mat2colcount = len(matrix2[0])
        mat2rowcount = len(matrix2)
        if (mat1colcount != mat2rowcount):
            return None  # can't multiply matrixes

        # check if matrixes can be multiplied
    except TypeError as e:
        return None  # input is not matrixes
    k = 0
    resmat = []
    while k < mat1rowcount:
        i = 0
        resrow = []
        while i < mat2colcount:
            sum = 0
            j = 0
            while j < mat1colcount:
                sum = sum + matrix1[k][j] * matrix2[j][i]
                j = j + 1
            resrow.append(sum)
            i = i + 1
        resmat.append(resrow)
        k = k + 1
    return resmat


def solutionOfMatrixByGaussElimination(matrix, resvec):
    """
    Return the solution vector of square matrix that has result vector resvec and also print the elementary matrixes that was using during Gauss elimination
    :param matrix: matrix
    :param resvec: result vector
    :return: solution vector of square matrix that has result vector resvec
    """
    if not isMatrix(matrix):
        return []
    try:
        iter(resvec)
    except TypeError as e:
        return []  # resvec is not vector
    rowscount = len(matrix)
    if rowscount != len(resvec):
        return []  # number of rows in matrix doesn't match number of rows in result vector
    elementaryMatrixesList = []
    i = 0
    while i < rowscount:
        colscount = len(matrix[i])
        if colscount != rowscount:
            return []  # matrix is not squared
        rowdivider = matrix[i][i]
        if rowdivider == 0:
            s = i + 1
            flag = 0
            while s < rowscount:
                if matrix[s][i] != 0:
                    # swap lines
                    temp = matrix[s].copy()
                    matrix[s] = matrix[i].copy()
                    matrix[i] = temp
                    flag = 1
                    rowdivider = matrix[i][i]
                    temp = resvec[s]
                    resvec[s] = resvec[i]
                    resvec[i] = temp
                    break
                s += 1
            if flag == 0:
                return []  # same lines in matrix
        elementaryMatrix = buildIdentityMatrix(rowscount)
        elementaryMatrix[i][i] /= rowdivider
        matrix = multiplyMatrix(elementaryMatrix, matrix)
        elementaryMatrixesList.insert(0, elementaryMatrix)
        resvec[i] = resvec[i] / rowdivider
        k = 0
        while k < rowscount:
            if k != i:
                n = matrix[k][i]
                elementaryMatrix = buildIdentityMatrix(rowscount)
                elementaryMatrix[k][i] = -n * matrix[i][i]
                matrix = multiplyMatrix(elementaryMatrix, matrix)
                elementaryMatrixesList.insert(0, elementaryMatrix)
                resvec[k] = resvec[k] - n * resvec[i]
            k += 1
        i += 1
    counter = len(elementaryMatrixesList)
    for i in elementaryMatrixesList: #print the elementary matrixes list
        #print("E" + str(counter) + ": " + str(i))
        counter -= 1
    return resvec



#main
print("Linear interpolation:")
result = linearInterpolation(points,XToFind)
if result is None:
    print("Can't find the result by linear interpolation because X to find is not between points in table")
else:
    print("The result is: " + str(result))

print("Polynomial interpolation:")
result = polynomialInterpolation(points,XToFind)
print("The result is: " + str(result))

print("Lagrange interpolation:")
result = lagrangeInterpolation(points,XToFind)
print("The result is: " + str(result))

print("Nivul interpolation:")
result = nivulInterpolation(points,XToFind)
print("The result is: " + str(result))


print("Spline interpolation:")
result = splineInterpolation(points, XToFind)
if result is None:
    print("Can't find the result by spline method interpolation because X to find is not between points in table")
else:
    print("The result is: " + str(result))