from typing import *
import sage.all as sagemath
import numpy as num

x = sagemath.var("x")
y = sagemath.var("y")

## <<==============================>> Utility <<================================ ##
def approx(x):
    return sagemath.numerical_approx(x, digits=4)
def diff(f, x):
    return sagemath.diff(f(x), x)
def vector(L):
    return sagemath.vector(L)

## <<==============================>> Functions Methods <<================================ ##

def numericalRootFinder(interval : tuple[int, int], f : Callable[[float], float], method : str, n : int):
    """
    Available methods : {"Dichotomy", "FixedPoint", "NewtonRaphson"}
    """
    (a, b) = (interval[0], interval[1])
    if f(a) * f(b) > 0: raise Exception("No root in initial interval.")
    
    def dichotomyMethod(a, b, f, n):
        i = 0
        
        # Tant que l'erreurr est trop grande, on divise les sous intervalles en deux et on choisit l'intervalle qui contient la racine
        while i <= n:
            middlePoint = (b + a) / 2
            i += 1
            
            if f(a) * f(middlePoint) <= 0:
                b = middlePoint
            elif f(b) * f(middlePoint) <= 0:
                a = middlePoint
        
        # L'intervalle final (a, b) contient la racine et on choisit (a+b)/2 comme approximation
        return approx((a + b) / 2)
    def fixedPointMethod(a, b, f, n):
        # On reformule la recherche des racines de f en terme de recherche de points fixe de g
        def g(x):
            return f(x) + x
        
        # Point initial
        x0 = (b + a) / 2
        
        # Constante de contraction de f
        derivative = diff(f(x), x)
        k = find_local_maximum(abs(derivative), a, b)[0] # type: ignore
        if k > 1:
            raise Exception("Function is not a contraction.")
        
        i = 0
        # Tant que l'erreur est trop grande, on réapplique g à x0
        while i <= n:
            i += 1
            x0 = g(x0)
        
        return approx(x0)
    def newtonRaphsonMethod(a, b, f, n):
        # Point initial
        x0 = (b + a) / 2
        
        derivative = diff(f(x), x)
        i = 0
        # Tant que l'erreur est trop grande, on trace la tangente et le nouveau point est l'intersection avec Ox
        while i <= n:
            i += 1
            x0 -= f(x0) / derivative(x=x0)
            if n > 100:
                raise Exception("Method did not converge after the maximum number of iterations allowed (100).")
        return approx(x0)
    
    methods = {"Dichotomy": dichotomyMethod, "FixedPoint": fixedPointMethod, "NewtonRaphson": newtonRaphsonMethod, }
    
    return methods[method](a, b, f, n)
def numericalIntegration(interval : tuple[int, int], f : Callable[[float], float], method : str, n : int):
    """
    Available methods : {"MiddlePoint", "Trapezoid", "Simpson"}
    """

    def middlePoint(f, a, b):
        return (b - a) * f((a + b) / 2)
    def trapezoid(f, a, b):
        return (b - a) * ((f(a) + f(b)) / 2)
    def simpson(f, a, b):
        return (b - a) * (f(a) + 4 * f((a + b) / 2) + f(b)) / 6
    
    methods = {
        "MiddlePoint": middlePoint, 
        "Trapezoid": trapezoid, 
        "Simpson": simpson, 
    }
    
    # Crée une subdivision régulière d e [a, b] et approxime l'intégrale par une méthode donnée sur tout les sous-segments
    def compositeMethode(interval, f, method, n):
        points = [interval[0] + k * (interval[1] - interval[0]) / n for k in range(n + 1)]
        area = 0
        for i in range(len(points) - 1):
            area += method(f, points[i], points[i + 1])
        return area
    
    return compositeMethode(interval, f, methods[method], n)
def interpolation(f : Callable[[float], float], points : Tuple[int, ...], method : str):
    """
    Available methods : {"Newton", "Lagrange"}
    """

    def lagrangeBasis(liste):  # Retourne la base de Lagrange associée aux points de L
        basis = []
        for i in range(len(liste)):
            P = 1
            for j in range(len(liste)):
                if i != j:
                    P *= (x - liste[j]) / (liste[i] - liste[j])
            basis.append(P)
        return basis
    def lagrangeMethod(f, points):
        return sum([f(points[i]) * lagrangeBasis(points)[i] for i in range(len(points))]) 
    
    # Différences divisées récursivement
    def dividedDiff(liste, f):
        n = len(liste)
        if n == 1:
            return f(liste[0])
        else:
            L1 = liste[-n + 1:]
            L2 = liste[0:n - 1]
        return (dividedDiff(L1, f) - dividedDiff(L2, f)) / (liste[n - 1] - liste[0])
    def newtonPolynomial(liste):     # On calcule le n-ième polynome de la base de Newton associée aux points
        polynomial = 1
        
        if len(liste) == 0: return 1
        for i in range(len(liste)):
            polynomial *= x - liste[i]
        return polynomial
    def newtonMethod(f, points):
        polynomial = 0
        
        for i in range(1, len(points) + 1):
            coeff = dividedDiff(points[0:i], f)     # On calcule la différence divisée associée à [x0, ..., xi]
            polynomial += coeff * newtonPolynomial(points[0:i - 1])     # Le polynome est egal a la f[x0, ..., xi] * N_xi
        return polynomial
    
    methods = {"Newton": newtonMethod, "Lagrange": lagrangeMethod}
    return methods[method](f, points)

## <<===============================>> Matrices Methods <<================================ ##

# Pivot de gauss
def pivotGauss(A):
    def findBiggestPivot(column, i):
        L = []
        for j in range(i, len(column)):
            L.append((j, abs(column[j])))
        return max(L, key=lambda x: abs(x[1]))

    for i in range(A.nrows()):  # Triangule la matrice
        bestPivotIndex = findBiggestPivot(A.column(i), i)[0]
        # Trouve le pivot le plus grand en valeur absolue de la colonne considérée

        if bestPivotIndex != i and i != A.nrows() - 1:
            A.swap_rows(i, bestPivotIndex)
        pivot = A[i, i]

        # Change le pivot courant pour le plus grand en valeur absolue, sauf pour la dernière ligne

        for j in range(i + 1, A.nrows()):
            coeff = -A[j, i] / pivot
            A[j] += coeff * A[i]
        # Elimine les valeurs en dessous du pivot

    for i in range(A.nrows()):  # Remonte le système
        currentPivotIndex = A.nrows() - i - 1
        A[currentPivotIndex] = A[currentPivotIndex] / float(
            A[currentPivotIndex, currentPivotIndex]
        )
        # Divise le coefficient diagonal pour avoir un 1 sur la diagonale.

        pivot = A[currentPivotIndex, currentPivotIndex]
        for j in range(currentPivotIndex):
            coeff = -A[j, currentPivotIndex] / pivot
            A[j] += coeff * A[currentPivotIndex]
        # Elimine les valeurs au dessus du pivot

    return A.column(A.nrows())

## <==================================> EDOs Methods <<================================== ##
   
# Take a vector field, and approximate a flow line given initial conditions, and over time t
def eulerMethod(F, time, stepSize, *initialConditions): 
   n = num.floor(time / stepSize)     # Number of iterations necessary for achieving time t

   L = []
   P0 = vector(initialConditions)
   for i in range(n):
      L.append(P0)
      P0 = vector(tuple(approx(P0[i], digits=4) for i in range(len(initialConditions))))  
      P0 = P0 + stepSize*F(*P0) # gamma(t0 + h) ~ gamma(t0) + hgamma'(t0) = gamma(t0) + hF(gamma(t0))
   return L

def lotkaVolterra(a, b, c, d): return (x*(a-b*y), y*(-c + d*x))
  
# 2D Pendulum data for modeling project
H = function('H')(x, y)
norm = function('norm')(x, y)
norm(x, y) = sagemath.sqrt(x**2+y**2)
H(x, y) = (y/norm(H(x, y)), -sagemath.sin(x)/norm(H(x, y)))

def pointListToPgfPlot(L, outputPath): # Procedure converting a list of points to a pgfplot readable data file
   with open(outputPath, 'w') as f:
      f.write('x y')
      f.write('\n')
      for line in L:
         f.write('{xvalue} {yvalue}'.format(xvalue = line[0], yvalue = line[1]))
         f.write('\n')
def computePendulumData():
   # Précision presque idéale 0.0035
   acc = 0.1
   # Trajectoires apériodiques
   aperiodicFlowPoints1 = eulerMethod(H, 12.5, acc, -3, 2) # Longueur d'arc 12.5
   aperiodicFlowPoints2 = eulerMethod(H, 13, acc, -3, 1) # Longueur d'arc 13
   aperiodicFlowPoints3 = eulerMethod(H, 13, acc, 9, -1) # Longueur d'arc 13
   aperiodicFlowPoints4 = eulerMethod(H, 12.5, acc, 9, -2) # Longueur d'arc 12.5

   # Trajectoires périodiques
   periodicFlowPoints1 = eulerMethod(H, 12, acc, 2.14, 0) # Longueur d'arc 12
   periodicFlowPoints2 = eulerMethod(H, 7, acc, 1.14, 0) # Longueur d'arc 7
   periodicFlowPoints3 = eulerMethod(H, 12, acc, 4.14, 0) # Longueur d'arc 12
   periodicFlowPoints4 = eulerMethod(H, 7, acc, 5.14, 0) # Longueur d'arc 7 
   return [
      aperiodicFlowPoints1, aperiodicFlowPoints2, aperiodicFlowPoints3, aperiodicFlowPoints4, periodicFlowPoints1, periodicFlowPoints2, periodicFlowPoints3, periodicFlowPoints4
   ]
def sendFormattedPendulumData(L):
   pointListToPgfPlot(L[0], 'C:\\Users\\Archimedean\\Desktop\\Maths\\NumericalAnalysisProjects\\data\\aperiodic1.dat')
   pointListToPgfPlot(L[1], 'C:\\Users\\Archimedean\\Desktop\\Maths\\NumericalAnalysisProjects\\data\\aperiodic2.dat')
   pointListToPgfPlot(L[2], 'C:\\Users\\Archimedean\\Desktop\\Maths\\NumericalAnalysisProjects\\data\\aperiodic3.dat')
   pointListToPgfPlot(L[3], 'C:\\Users\\Archimedean\\Desktop\\Maths\\NumericalAnalysisProjects\\data\\aperiodic4.dat')

   pointListToPgfPlot(L[4], 'C:\\Users\\Archimedean\\Desktop\\Maths\\NumericalAnalysisProjects\\data\\periodic1.dat')
   pointListToPgfPlot(L[5], 'C:\\Users\\Archimedean\\Desktop\\Maths\\NumericalAnalysisProjects\\data\\periodic2.dat')
   pointListToPgfPlot(L[6], 'C:\\Users\\Archimedean\\Desktop\\Maths\\NumericalAnalysisProjects\\data\\periodic3.dat')
   pointListToPgfPlot(L[7], 'C:\\Users\\Archimedean\\Desktop\\Maths\\NumericalAnalysisProjects\\data\\periodic4.dat')