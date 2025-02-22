# Méthode du gradient pour les fonctions QUADRATIQUES. Sinon la fonction gradf ne représente rien
# Et surtout ne représente pas le gradient.
def gradientMethods(x0, b, error, maxIt, method, prodMatrix):  
    gradf = lambda x : prodMatrix(x)-b
  
    if method == "pasFixe":
        iterFunction = lambda x : x-0.25*gradf(x)
    elif method == "pasOptimal":
        # Formule donnée par le cours pour le pas optimal
        pasOptimal = lambda x : (gradf(x).transpose() * gradf(x)) / (gradf(x).transpose() * prodMatrix(gradf(x)))
        # On recalcule le pas optimal au point courant à chaque itération.
        # C'est un matrice (1, 1) donc on doit prendre sa seule composante car sage ne gère pas l'opération sinon.
        iterFunction = lambda x : x-pasOptimal(x)[0,0]*gradf(x)

    i = 0
    while(i < maxIt and (prodMatrix(x0)-b).norm(2) > error):
        x0 = iterFunction(x0)
        i = i + 1
        
    return (x0, i)

def testMethods(dimension, error, maxIt, method, textOutput):
    # Vecteur initial
    x0 = Matrix(RR, [[0] for i in range(dimension)])
    # Second membre
    b = Matrix(RR, [[i] for i in range(dimension)])

    def prodA(u):
        v = Matrix(RR,dimension,1)
        v[0,0] = 2.0*u[0][0] - u[1][0]
        for i in range(1,dimension-1):
            v[i,0] = -u[i-1][0] +2.0*u[i][0] - u[i+1][0]
        v[dimension-1,0] = -u[dimension-2][0] +2.0*u[dimension-1][0]
        return v
    
    sol, itNumber = gradientMethods(x0, b, error, maxIt, method, prodA)
    if textOutput == True:
        print(
            sol, 
            "\n\n", "La solution est valide: {bool} \n La solution est atteinte en {count} itérations.".format(
                bool = True if (prodA(sol) - b).norm() < 10^(-5) else False, 
                count = itNumber
            )
        )
    return itNumber

print('<<============= GRADIENT PAS FIXE ==============>\n')
testMethods(10, 10^(-3), 10000, "pasFixe", textOutput = True)
print('\n<<============= GRADIENT PAS OPTIMAL ==============>\n')
testMethods(10, 10^(-3), 10000, "pasOptimal", textOutput = True)


print('\n<<============= COMPARAISON DES METHODES ==============>\n')
L1 = [testMethods(i, 10^(-3), 5000, "pasFixe", textOutput = False) for i in [10, 20, 30, 40, 50]]
L2 = [testMethods(i, 10^(-3), 5000, "pasOptimal", textOutput = False) for i in [10, 20, 30, 40, 50]]

print(L1, '\n', L2)