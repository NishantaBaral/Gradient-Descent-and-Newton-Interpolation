import numpy as np

#constructing a divided difference table following Algorithm 7.3 in the book to a tee. Nothing fancy here.
def divided_difference(x,y):
    n = len(x)
    
    table = np.empty((n,n))
    table[:,0] = y

    for k in range(1,n):
        for i in range (n-k):
            table[i,k] = (table[i+1,k-1] - table[i,k-1]) / (x[i+k]-x[i])
   
    return table[0,:]

#Scaled chevbyshv nodes and reversing the order to get em in increasing order.
#Just a note, there are two formulas for chebyshev nodes, one is the one used here and the other is x_k = scaling terms*cos((2k-1)/(2n)*pi)
def chebyshev_nodes(a,b,n):
    nodes = []
    for k in range(n):
        node = 0.5*(a+b) + 0.5*(b-a)*np.cos((2*k+1)/(2*n+2)*np.pi)
        nodes.append(node)
    return nodes[-1::-1]

#newton interpolation following Algorithm 7.2 in the book. The pseudocode has \prod_{j=0}^{i-1} but in python that changes to loop for j in range 0 to i
def newton_interpolation(x,y,x0):
    coeff = divided_difference(x,y)
    polynomial = coeff[0]

    for i in range(1,len(coeff)):
        term = coeff[i]
        for j in range(0,i):
            term *= (x0 - x[j])
        polynomial += term
    return polynomial


if __name__ == "__main__":
    #Doing the interpolation for f(x) = log(1+x) at x0 = 1 using 14 Chebyshev nodes in the interval [0, e-1]
    nodes = chebyshev_nodes(0,np.exp(1)-1,14)
    values = np.log(1+np.array(nodes))
    x0 = 1
    result = newton_interpolation(nodes, values, x0)
    print(f"Interpolated value at x={x0}: {result}")
    error = np.abs(np.log(1+x0) - result)
    print(f"Error: {error}")
