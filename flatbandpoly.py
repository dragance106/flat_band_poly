from sympy import *

def fbpoly(*args):
    """
    args is a tuple where the first entry denotes the matrix describing the intracell links in a lattice,
    while the following entries denote intercell link amplitudes and their types.
    To simplify the input for various possibilities of finite-range hopping,
    each intercell matrix is followed by the hopping type
    which is encoded by a combination of the dimensional variables z1...zk.
    For example, after the intercell matrix Matrix([[0,1],[0,0]])
    we may specify z1*z2 to denote that it represents
    simultaneous one step hoppings in both the first (z1) and the second (z2) dimension.
    Or we may specify z1**2 to denote that it represents
    two step hopping along the first (z1) dimension.
    Or we may simply use z1 or z2 to denote the one step hopping
    along the first or the second dimension.

    Do note that the Bloch Hamiltonian is created directly from these arguments
    in the form

    H_B = args[0] + args[1] * args[2] + args[1].H (or .T) * args[2]^(-1)
                  + args[3] * args[4] + args[3].H (or .T) * args[4]^(-1)
                  + ...

    so that you can also use other variables to denote directions (such as x, y or z, for example).
    However, to obtain valid and interpretable results
    you should necessarily use exponents to denote several steps along one dimension
    and the products to denote simultaneous steps along different dimensions
    (such as x**3 or x**2*y).

    The method rewrites the characteristic polynomial det(alpha I - H_B) of the Bloch Hamiltonian
    in terms of monomials of the intercell hopping variables, and
    then returns the list of these coefficients.
    Since the coefficients are polynomials in alpha by themselves,
    alpha should not be used as one of the intercell hopping variables.

    :param args: The tuple of arguments
                 the first one being the matrix that describes the intracell links,
                 the remaining ones come in pairs of
                 the matrix that denotes the intercell hopping amplitudes and
                 the variable symbol that describes the type of the hopping step.
    :return: polynomial coefficients in alpha when the characteristic polynomial of the Bloch Hamiltonian
             det(alpha I - H_B) is rewritten in terms of monomials of the intercell hopping variables.
    """

    # Constructing the Bloch Hamiltonian out of the intracell matrix
    # and the pairs of (intercell matrix, hopping encoding variable)
    hb = args[0]
    dimvars = set()         # empty set of dimensional variables
    for i in range(1, len(args), 2):
        # should use .H below if not all entries of args[i] are real numbers, otherwise .T is fine
        hb = hb + args[i] * args[i+1] + args[i].T / args[i+1]
        # pick up the dimensional variables that appear among the arguments
        dimvars = dimvars.union(args[i+1].free_symbols)

    # Compute the characteristic polynomial in the ring of variables from args[2], args[4], ...
    # and then convert it into an ordinary expression
    cp = hb.charpoly(x='alpha').as_expr()

    print(f'original charpoly was: {cp}')

    # Get rid of the negative exponents by multiplying with the denominator
    cp *= fraction(together(cp))[1]

    # Convert it into a multivariate polynomial in terms of dimensional variables
    cp = Poly(cp, *dimvars)

    # Return the obtained multivariate polynomial
    return cp


def do_kagome():
    # kagome lattice
    z1, z2 = symbols('z1 z2')
    p = fbpoly(Matrix([[0,1,1],[1,0,1],[1,1,0]]),
                Matrix([[0,0,0],[1,0,0],[0,0,0]]), z1,
                Matrix([[0,0,0],[0,0,0],[1,0,0]]), z2,
                Matrix([[0,0,0],[0,0,0],[0,1,0]]), z1**(-1)*z2)

    coeffs = p.coeffs()
    g = gcd_list(coeffs)

    print(f'the flat band polynomial is {p};')
    print(f'its coefficients are {coeffs};')
    print(f'and their gcd is {g}.')


def do_octahedron():
    B = zeros(6,6)
    # two cones
    for i in range(1,5):
        B[0,i]=1; B[i,0]=1;
        B[5,i]=1; B[i,5]=1;
    # the square
    B[1,2]=1; B[2,1]=1
    B[2,3]=1; B[3,2]=1
    B[3,4]=1; B[4,3]=1
    B[4,1]=1; B[1,4]=1

    w = symbols('w')
    Aw = zeros(6,6)
    Aw[5,0] = 1

    p = fbpoly(B, Aw, w)
    coeffs = p.coeffs()
    g = gcd_list(coeffs)

    print(f'the flat band polynomial is {p};')
    print(f'its coefficients are {coeffs};')
    print(f'and their gcd is {g}.')


def do_dice_twisted():
    z1, z2 = symbols('z1 z2')
    p=fbpoly(Matrix([[0,1,1],[1,0,1],[1,1,0]]),
              Matrix([[0,1,0],[0,0,0],[0,0,0]]), z1,
              Matrix([[0,0,0],[1,0,0],[1,0,0]]), z2,
              Matrix([[0,0,0],[0,0,0],[1,0,0]]), z1*z2)

    coeffs = p.coeffs()
    g = gcd_list(coeffs)

    print(f'the flat band polynomial is {p}.')
    print(f'its coefficients are:')
    for ct in coeffs:
        print(f'{ct}')
    print(f'their gcd is {g}')


def do_kagome_system():
    # kagome lattice with unknown edge weights - x for intracell, y for intercell links
    z1, z2 = symbols('z1 z2')
    x, y = symbols('x y')
    p = fbpoly(Matrix([[0,x,x],[x,0,x],[x,x,0]]),
                Matrix([[0,0,0],[y,0,0],[0,0,0]]), z1,
                Matrix([[0,0,0],[0,0,0],[y,0,0]]), z2,
                Matrix([[0,0,0],[0,0,0],[0,y,0]]), z1**(-1)*z2)

    coeffs = p.coeffs()
    g = gcd_list(coeffs)

    print(f'the flat band polynomial is {p}.')
    print(f'its coefficients are {coeffs}')
    print(f'their gcd is {g}.')

    # solve the systems of equations
    sols = solve(coeffs, [x, y], dict=True)
    print(f'solutions in general: {sols}')


def do_dice_system():
    z1, z2 = symbols('z1 z2')
    a, b, c, d, e, f, g = symbols('a b c d e f g')
    p = fbpoly(Matrix([[0,a,b],[a,0,c],[b,c,0]]),
              Matrix([[0,f,0],[0,0,0],[0,0,0]]), z1,
              Matrix([[0,0,0],[d,0,0],[e,0,0]]), z2,
              Matrix([[0,0,0],[0,0,0],[g,0,0]]), z1*z2)

    coeffs = p.coeffs()
    g = gcd_list(c)

    print(f'the char.poly is {p}.')
    print(f'its coefficients are:')
    for ct in coeffs:
        print(f'{ct}')
    print(f'their gcd is {g}')

    alpha = symbols('alpha')
    newcoeffs = [ct.subs([(c,0)]) for ct in coeffs]
    print(f'after substituting c=0:')
    for ct in newcoeffs:
        print(f'{ct}')

    sols = solve(newcoeffs, [a, b, d, e, f, g], dict=True)
    print(f'solutions for general alpha: {sols}')


def do_lieb_extension():
    # extended lieb lattice
    z1, z2 = symbols('z1 z2')
    t1, t2, t3, ea, eb, ec, ed = symbols('t1 t2 t3 ea eb ec ed')
    cp = fbpoly(Matrix([[ea,t1,t1,0],[t1,eb,t2,t3],[t1,t2,ec,t3],[0,t3,t3,ed]]),
                Matrix([[0,0,0,0],[t1,0,0,0],[0,0,0,0],[0,0,0,0]]), z1,
                Matrix([[0,0,0,0],[0,0,0,0],[t1,0,0,0],[0,0,0,0]]), z2)

    coeffs = cp.coeffs()

    # print(f'the char.poly is {cp}.')
    # print(f'its coefficients are:')
    # for c in coeffs:
    #     print(f'{c}')

    # substitute a particular value for alpha
    alpha = symbols('alpha')
    newcoeffs = [c.subs([(t3, sqrt(-t2*(alpha-ed))),
                         (eb, ec),
                         (t2, -alpha+ec)]) for c in coeffs]

    print(f'newcoeffs are:')
    for c in newcoeffs:
        print(f'{expand(c)}')

    # solve the systems of equations
    # sols = solve(coeffs, [t2, t3], dict=True)
    # sols_at_minus_two = solve(newcoeffs, [x, y], dict=True)

    # print(f'solutions at alpha={alpha}: {sols}')
    # print(f'solutions at alpha=-2: {sols_at_minus_two}')


def do_diamond_general():
    z = symbols('z')
    a, b, c, d, e, f = symbols('a b c d e f')
    cp = fbpoly(Matrix([[0,a,b,0],[a,0,c,d],[b,c,0,e],[0,d,e,0]]),
                Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[f,0,0,0]]), z)

    coeffs = cp.coeffs()

    print(f'the char.poly is {cp}.')
    print(f'its coefficients are:')
    for c in coeffs:
        print(f'{c}')

    # substitute a particular value for alpha
    alpha = symbols('alpha')
    newcoeffs = [c.subs(alpha, -1) for c in coeffs]

    print(f'newcoeffs are:')
    for c in newcoeffs:
        print(f'{expand(c)}')
    print(f'their gcd is: {gcd_list(newcoeffs)}')

    # solve the systems of equations
    sols = solve(newcoeffs, [a, b, c, d, e, f], dict=True)
    # sols_at_minus_two = solve(newcoeffs, [x, y], dict=True)

    print(f'solutions at alpha={alpha}: {sols}')
    # print(f'solutions at alpha=-2: {sols_at_minus_two}')



if __name__=="__main__":
    # do_kagome()
    # do_octahedron()
    # do_dice_twisted()
    # do_kagome_system()
    # do_dice_system()

    # do_lieb_extension()
    # do_diamond_general()
    # do_diamond()

    # uncomment one of the above lines to run the appropriate method
    pass
