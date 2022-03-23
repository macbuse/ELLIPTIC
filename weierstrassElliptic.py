from mpmath import jtheta, pi, exp, sqrt, polyroots, agm, log

def p_weierstrass_from_g2_g3(g2, g3, derivative = 0):
    if derivative < 0 or derivative > 3 or not isinstance(derivative, int):
        raise ValueError("`derivative` must be an integer between 0 and 3.")
    r1, r2, r3 = polyroots([4, 0, -g2, -g3])
    e3 = r3
    a = sqrt(r1 - r3)
    b = sqrt(r1 - r2)
    c = sqrt(r2 - r3)
    w1 = None
    if abs(a + b) < abs(a - b):
        b *= -1
    if abs(a + c) < abs(a - c):
        c *= -1
    if abs(c + 1j*b) < abs(c - 1j*b):
        e3 = r1
        a = sqrt(r3 - r1)
        b = sqrt(r3 - r2)
        c = sqrt(r2 - r1)
        w1 = 1 / agm(1j*b, c) 
    else:
        w1 = 1 / agm(a, b)
        
    w3 = 1j / agm(a, c)
    q = exp(1j * pi * w3/w1)
    if derivative != 1:
        pw0 = lambda z: (
            e3 + (pi * jtheta(2, 0, q) * jtheta(3, 0, q) * jtheta(4, z/w1, q)
              / (pi * w1 * jtheta(1, z/w1, q)))**2    
        )
        if derivative == 0:
            return pw0
        if derivative == 2:
            return lambda z: 6*pw0(z)**2 - g2/2
    f = (jtheta(1, 0, q, 1)**3 
        / (jtheta(2, 0, q) * jtheta(3, 0, q) * jtheta(4, 0, q)))
    pw1 = lambda z: (
        -2*(1/w1)**3 * jtheta(2, z/w1, q) * jtheta(3, z/w1, q) 
            * jtheta(4, z/w1, q) * f / jtheta(1, z/w1, q)**3
    )
    if derivative == 1:
        return pw1
    if derivative == 3:
        return lambda z: 12 * pw0(z) * pw1(z)  



def p_weierstrass_from_w1_w2(w1, w2):
    print(w1,w2, w2/w1)
    if w2.imag*w1.real < w1.imag*w2.real:
        raise ValueError("No solution. Do you want to exchange `w1 and `w2`?")
    ratio = w2 / w1
    if ratio.imag <= 0:
        raise ValueError("The ratio `w2/w1` must have a positive imaginary part.")
    q = exp(1j * pi * ratio)
    j2 = jtheta(2, 0, q)
    j3 = jtheta(3, 0, q)
    g2 = 4/3 * (pi/2/w1)**4 * (j2**8 - (j2*j3)**4 + j3**8) 
    g3 = 8/27 * (pi/2/w1)**6 * (j2**12 - (
        (3/2 * j2**8 * j3**4) + (3/2 * j2**4 * j3**8) 
    ) + j3**12)
    return p_weierstrass_from_g2_g3(g2, g3)


def p_weierstrass_from_tau(tau):
    return p_weierstrass_from_w1_w2(1.0, tau)

if __name__ == "___main__":
    z = 1+1j
    w1 = 1.0 + 4.0j
    w2 = 1.0 + 15.0j
    print(p_weierstrass_from_w1_w2(w1, w2)(z))
    # (0.0076936368424553 - 0.498821838483149j)

    print(p_weierstrass_from_tau(1 + 4j)(z))
    #(-0.430565782630798 - 3.62705469588323e-16j)

