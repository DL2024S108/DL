from sage import all

def compute_ddt4(sb):
    """
    compute DDT4
    """

    n = sb.input_size()
    ddt = sb.difference_distribution_table()
    ddt4 = [[0 for _ in range(2**n)] for _ in range(2**n)]
    for D1 in range(2**n):
        for D2 in range(2**n):
            for D3 in range(2**n):
                for D4 in range(2**n):
                    for D5 in range(2**n):
                        ddt4[D1][D5] += ddt[D1][D2] * ddt[D2][D3] * ddt[D3][D4] * ddt[D4][D5]
    return ddt4

def compute_corr_13r_v0(ddt4, dlct, ddt, D1, L):
    """
    Compute the correlation for 13-round DLD-v0 
    """

    corr = 0
    for D2 in range(2**n):
        for D6 in range(2**n):
            for D7 in range(2**n):
                corr += ddt[D1][D2] * ddt4[D2][D6] * ddt[D1][D7] * dlct[D6 ^ D7][L] * dlct[D2][L]
    corr = corr/(2**(8*n - 2))
    return corr

if __name__ == '__main__':
    from sage.crypto.sboxes import TWINE as sb
    from sage.crypto.sboxes import SBox
    from sboxanalyzer import *
    sa = SboxAnalyzer(sb)    
    n = sa.input_size()    
    ddt = sa.difference_distribution_table()
    lat = sa.linear_approximation_table(scale='correlation')
    dlct = sa.differential_linear_connectivity_table()
    ddt4 = compute_ddt4(sa)

    for D in range(2**n):
        for L in range(2**n):
            corr1 = abs(compute_corr_13r_v0(ddt4, dlct, ddt, D, L))
            if corr1 > 0:
                abscorr1 = float(log(corr1, 2))
                print("({}, {}): Corr = 2^{:0.2f}".format(hex(D), hex(L), abscorr1))
                print("#"*30)
            else:
                print("({}, {}): Corr = 0".format(hex(D), hex(L)))
                print("#"*30)
