from sage import all

def doct_product(a, b, n=4):
    """
    Compute the dot product of two inetgers as a list of bits
    """

    output = 0
    for i in range(n):
        output ^= ((a >> i) & 1) * ((b >> i) & 1)
    return output


def compute_corr_13r_v0(sa, D, L):
    """
    Compute the correlation for 13-round DLD-v0 
    """

    lat = sa.linear_approximation_table(scale='correlation')
    dlct = sa.differential_linear_connectivity_table()
    n = sa.input_size()
    corr = 0
    for LM in range(2**n):
        corr += (lat[LM, L]**2) * (dlct[D][LM])
    corr = corr/(2**(n - 1))
    return corr

def compute_corr_13r_v1(sa, D, L):
    """
    Compute the correlation for 13-round DLD-v1 
    """

    corr = 0
    ddt = sa.difference_distribution_table()
    dlct = sa.differential_linear_connectivity_table()
    for DM in range(2**n):
        corr += ddt[D, DM] * dlct[DM][L]
    corr = corr/(2**(2*n - 1))
    return corr

if __name__ == '__main__':
    from sage.crypto.sboxes import TWINE as sb
    from sboxanalyzer import *
    sa = SboxAnalyzer(sb)
    n = sa.input_size()

    for D in range(2**n):
        for L in range(2**n):
            corr1 = abs(compute_corr_13r_v0(sa, D, L))
            corr2 = abs(compute_corr_13r_v1(sa, D, L))
            if corr1 != corr2:
                raise ValueError("Correlation mismatch")
            if corr1 > 0:
                abscorr1 = float(log(corr1, 2))
                abscorr2 = float(log(corr2, 2))
                print("({}, {}): Correlation = 2^{:0.2f}".format(hex(D), hex(L), abscorr1))
                print("({}, {}): Correlation = 2^{:0.2f}".format(hex(D), hex(L), abscorr2))
                print("#"*10)
            else:
                print("({}, {}): Corr = 0".format(hex(D), hex(L)))
                print("#"*10)
    