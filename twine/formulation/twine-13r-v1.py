from sage import all


def compute_double_dlct_1(sa):
    """
    Compute the double differential-linear connectivity table for the S-box of TWINE
    """

    n = sa.input_size()
    dlct = sa.differential_linear_connectivity_table()
    lat = sa.linear_approximation_table(scale='correlation')
    double_dlct = [[0 for j in range(2**n)] for i in range(2**n)]
    for D in range(2**n):
        for L in range(2**n):
            for LM in range(2**n):
                double_dlct[D][L] += (dlct[D][LM] / 2**(n - 1)) * (lat[LM, L]**2)
    return double_dlct

def compute_double_dlct_2(sa):
    """
    Compute the double differential-linear connectivity table for the S-box of TWINE
    """

    n = sa.input_size()
    dlct = sa.differential_linear_connectivity_table()
    ddt = sa.difference_distribution_table()
    double_dlct = [[0 for j in range(2**n)] for i in range(2**n)]
    for D in range(2**n):
        for L in range(2**n):
            for DM in range(2**n):
                double_dlct[D][L] += (ddt[D][DM]/2**n) * (dlct[DM][L] / 2**(n - 1))
    return double_dlct


def check_two_sides(n, dlct, lat, ddt):
    ddlct_1 = [[0 for _ in range(2**n)] for _ in range(2**n)]
    ddlct_2 = [[0 for _ in range(2**n)] for _ in range(2**n)]
    lat = [[2**n * lat[i][j] for j in range(2**n)] for i in range(2**n)]
    for D in range(2**n):
        for L in range(2**n):
            for DM in range(2**n):
                ddlct_1[D][L] += ddt[D][DM] * (dlct[DM][L])
    for D in range(2**n):
        for L in range(2**n):
            for LM in range(2**n):
                ddlct_2[D][L] += dlct[D][LM] * (lat[LM][L]**2) * 2**(-n)
    if ddlct_1 != ddlct_2:
        raise ValueError("Mismatch")
    else:
        print("Match")




def compute_corr_13r_v0(D, L, double_dlct, lat):
    """
    Compute the correlation for 13-round DLD-v0 
    """

    corr = 0
    for LM in range(2**n):
        corr += (double_dlct[D][LM]) * (lat[LM, L]**2)
    return corr

def compute_corr_13r_v1(D1, L4, ddt, dlct):
    """
    Compute the correlation for 13-round DLD-v1
    """

    corr = 0
    for D2 in range(2**n):
        for D3 in range(2**n):
            corr += ddt[D1][D2] * ddt[D2][D3] * dlct[D3][L4]
    corr = corr/(2**(3*n - 1))
    return corr

def compute_corr_13r_v2(D, L, double_dlct, ddt):

    corr = 0
    for DM in range(2**n):
        corr += ddt[D][DM] * double_dlct[DM][L]
    corr = corr/2**(n)
    return corr

if __name__ == '__main__':
    from sage.crypto.sboxes import TWINE as sb
    from sboxanalyzer import *
    sa = SboxAnalyzer(sb)    
    n = sa.input_size()    
    ddt = sa.difference_distribution_table()
    lat = sa.linear_approximation_table(scale='correlation')
    dlct = sa.differential_linear_connectivity_table()
    double_dlct = compute_double_dlct_1(sa)
    double_dlct_temp = compute_double_dlct_2(sa)
    check_two_sides(n, dlct, lat, ddt)
    if double_dlct != double_dlct_temp:
        raise ValueError("Double DLCT table mismatch")
    dict_of_results = {} 

    for D in range(2**n):
        for L in range(2**n):
            corr1 = abs(compute_corr_13r_v0(D, L, double_dlct, lat))         
            corr2 = abs(compute_corr_13r_v1(D, L, ddt, dlct))     
            corr3 = abs(compute_corr_13r_v2(D, L, double_dlct, ddt))
            if corr1 != corr2:
                raise ValueError("Correlation mismatch")
            if corr1 != corr3:
                raise ValueError("Correlation mismatch")
            if corr1 > 0:
                abscorr1 = float(log(corr1**2, 2))
                print("({}, {}): Corr = 2^{:0.2f}".format(hex(D), hex(L), abscorr1))
                print("#"*30)
                dict_of_results[(D, L)] = abscorr1
            else:
                print("({}, {}): Corr = 0".format(hex(D), hex(L)))
                print("#"*30)

    # D = 0xd
    # L = 0xa
    # corr1 = abs(compute_corr_13r_v0(D, L, double_dlct, lat))
    # sn = " - " if corr1 < 0 else " + "
    # abscorr1 = float(log(corr1**2, 2))
    # print("({}, {}): Corr = {}2^{:0.2f}".format(hex(D), hex(L), sn, abscorr1))


    # dict_values = dict_of_results.values()
    # dict_values_except_zero = [v for v in dict_values if v != 0]
    # max_value = max(dict_values_except_zero)
    # max_keys = [k for k, v in dict_of_results.items() if v == max_value]
    # print("Maximum correlation: {}".format(max_value))
    # print("The input-output pairs with maximum correlation are:")
    # for k in max_keys:
    #     print(k)
