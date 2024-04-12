from sage.all import *
from sboxanalyzer import *
import pickle

from sage.crypto.sboxes import AES as sb
sa = SboxAnalyzer(sb)

try:
    with open('ddt.pkl', 'rb') as file:
        ddt = pickle.load(file)
    print("ddt.pkl exists")
except:
    print("ddt.pkl does not exist")
    print("Generating ddt.pkl")
    ddt = sa.difference_distribution_table()
    ddt = [[ddt[i][j] for j in range(256)] for i in range(256)]
    with open('ddt.pkl', 'wb') as file:
        pickle.dump(ddt, file)

try:
    with open('dlct.pkl', 'rb') as file:
        dlct = pickle.load(file)
    print("dlct.pkl exists")
except:
    print("dlct.pkl does not exist")
    print("Generating dlct.pkl")
    dlct = sa.differential_linear_connectivity_table()
    with open('dlct.pkl', 'wb') as file:
        pickle.dump(dlct, file)