#!/usr/env/bin python3
#-*- coding: UTF-8 -*-

"""
MIT License

Copyright (c) 2023 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import logging
from pathlib import Path
from random import randint
logging.basicConfig(filename="minizinc-python.log", level=logging.DEBUG)
import time
import minizinc
import datetime
from argparse import ArgumentParser, RawTextHelpFormatter
from draw import *
import copy


class DiffLin:
    DL_counter = 0

    def __init__(self, param) -> None:
        DiffLin.DL_counter += 1
        self.id = DiffLin.DL_counter
        self.name = "DiffLin" + str(self.id)
        self.type = "DiffLin"
        self.RU = param["RU"]
        self.RM = param["RM"]
        self.RL = param["RL"]
        self.RMU = param["RMU"]
        self.RML = param["RML"]
        self.offset = param["offset"]
        self.WU = param["WU"]
        self.WM = param["WM"]
        self.WL = param["WL"]   
        self.NC = param["NC"]
        self.RD = self.RU + self.RM + self.RL
        self.cp_solver_name = param["solver"]
        ##################################################
        # # Use this block if you install Or-Tools bundeled with MiniZinc
        # if self.cp_solver_name == "ortools":
        #     self.cp_solver_name = "com.google.ortools.sat"
        ##################################################        
        self.cp_solver = minizinc.Solver.lookup(self.cp_solver_name)
        self.time_limit = param["timelimit"]
        self.num_of_threads = param["np"]
        self.mzn_file_name = None
        self.output_file_name = param["output"]
        self.mzn_file_name = "attack.mzn"
    
    #############################################################################################################################################
    #############################################################################################################################################    
    #  ____                           _        __                        ____   _       _    _                       _       _                 
    # / ___|   ___   __ _  _ __  ___ | |__    / _|  ___   _ __    __ _  |  _ \ (_) ___ | |_ (_) _ __    __ _  _   _ (_) ___ | |__    ___  _ __ 
    # \___ \  / _ \ / _` || '__|/ __|| '_ \  | |_  / _ \ | '__|  / _` | | | | || |/ __|| __|| || '_ \  / _` || | | || |/ __|| '_ \  / _ \| '__|
    #  ___) ||  __/| (_| || |  | (__ | | | | |  _|| (_) || |    | (_| | | |_| || |\__ \| |_ | || | | || (_| || |_| || |\__ \| | | ||  __/| |   
    # |____/  \___| \__,_||_|   \___||_| |_| |_|   \___/ |_|     \__,_| |____/ |_||___/ \__||_||_| |_| \__, | \__,_||_||___/|_| |_| \___||_|   
    #                                                                                                  |___/                                   
    # Search for a distinguisher using MiniZinc

    def search(self):
        """
        Search for a distinguisher
        """

        if self.time_limit != -1:
            time_limit = datetime.timedelta(seconds=self.time_limit)
        else:
            time_limit = None
    
        start_time = time.time()
        #############################################################################################################################################
        print(f"Searching a distinguisher for {self.RD} rounds of SERPENT ...")
        self.cp_model = minizinc.Model()
        self.cp_model.add_file(self.mzn_file_name)
        self.cp_inst = minizinc.Instance(solver=self.cp_solver, model=self.cp_model)
        self.cp_inst["RU"] = self.RU
        self.cp_inst["RM"] = self.RM
        self.cp_inst["RL"] = self.RL
        self.cp_inst["RMU"] = self.RMU
        self.cp_inst["RML"] = self.RML
        self.cp_inst["NC"] = self.NC        
        self.cp_inst["offset"] = self.offset
        self.cp_inst["WU"] = self.WU
        self.cp_inst["WM"] = self.WM
        self.cp_inst["WL"] = self.WL
        self.result = self.cp_inst.solve(timeout=time_limit, 
                                         processes=self.num_of_threads, 
                                         verbose=False, 
                                         debug_output=Path("./debug_output.txt",
                                         intermediate_solutions=True),
                                         random_seed=randint(0, 100),
                                         optimisation_level=2)
        #############################################################################################################################################
        elapsed_time = time.time() - start_time
        print("Time used to find a distinguisher: {:0.02f} seconds".format(elapsed_time))
        print(f"Solver status: {self.result.status}")
        if minizinc.Status.has_solution(self.result.status) or self.result.status == minizinc.Status.ERROR:
            self.attack_summary, self.upper_trail, self.lower_trail = self.parse_solution()
            print(self.attack_summary)
            draw = DrawDL(self, output_file_name=self.output_file_name)
            draw.generate_distinguisher_shape()
        elif self.result.status == minizinc.Status.UNSATISFIABLE:
            print("Model is unsatisfiable") 
        elif self.result.status == minizinc.Status.UNKNOWN:
            print("Unknown error!")
        else:
            print("Solving process was interrupted")

    #############################################################################################################################################
    #############################################################################################################################################
    #  ____                           _    _             ____          _         _    _               
    # |  _ \  __ _  _ __  ___   ___  | |_ | |__    ___  / ___|   ___  | | _   _ | |_ (_)  ___   _ __  
    # | |_) |/ _` || '__|/ __| / _ \ | __|| '_ \  / _ \ \___ \  / _ \ | || | | || __|| | / _ \ | '_ \ 
    # |  __/| (_| || |   \__ \|  __/ | |_ | | | ||  __/  ___) || (_) || || |_| || |_ | || (_) || | | |
    # |_|    \__,_||_|   |___/ \___|  \__||_| |_| \___| |____/  \___/ |_| \__,_| \__||_| \___/ |_| |_|
    # Parse the solution and print the distinguisher's specifications

    def parse_solution(self):
        """
        Parse the solution and print the distinguisher's specifications
        """
        
        upper_trail = {"x": [[[0 for _ in range(self.NC)] for _ in range(4)] for _ in range(self.RU + self.RM + 1)],
                       "y": [[[0 for _ in range(self.NC)] for _ in range(4)] for _ in range(self.RU + self.RM)]}
        for r in range(self.RU):
            for row in range(4):
                upper_trail["x"][r][row] = self.result["xu"][r][row]
                upper_trail["y"][r][row] = self.result["yu"][r][row]
        for r in range(self.RU, self.RU + self.RM + 1):
            for row in range(4):
                upper_trail["x"][r][row] = self.result["xmu"][r - self.RU][row]
                if r < self.RU + self.RM:
                    upper_trail["y"][r][row] = self.result["ymu"][r - self.RU][row]
        lower_trail = {"x": [[[0 for _ in range(self.NC)] for _ in range(4)] for _ in range(self.RM + self.RL + 1)],
                       "y": [[[0 for _ in range(self.NC)] for _ in range(4)] for _ in range(self.RM + self.RL)]}
        for r in range(self.RM):
            for row in range(4):
                lower_trail["x"][r][row] = self.result["xml"][r][row]
                lower_trail["y"][r][row] = self.result["yml"][r][row]
        for r in range(self.RM, self.RM + self.RL + 1):
            for row in range(4):
                lower_trail["x"][r][row] = self.result["xl"][r - self.RM][row]
                if r < self.RM + self.RL:
                    lower_trail["y"][r][row] = self.result["yl"][r - self.RM][row]
        input_diff = ""
        for row in range(4):
            input_diff = hex(int("".join(list(map(str, self.result["xu"][0][row]))), 2))[2:].zfill(8) + input_diff
        input_diff_middle = ""
        for row in range(4):
            input_diff_middle =  hex(int("".join(list(map(str, self.result["xmu"][0][row]))), 2))[2:].zfill(8) + input_diff_middle
        output_mask_middle = ""
        for row in range(4):
            output_mask_middle = hex(int("".join(list(map(str, self.result["xml"][self.RM][row]))), 2))[2:].zfill(8) + output_mask_middle
        output_mask = ""
        for row in range(4):
            output_mask = hex(int("".join(list(map(str, self.result["xl"][self.RL][row]))), 2))[2:].zfill(8) + output_mask
        
        attack_summary = f"Attack summary:\n"
        attack_summary += f"Setting: RU: {self.RU}, RM: {self.RM}, RL: {self.RL}, RMU: {self.RMU}, RML: {self.RML}, offset: {self.offset}\n"
        attack_summary += f"Weights: WU: {self.WU}, WM: {self.WM}, WL: {self.WL}\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"input diff.: \n{input_diff}\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"input diff. middle: \n{input_diff_middle}\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"output mask middle: \n{output_mask_middle}\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"output mask: \n{output_mask}\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"PU   :  {self.result['PU']}\n"
        attack_summary += f"Q^2  :  {self.result['QL']}\n"
        attack_summary += f"Number of effective S-boxes in the middle:       {self.result['CMSbox']}\n"
        attack_summary += f"Number of effective bit-positions in the middle: {self.result['CMBit']}\n"
        attack_summary += "#"*50 + "\n"
        return attack_summary, upper_trail, lower_trail

    #############################################################################################################################################
    #############################################################################################################################################
    #############################################################################################################################################
    #   ____                                 _           ____  _              _               _                 _____   __   __              _   
    #  / ___| ___   _ __ ___   _ __   _   _ | |_  ___   / ___|| | _   _  ___ | |_  ___  _ __ (_) _ __    __ _  | ____| / _| / _|  ___   ___ | |_ 
    # | |    / _ \ | '_ ` _ \ | '_ \ | | | || __|/ _ \ | |    | || | | |/ __|| __|/ _ \| '__|| || '_ \  / _` | |  _|  | |_ | |_  / _ \ / __|| __|
    # | |___| (_) || | | | | || |_) || |_| || |_|  __/ | |___ | || |_| |\__ \| |_|  __/| |   | || | | || (_| | | |___ |  _||  _||  __/| (__ | |_ 
    #  \____|\___/ |_| |_| |_|| .__/  \__,_| \__|\___|  \____||_| \__,_||___/ \__|\___||_|   |_||_| |_| \__, | |_____||_|  |_|   \___| \___| \__|
    #                         |_|                                                                       |___/                                    
    # Take the clustering effect for the differential trail into account

    def compute_clustering_effect(self):
        """
        Take the clustering effect for the differential trail into account
        """
        
        print("#"*50)
        print("Computing the clustering effect for the upper trail ...")
        params_default = {"rounds" : self.RU,
                          "ncolumns" : self.NC,
                          "mode" : 2,
                          "sweight" : 0,
                          "endweight" : 4*self.NC,                            
                          "timelimit" : -1,
                          "fixedVariables" : {}}
        
        params = copy.deepcopy(params_default)
        for row in range(4):
            params["fixedVariables"][f"x_0_{row}"] = "".join(list(map(str, self.result["xu"][0][row])))
        for row in range(4):
            params["fixedVariables"][f"x_{self.RU}_{row}"] = "".join(list(map(str, self.result["xu"][self.RU][row])))
        UDiffEffect = Differential(params, exact=True)
        UDiffEffect.make_model()
        udiff_effect = UDiffEffect.solve(log=0)
        return udiff_effect

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#  _   _                    ___         _                __                   
# | | | | ___   ___  _ __  |_ _| _ __  | |_  ___  _ __  / _|  __ _   ___  ___ 
# | | | |/ __| / _ \| '__|  | | | '_ \ | __|/ _ \| '__|| |_  / _` | / __|/ _ \
# | |_| |\__ \|  __/| |     | | | | | || |_|  __/| |   |  _|| (_| || (__|  __/
#  \___/ |___/ \___||_|    |___||_| |_| \__|\___||_|   |_|   \__,_| \___|\___|
                                                                            
def loadparameters(args):
    '''
    Extract parameters from the argument list and input file
    '''

    # Load default values
    params = {"RU": 1,
              "RM": 2,
              "RL": 1,
              "RMU": 0,
              "RML": 0,
              "offset": 0,
              "WU": 1,
              "WM": 1,
              "WL": 1,
              "NC": 64,
              "np" : 8,
              "tl"  : -1,
              "solver"  : "ortools",
              "output"  : "output.tex"}

    # Override parameters if they are set on command line
    if args.RU is not None:
        params["RU"] = args.RU
    if args.RM is not None:
        params["RM"] = args.RM
    if args.RL is not None:
        params["RL"] = args.RL
    if args.RMU is not None:
        params["RMU"] = args.RMU
    if args.RML is not None:
        params["RML"] = args.RML
    if args.offset is not None:
        params["offset"] = args.offset
    if args.WU is not None:
        params["WU"] = args.WU
    if args.WM is not None:
        params["WM"] = args.WM
    if args.WL is not None:
        params["WL"] = args.WL
    if args.NC is not None:
        params["NC"] = args.NC
    if args.np is not None:
        params["np"] = args.np
    if args.timelimit is not None:
        params["timelimit"] = args.timelimit
    if args.solver is not None:
        params["solver"] = args.solver
    if args.output is not None:
        params["output"] = args.output

    return params

def main():
    '''
    Parse the arguments and start the request functionality with the provided
    parameters.
    '''
    
    parser = ArgumentParser(description="This tool finds a nearly optimum boomerang"
                                        "distinguisher for SKINNY family of block ciphers.",
                            formatter_class=RawTextHelpFormatter)
    
    parser.add_argument("-RU", type=int, default=0, help="Number of rounds for EU")
    parser.add_argument("-RM", type=int, default=3, help="Number of rounds in the middle")
    parser.add_argument("-RL", type=int, default=0, help="Number of rounds for EL")
    parser.add_argument("-RMU", type=int, default=0, help="Number of rounds passed probabilistically at the beginning of EM")
    parser.add_argument("-RML", type=int, default=0, help="Number of rounds passed probabilistically at the end of EM")
    parser.add_argument("-offset", type=int, default=0, help="Offset for the starting round")
    parser.add_argument("-WU", type=int, default=4, help="Weight factor for the upper trail")
    parser.add_argument("-WM", type=int, default=1, help="Weight factor for the middle part")
    parser.add_argument("-WL", type=int, default=4, help="Weight factor for the lower trail")


    parser.add_argument("-NC", type=int, default=32, help="Number of columns in each row, e.g., 64 for SERPENT-256")
    parser.add_argument("-np", type=int, default=8, help="Number of parallel threads")

    parser.add_argument("-tl", "--timelimit", type=int, default=36000, help="Time limit in seconds")
    parser.add_argument("-sl", "--solver", default="ortools", type=str,
                        choices=['gecode', 'chuffed', 'coin-bc', 'gurobi', 'picat', 'scip', 'choco', 'ortools', 'cplex', 'cbc'],
                        help="choose a cp solver\n")
    parser.add_argument("-o", "--output", default="output.tex", type=str, help="Output file name")

    # Parse command line arguments and construct parameter list
    args = parser.parse_args()
    params = loadparameters(args)
    dld = DiffLin(params)
    dld.search()
    # udiff_effect = dld.compute_clustering_effect()
    # print(f"Differential effect for upper trail: {udiff_effect}")

if __name__ == "__main__":
    main()
