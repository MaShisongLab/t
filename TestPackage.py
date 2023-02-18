class TestPackage:
    

    gene_num = [];
    gene_name = [];
    pcor_all = [];
    pcor_sampling_num = [];
    coexpressed_cell_num = [];
    samples_num = [];
    RoundNumber = [];
    SigEdges = [];
    
    def __init__(self, round_num , dataset_name = 'na'):
        
        import numpy as np
        import pandas as pd
        import random
        import time
        
        print('Calculating pcor in %d iterations.' %round_num)
       
