class SingleCellGGM:
    

    gene_num = [];
    gene_name = [];
    pcor_all = [];
    pcor_sampling_num = [];
    coexpressed_cell_num = [];
    samples_num = [];
    RoundNumber = [];
    SigEdges = [];
    
    def __init__(self, x, round_num , gene_name, dataset_name = 'na'):
        
        import numpy as np
        import pandas as pd
        import random
        import time

        seed = 98
        cut_off_pcor = 0.03
        cut_off_coex_cell = 10
        selected_num = 2000
        data_name = dataset_name
        
        nrow = x.shape[0]
        ncol = x.shape[1]
        
        aa = x > 0
        aa = aa.astype(float)
        cellnum = sum(aa)
        cellnum = cellnum.A.astype(int)
        cellnum = cellnum.flatten()

        coex = np.transpose(aa) * aa
        coex = coex.A
        
        cov_all = np.cov(x,rowvar = False)
        
        rho = np.corrcoef(x, rowvar = False)
        del x
        
        gene_id = gene_name
        
        pcor_all = np.ones((ncol,ncol))
        pcor_sampling_num = np.zeros((ncol,ncol))
        time_trend = np.zeros((100,1))
        
        print('Calculating pcor in %d iterations.\n' %round_num)
        random.seed(seed)
        for i in range(1,round_num+1) :
            loop_start_t = time.time()
            j = random.sample(range(0,ncol),selected_num)
            j = np.asarray(j)
            cov_x = cov_all[j,].take(j,1)
            ix = np.linalg.inv(np.mat(cov_x))
            d = np.diag(np.sqrt(np.diag(ix)))
            d = np.linalg.inv(d)
            pc = - d * ix * d
            pc = np.eye( selected_num ) * 2 + pc

            def loopjit(j,selected_num,pcor_sampling_num,pc,pcor_all):
                for m in range(0,selected_num) :
                            for n in range(0,selected_num) :
                                r = j[m]
                                s = j[n]
                                if r > s :
                                    pcor_sampling_num[r,s] = pcor_sampling_num[r,s] + 1
                                    if abs(pc[m,n]) < abs(pcor_all[r,s]) :
                                        pcor_all[r,s] = pc[m,n]
            loopjit(j,selected_num,pcor_sampling_num,pc,pcor_all)
            loop_time = time.time() - loop_start_t
            idx_time  = np.mod(i,100) - 1
            time_trend[idx_time] = loop_time
            average_loop_time = np.mean(time_trend)
            time_left = (round_num - i) * average_loop_time / 3600

            if i == 100 :
                print('Estimated to complete in %.2f hours.\n' %time_left)
            if np.mod(i,1000) == 0 & i < round_num :
                print('%d iterations done. %.2f sec/iteration in average. %.2f hours to go.\n'  %(i, average_loop_time, time_left))
            if i == round_num :
                print('%d iterations done.\n' %i)
                
        
        self.samples_num = nrow
        self.gene_num = ncol
        self.RoundNumber = round_num
        self.gene_name = gene_id
        self.pcor_sampling_num = pcor_sampling_num.astype(int)
        
        pcor_all[pcor_sampling_num == 0] = 0
        self.pcor_all = pcor_all.astype(float)
        self.coexpressed_cell_num = coex.astype(int)

        idx = np.where((pcor_all >= cut_off_pcor)&(pcor_all < 1)&(coex >= cut_off_coex_cell))
        e1 = list(gene_id[idx[0]])
        e1n = list(cellnum[idx[0]].astype(int))
        e2 = list(gene_id[idx[1]])
        e2n = list(cellnum[idx[1]].astype(int))
        e3 = list(pcor_all[idx])
        e4 = list(pcor_sampling_num[idx].astype(int))
        e5 = list(rho[idx])
        e6 = list(coex[idx].astype(int))
        e7 = [data_name]*len(e1)
        self.SigEdges = {'GeneA':e1,'GeneB':e2,'Pcor':e3,'SamplingTime':e4,\
                        'r':e5,'Cell_num_A':e1n,'Cell_num_B':e2n,\
                        'Cell_num_coexpressed':e6,'Dataset':e7}
        self.SigEdges = pd.DataFrame(self.SigEdges)
