# explicit

**Got a list of genes and you wonder what TFs regulate their expression? EXPLICIT is the right tool to try.**

The EXPLICIT approach has been developed to construct a gene expression predictor model for the plant species *Arabidopsis thaliana*. The predictor uses the expression of 1,678 transcription factor (TF) genes to predict the expression of 29,182 non-TF genes. It further enables downstream inference of TF regulators for genes and gene modules functioning in diverse plant pathways. Please check the [original paper](https://github.com/MaShisongLab/explicit#Reference) by Geng *et al.* for more details. The EXPLICIT package presented here enables users to 1. Infer TF regulators for their own gene modules; 2. Draw chord diagrams showing TF-target genes regulation for the modules; 3. Create custom gene expression predictors using their own gene expression data. (*Note: below is an example showing the analysis flow-chart for a gene module involved in vascular system development.*)

<a><img alt="The Analysis Work-flow" src="./data/working_flow.png" align="center" width="700" ></a>

## Table of Contents
- [Install](https://github.com/MaShisongLab/explicit#Install)
- [Usage](https://github.com/MaShisongLab/explicit#Usage)
   - Infer TF regulators for gene modules
   - Draw chord diagrams showing TF-target genes regulation for the modules
   - Create custom gene expression predictor
- [Reference](https://github.com/MaShisongLab/explicit#Reference)

## Install
This package requires [Perl](https://www.activestate.com/products/perl/downloads/), [R](https://www.r-project.org/), and the [circlize](https://www.rdocumentation.org/packages/circlize/) package in R. 

`circlize` can be installed within an R console via the command:

```R
install.packages("circlize")
```
[MATLAB](https://www.mathworks.com/products/matlab.html) is optional. Only required if you want to create your own predictor model using custom expression data. 

Once the required software is installed, just download or clone the whole package to a local computer and start using it from the package's home directory. 

## Usage

### 1. Infer TF regulators for gene modules

#### a. Prepare the module file
The file used to store gene modules information is `modules_to_analyze.txt`. It is preloaded with 1,085 gene modules identified from a GGM gene co-expression network described in the paper by [Geng et al.](https://github.com/MaShisongLab/explicit#Reference). The following analysis will proceed with these preloaded modules. On the other hand, you can also edit the file, replacing these modules with your own ones. The file has the following format, with the first column being gene ids and the second column being module names. The two columns are separated by a tab. For gene ids, only standard Arabidopsis AGI ids are currently supported. Multiple modules can be analyzed at the same time. <i>Once finish editing, save the file without changing its name</i>.
```
Gene_Name   ModuleID
AT1G25360   Module138
AT2G22340   Module138
AT5G75660   Moudle138
AT2G22130   Module139
AT4G12350   Module139
.........   ........
```
#### b. Conduct enrichment assay to identify TF regulators for the modules
The Perl script `getArabidopsisRegulatorTFs.pl` will do the job. It takes the modules from the file `modules_to_analyze.txt` to conduct an enrichment assay to identify potential TF regulators. Results are outputted to a file named `results.regulator.tfs.txt`, which can be open in EXCEL.

The command to use:
```shell
perl getArabidopsisRegulatorTFs.pl
```
Here is an example of the output results:
![alt-text](./data/results_sample.png "Output result example")


### 2. Draw chord diagrams showing TF-target genes regulation for the modules
#### a. Obtain the chord-list file for a module of interest
The Perl script `getChordLists.pl` will extract the TF-target gene pairs from the `results.regulator.tfs.txt` for the module specified. The results are outputted to a file named `chord.lists.txt`, which will be used in the next step to draw a chord diagram. 
```shell
# The command line format is :  perl getChordLists.pl XXXXXXX, where XXXXXX is the name of the module.
# For example, to extract chord-list file for Module0105
perl getChordLists.pl Module0105
```
For other module, just replace `Module0105` with the name of that module.

#### b. Draw the chord diagram according to the chord-list in R
The function `drawChordDiagram` in R will be used to draw the chord diagram. It has two inputs: `chordfile` specifyies the name of the chord-list file, `ratio` specifies the relative size of target gene area occupies. The function requires the `circlize` package. <br><br>
Open an R console and navigate to the home directory of the explicit package, which contains the `chord.lists.txt` file. Within the R console, type the following commands:
```R
source("Rscripts.R")   # Contains scripts for the drawChordDiagram function.
library("circlize")
drawChordDiagram(chordfile = "chord.lists.txt", ratio = 1)
drawChordDiagram(chordfile = "chord.lists.txt", ratio = 0.6)
```
Repeat step <b>a</b> and <b>b</b> to draw diagrams for another module.

#### c. Use a single command in R to draw the diagram
You can directly issue the following commands within an R console to draw Chord diagrams for modules:
```R
source("Rscripts.R")  
library("circlize")

# To draw diagram for Module0105
directChordDiagram( module="Module0105", ratio = 1, tfnum = 50, targetnum = 15)

# For other modules
directChordDiagram( module="Module0084", ratio = 1, tfnum = 50, targetnum = 15)
directChordDiagram( module="Module0081", ratio = 1, tfnum = 50, targetnum = 15)
directChordDiagram( module="Module0105", ratio = 1, tfnum = 50, targetnum = 15)
```
`tfnum` and `targetnum` specify, respectively, the maximum number of TFs and target genes to be included within the chord Diagram.

### 3. Create custom gene expression predictor 
Currently, we have only the gene expression predictor model for *Arabidopsis thaliana*. We are working on predictor models for other species. At the same time, you can also create your won custom gene expression predictor. However, a large number of training samples are required for training the model. The number should be at least 5 - 10 times larger than the number of input TFs.

The MATLAB function <B>explicit</B>, as specified within the file "explicit.m", will be used to create the model. The file can be found within the root folder of the package. MATLAB is required for this analysis.
```matlab
mdl = explicit( TF_expression, TG_expression, TF_name, TG_name)
```
`TF_expression`: the expression matrix for TF, with rows representing samples and columns representing genes <br>
`TG_expression`: the expression matrix for target genes, with rows representing samples and columns representing genes <br> 
`TF_name`: the names of the TF genes <br>
`TG_name`: the names of the target genes <br>

Here is a detailed procedure to create the Arabidopsis predictor.

#### a. Obtain the Arabidopsis gene expression matrix
Download two matrices `At.matrix.demo.h5` and `At.matrix.full.h5` from [Figshare](https://figshare.com/s/0c838ad4ef6a764daf53) (https://figshare.com/s/0c838ad4ef6a764daf53) , and place them within the <b>root directory of the EXPLICIT package </b>. `At.matrix.full.h5` is a full matrix with 24545 samples, while `At.matrix.demo.h5` has 5000 randomly selected samples from the full matrix. <i>We recommend to work with the smaller matrix `At.matrix.demon.h5` first</i>, as it requires less computational resources. Both are hdf5 format files with the following data structure: 
```bash
At.matrix.demo.h5
├─expression_log2cpm  		(5000 samples [row] X 38194 genes [column])
├─gene_name			(38194 genes)
├─rnaseq_id			(5000 samples)
├─idx_tf_gene			(specifying TF genes used for model construction)
├─idx_target_gene 		(specifying target genes used for model construction)
└─independent_samples_for_validation
   ├───expression_log2cpm  	(2 samples [row] X 39184 genes [column])
   ├───gene_name		(38194 genes)
   └───sample_id		(2 samples)
   
At.matrix.full.h5
├─expression_log2cpm  		(24545 samples [row] X 38194 genes [column])
├─gene_name			(38194 genes)
├─rnaseq_id			(24545 samples)
├─idx_tf_gene			(specifying TF genes used for model construction)
├─idx_target_gene 		(specifying target genes used for model construction)
└─independent_samples_for_validation
   ├───expression_log2cpm  (2 samples [row] X 38194 genes [column])
   ├───gene_name		(38194 genes)
   └───sample_id		(2 samples)
```
#### b. Create the expression predictor
Conduct the following analysis within a MATLAB console.
```matlab
%navigate to and start within the root directory of the explicit package
mtx_demo = h5read("At.matrix.demo.h5","/expression_log2cpm");
gene_name = h5read("At.matrix.demo.h5","/gene_name");

% which of the 38194 genes are TFs to be used
itf = h5read("At.matrix.demo.h5","/idx_tf_gene") == 1; 	

% target genes to be used
itarget = h5read("At.matrix.demo.h5","/idx_target_gene") == 1; 	

% obtain the TF expression matrix
tf_mtx_demo = mtx_demo(:,itf); 		

% obtain the target gene expression matrix
target_mtx_demo = mtx_demo(:,itarget); 				
tf_name = gene_name(itf);
target_name = gene_name(itarget);

% this produces the predictor model
mdl_demo = explicit( tf_mtx_demo, target_mtx_demo, tf_name, target_name);

% A look into the predictor model
mdl_demo 	

% The first 5 significant TF-target gene pairs
mdl_demo.SigEdges(1:5,:)  
```

The output is:
```bash
mdl_demo = 

  explicit with properties:

                          beta: [1679x29182 double] 	% beta coefficients (1 intercept + 1678 TFs) X 29182 genes
                   beta_pvalue: [1679x29182 double] 	% pValues for the beta coefficiets
                       TF_name: [1x1679 string]     	% names for 1 intercept + 1678 TFs
                   Target_name: [1x29182 string]    	% names for 29182 gene
                         NRMSE: 0.0660			% NRMSE for all the training samples
         Correlation_by_sample: [5000x1 double]		% R of predicted &acutal expression for 5000 samples
    Correlation_by_target_gene: [29182x1 double]    	% R of predicted & actuall exp of each gene across all samples
                           SST: [29182x1 double]    	% SST for regression model of every gene
                           SSR: [29182x1 double]    	% SSR for regression model of every gene
                           SSE: [29182x1 double]	% SSE for regression model of every gene
                         Fstat: [29182x1 double]	% Fstat for regression model of every gene
                       Fpvalue: [29182x1 double]	% Fpvalue for regression model of every gene
                      SigEdges: [252595x4 table]	% Significant TF-target gene pairs with pValue <= 0.00001

mdl_demo.SigEdges(1:5,:) =
       Gene            TF          beta      beta_pvalue
    ___________    ___________    _______    ___________

    "AT1G01020"    "AT1G01010"     0.1231     1.373e-23 
    "AT1G01020"    "AT1G15790"    -0.0596     1.074e-06 
    "AT1G01020"    "AT1G25580"    -0.0828     2.755e-06 
    "AT1G01020"    "AT1G30490"     0.0908     1.053e-06 
    "AT1G01020"    "AT2G01060"     0.1191     7.535e-06 
```
#### c. Use the predictor model to predict independent samples
We have generated two indepdent RNA-Seq datasets from Arabidopsis shoot and root samples. These two datasets were not used in the model training. Below we see how well the predictor model predict these two samples. <i>Note: the order of the gene names within the test samples are the same as those within demo matrix. </i> 
```matlab
test_mtx = h5read("At.matrix.demo.h5","/independent_samples_for_validation/expression_log2cpm");

% test_mtx contains two rows, for 'root' and 'shoot' samples, respetively.
test_sample_id = h5read("At.matrix.demo.h5","/independent_samples_for_validation/sample_id");
test_tf_mtx = test_mtx(:,itf);
actual_target_mtx = test_mtx(:,itarget);  

% [intercept (values of 1) + TF's expression matrix] X predictor's beta coefficient matrix to
% generate the predicted expression matrix for the target genes.
predicted_target_mtx = [ones(size(test_tf_mtx,1),1) test_tf_mtx] * mdl_demo.beta ;

% calculate the correlation between predicted and actual expression
corr( actual_target_mtx(1,:)', predicted_target_mtx(1,:)') % The correlation for root is 0.9920
corr( actual_target_mtx(2,:)', predicted_target_mtx(2,:)') % The coorelation for shoot is 0.9900

% Calculate NRMSE (Normalized Root Mean Square Error)
% NRMSE for root & shoot are 0.0858 & 0.1016
residual_mtx = predicted_target_mtx - actual_target_mtx ;
NRMSE = sqrt(sum(residual_mtx.^2, 2) ./ sum( actual_target_mtx.^2, 2)) 
```
<b>Note: If you want to use your own RNA-Seq datasets to test the model, make sure the gene names matched in the same order to the gene names of the TF and target gene matrices used above, and the gene expression value should be log2 transformed (log2(CPM + 1)). Alternatively, you can arrange the TFs and target genes of the demo matrix with the same order as your RNA-Seq datasets, and create the gene expression predictor accordingly. </b>

#### d. Investigate how the number of training samples affects the predictor power
The number of training samples affects the predictor's predicting power. The function <b>`explicit_eosn`</b>, standing for effect of sample number, investigates such effects. Its inputs are `(TF_expression, Target_expression, TestSampleNum)`, with `TestSampleNum` being the number of samples randomly selected  and hold out as test samples.
```matlab
% hold out 500 samples as test samples.
mdl_eosn = explicit_eosn( tf_mtx_demo, target_mtx_demo, 500) 
mdl_eosn
mdl_eosn.stat
```
The output is:
```bash
    Run    TrainingSampleNum    R_training    NRMSE_training    R_test     NRMSE_test
    ___    _________________    __________    ______________    _______    __________

     1           1700            0.99995        0.0071425       0.68378       0.931  
     2           1725            0.99988         0.010676       0.79395     0.63758  
     3           1750            0.99983         0.013001       0.82836     0.54658  
     4           1775            0.99975         0.015839       0.86977     0.43926  
     5           1800            0.99969         0.017487       0.88968     0.38962  
     6           1850            0.99958         0.020327       0.91018     0.34773  
     7           1900            0.99943         0.023749       0.92621     0.30313  
     8           1950             0.9993         0.026351       0.93624     0.27771  
     9           2000             0.9992         0.028108       0.94274     0.26096  
    10           2100            0.99897         0.031836       0.95118      0.2338  
    11           2200            0.99877         0.034816       0.95963     0.20954  
    12           2300            0.99854         0.037856       0.96267     0.20082  
    13           2400            0.99838         0.039831       0.96647     0.18906  
    14           2500            0.99819         0.042092       0.96992     0.17783  
    15           3000            0.99737          0.05082       0.97769     0.15087  
    16           3500            0.99679         0.056046         0.981     0.13833  
    17           4000            0.99629         0.060238       0.98307        0.13  
    18           4500            0.99587          0.06356       0.98441     0.12439  
```
Eighteen predictor models were built with between 1700 and 4500 training samples, and the predicting accuracy on test samples (R_test) increased along with the number of training samples.
#### e. Perform K-fold Cross-Validation
K-fold Cross-Validation can be also used to test the predictor's performance. The function <b>`explicit_kfcv`</b> does the job. Its inputs are `(TF_expression, Target_expression, tf_name, target_name, repeats, folds)`, with `repeats` and `folds` being the number of repeats and the folds for the analysis.
```matlab
mdl_kfcv = explicit_kfcv(tf_mtx_demo, target_mtx_demo, tf_name, target_name, 5, 10) % 5 repeats of 10-fold CV
mdl_kfcv
mdl_kfcv.CV_Stat
mdl_kfcv.AllEdges(1:5,:)
```
The output is:
```bash
mdl_kfcv = 

  explicit_kfcv with properties:

                 Total_repeats: 5
                   Fold_number: 10
    Correlation_of_target_gene: [29182x50 double]	% R of target genes in test samples in all 50 CV runs
              Target_gene_name: {29182x1 cell}
                      AllEdges: [48996578x9 table]	% The statistics for each coefficient across 50 CV runs
                       CV_Stat: [50x6 table]		% R and NRMSE for training and test samples in each CV run
                       
mdl_kfcv.CV_Stat =
   Repeat    Fold_No    NRMSE_training    NRMSE_test    R_training    R_test 
    ______    _______    ______________    __________    __________    _______

      1          1          0.063782        0.12305       0.99584      0.98481
      1          2          0.063719        0.12509       0.99585      0.98418
      1          3          0.063326        0.13161        0.9959      0.98232
      1          4          0.063505        0.13506       0.99588      0.98209
      1          5          0.063122        0.14025       0.99592      0.98067
      1          6          0.063771        0.11977       0.99584      0.98526
      1          7          0.063166        0.13044       0.99591      0.98294
      1          8          0.063427         0.1368       0.99589      0.98111
      1          9          0.063389        0.13608       0.99589      0.98168
      1         10           0.06349        0.13053       0.99588      0.98278
      2          1          0.063595        0.13268       0.99586      0.98234
      2          2          0.063495        0.12829       0.99587      0.98335
      2          3          0.063617        0.12701       0.99586      0.98347
      2          4           0.06349        0.12768       0.99588      0.98357
      2          5           0.06351        0.13197       0.99587      0.98255
     ...	...		...	       ...	     ...	  ...

mdl_kfcv.AllEdges(1:5,:) =
     Gene            TF          beta      beta_pvalue    CV_beta_mean    CV_beta_std    beta_bias    beta_relative_bias    relative_std
    ___________    ___________    _______    ___________    ____________    ___________    _________    __________________    ____________

    "AT1G01020"    "intercept"    0.58075       0.02262        0.58822        0.13633       0.00747           0.0129             0.2318   
    "AT1G01020"    "AT1G01010"    0.12307     1.373e-23        0.12226        0.01281      -0.00081           0.0066             0.1048   
    "AT1G01020"    "AT1G01030"    0.01859       0.08915        0.01981        0.00563       0.00122           0.0656             0.2842   
    "AT1G01020"    "AT1G01060"    -0.0292      0.001083       -0.02862        0.00599       0.00058           0.0199            -0.2093   
    "AT1G01020"    "AT1G01250"    0.02899       0.00919        0.02742        0.00629      -0.00157           0.0542             0.2294  

```
`mdl_kfcv.AllEdges` can be used to check if any of the coefficient is stable or not across the CV runs.
#### f. Perform Cross-Validation on independent samples
The function <b>`explicit_cv`</b> can be used to conduct Cross-Validation. Its inputs are `(Training_TF_expression, Training_traget_expression, Test_TF_expression, Test_target_expression, Test_sample_ids)`. The function uses Training TF and target expression matrices to build the predictor model, and then test the model on the Test samples. It provides a convinient way to test different models, which are specified by different input training samples.
```matlab
test_target_mtx = test_mtx(:,itarget);

% Use all 5000 samples in the demo matrix to build the preidctor model, 
% and test on the shoot and leaf RNA-Seq samples.
mdl_cv_5000 = explicit_cv( tf_mtx_demo, target_mtx_demo, test_tf_mtx, test_target_mtx, test_sample_id)

% Use only 4000, 3000, or 2000 samples to build the predictor model
mdl_cv_4000 = explicit_cv( tf_mtx_demo(1:4000,:), target_mtx_demo(1:4000,:), test_tf_mtx, test_target_mtx, test_sample_id)
mdl_cv_3000 = explicit_cv( tf_mtx_demo(1:3000,:), target_mtx_demo(1:3000,:), test_tf_mtx, test_target_mtx, test_sample_id)
mdl_cv_2000 = explicit_cv( tf_mtx_demo(1:2000,:), target_mtx_demo(1:2000,:), test_tf_mtx, test_target_mtx, test_sample_id)
mdl_cv_5000
mdl_cv_5000.Test_Sample_Stat
```
The output is:
```bash
mdl_cv_5000 = 

  explicit_cv with properties:

             Training_Sample_Num: 5000
               Actual_Target_Exp: [2x29182 double]
            Predicted_Target_Exp: [2x29182 double]
                 Test_Sample_Num: 2
                Test_Sample_Stat: [2x3 table]
           Test_Sample_NRMSE_all: 0.0933
    Test_Sample_Mean_Correlation: 0.9910

mdl_cv_5000.Test_Sample_Stat =
    Sample     R_test     NRMSE_test
    _______    _______    __________

    'root'     0.99196     0.085795 
    'shoot'    0.99002      0.10157 
```
When using 5000, 4000, 3000, or 2000 samples to train the predictor model, the mean correlation for test samples are 0.9910, 0.9891, 0.9841, 0.9572, respetively.
#### g. Use the full matrix to perform the analysis
Next we will analyze the full matrix. <i>Note: Since the matrix is large, it requires a large amount of computational resource. It is recommended to have at least 80G memory available.</i>
```matlab
% clear all previous variables to free memory usage.
clearvars ; 


% Step A. Create the predictor
mtx_full = h5read("At.matrix.full.h5","/expression_log2cpm");
gene_name = h5read("At.matrix.full.h5","/gene_name");

% which of the 31984 genes are TFs to be used
itf = h5read("At.matrix.full.h5","/idx_tf_gene") == 1;   

% target genes to be used
itarget = h5read("At.matrix.full.h5","/idx_target_gene") == 1;  

% obtain the TF expression matrix
tf_mtx_full = mtx_full(:,itf);    

% obtain the target gene expression matrix
target_mtx_full = mtx_full(:,itarget);                      
tf_name = gene_name(itf);
target_name = gene_name(itarget);

% this produces the predictor model.
mdl_full = explicit( tf_mtx_full, target_mtx_full, tf_name, target_name); 
mdl_full    

% The predictor has 3298936  SigEdges (TF-target gene pairs) with pValue <= 0.00001
% Next the SigEdges with pValue <= 1e-9 are extracted and saved to a file named "Arabidopsis.SigEdges.1e-9.txt". 
% This file is the same as the file "At.SigEdges.txt" within the data directory.
% There are 980736 SigEdges with pValue <= 1e-9
i = mdl_full.SigEdges{:,4} <= 1e-9 ;
sum(i)       
writetable( mdl_full.SigEdges(i,:), "Arabidopsis.SigEdges.1e-9.txt", "Delimiter","tab")


% Step B. Investigate how the number of training samples affects the predictor power
mdl_eosn = explicit_eosn( tf_mtx_full, target_mtx_full, 3000)
mdl_eosn
mdl_eosn.stat

% Step C. K-fold Cross-Validation
mdl_kfcv = explicit_kfcv(tf_mtx_full, target_mtx_full, tf_name, target_name, 5, 10)
mdl_kfcv
mdl_kfcv.CV_Stat
mdl_kfcv.AllEdges(1:5,:)

% Step D. Cross-Validation on independent samples
test_mtx = h5read("At.matrix.full.h5","/independent_samples_for_validation/expression_log2cpm");
test_sample_id = h5read("At.matrix.full.h5","/independent_samples_for_validation/sample_id");
test_tf_mtx = test_mtx(:,itf);
test_target_mtx = test_mtx(:,itarget);
mdl_cv_full = explicit_cv( tf_mtx_full, target_mtx_full, test_tf_mtx, test_target_mtx, test_sample_id)
mdl_cv_full
mdl_cv_full.Test_Sample_Stat
```
The output is:
```bash
mdl_cv_full = 

  explicit_cv with properties:

             Training_Sample_Num: 24545
               Actual_Target_Exp: [2x29182 double]
            Predicted_Target_Exp: [2x29182 double]
                 Test_Sample_Num: 2
                Test_Sample_Stat: [2x3 table]
           Test_Sample_NRMSE_all: 0.0796
    Test_Sample_Mean_Correlation: 0.9934

mdl_cv_full.Test_Sample_Stat =
    Sample     R_test     NRMSE_test
    _______    _______    __________

    'root'     0.99449     0.070692 
    'shoot'    0.99224     0.089135 
```
The mean correlation for the test samples is 0.9934 for the full model, which is better than the model constructed with 5000 samples within the demo matrix. 

## Reference

Will update soon.
