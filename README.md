<h1 align="center">Psychometric Network Inference</h1>

Network-based approaches have provided an important contribution for the understanding of the study of mental disorders. A growing number of statistical models, developed in the context of continuous variables in high-dimensional settings, are currently being used to infer dependencies between network elements (e.g., symptoms or behavioral elements) in psychometrics. However, psychometric datasets typically correspond to **low-dimensional statistical settings**, namely with a low number of variables collected from a large enough sample size and the variables collected are **ordinal** rather than Gaussian. 

In a large-scale simulation study [1], we tested and compared the performance of 14 methodological approaches including some still limitedly used in the context of psychometrics network inference. We assessed  the impact of various factors such as the sample size, the number of variables (i.e., network elements), the density of the true underlying graph and the number of ordinal levels. 

This page contains the scripts (`Simulations_ScenarioXX_script.R`) to perform the experiments from [1]. It also contains the complete set of output figures (`ScenXX_Figs.tgz`) for each of the 4 scenarios tested. 
The original dataset utilized for the simulation in Scenario 2 are accessible through a formal request to the National Institute of Mental Health (NIMH) at https://nda.nih.gov/.

## Requirements 
To run those scripts and perform our experiments, you will need R and the packages `optparse`, `igraph`, `BDgraph`, `qgraph`, `bootnet`, `GGMnonreg`, `GeneNet`,   `PAsso`,   `psychonetrics`,   `BGGM`, `e1071`.   

## Methods included in the comparison
Table with the 14 methods compared: 

| Method |  R Package | Function | Input | Parameters | 
| ----- | ----- | ----- | ----- | ----- |  
|  poly.mle |  qgraph   |  cor\_auto   | raw data| -  | 
|  poly.wls|  psychonetrics    |  ggm   \%$>$\%  prune   \%$>$\%  modelsearch   | raw data | - |
|  pears |  stats   |  cor   | raw data | - |
|  Glasso.poly |  qgraph  | EBICglasso  |  poly.mle   output | default| 
|  Glasso.pears |  qgraph  | EBICglasso  | pearson   output| default|  
|  ggmMS.poly |  qgraph  | ggmModSelect  | poly.mle   output| default| 
|  ggmMS.pears |  qgraph  | ggmModSelect  | pearson   output| default| 
|  GGMnr.neighsel |  GGMnonreg  | ggm\_inference  | raw data | boot=FALSE| 
|  GGMnr.boot.poly |  GGMnonreg  | ggm\_inference  | raw data | default| 
|  GGMnr.boot.pears |  GGMnonreg  | ggm\_inference  | raw data| method=  "polychoric"|  
|  BGGM.explore |  BGGM  | explore  | raw data |type="ordinal", impute = FALSE| 
|  BGGM.estimate |  BGGM  | estimate  | raw data |type="ordinal", impute = FALSE| 
|  ggmSS |  GeneNet | ggm.estimate    .pcor  | raw data| default| 
|  PAsso |  PAsso | PAsso  | raw data |default| 

## References: 
1. Claudia Delli Colli, Blerina Sinaimeri, Giuseppe F. Italiano, Igor Branchi and Catherine Matias (2024). Psychometric network inference: A comparative analysis. Submitted. [PsyArXiv preprint](https://osf.io/preprints/psyarxiv/mcj9a)
