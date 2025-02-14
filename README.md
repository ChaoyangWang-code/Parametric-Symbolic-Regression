# A Parametric Symbolic Regression Framework for Discovering Unified Crack Growth Models from Diverse Experimental Data  

**Code License:** GPL 3.0  
**Data License:** CC BY 4.0  

## Overview  
We introduce a data-driven approach for constructing crack growth models, termed the Parametric Symbolic Regression (PSR) framework. The PSR framework directly leverages extensive experimental data encompassing multiple influencing factors to derive a unified crack growth model that is simple, well-generalized, and stable. Our PSR framework was implemented through extensive development and modification of GPTIPS 2.0 ([GitHub link](https://github.com/domsearson/gptips-2-0)).  

---

## 💡Key Features  
1. **Parametric Modeling Framework:** The PSR framework constructs a symbolic regression framework in which all candidate models are parameterized rather than specific.  
2. **Comprehensive Multi-Criteria Evaluation:** We have developed multi-criteria evaluation metrics to assess candidate models. Additionally, we introduced a multi-objective genetic algorithm to optimize these models.  
3. **Parallel Computing Acceleration:** This program includes parallel computing functionality, which facilitates the use of multi-core CPUs on both Windows and Linux (especially on high-performance computing platforms) to accelerate the symbolic regression process.  

---

## 🚀Getting Started  

### Software Information  
- This program is designed to run on MATLAB versions 2018a and above in both Windows and Linux environments.  
  - Note: Using versions earlier than 2018a may result in errors.  
- Existing results in the `Configuration and Results` folder were generated under: Linux, MATLAB 2021a.  

### Basic Usage Guide  
1. **Select the Configuration File**  
   - Copy `gpdemo_crack_growth_config.m` from the `Configuration and Results` folder to the root directory.  
   - Existing configuration files in this folder can be executed directly or modified as needed.  

2. **Running on Windows Platform**  
   - Set `gp.fitness.savepath` in `gpdemo_crack_growth_config.m`.  
   - Execute `gpdemo_crack_growth.m` to run the PSR code.  

3. **Running on Linux Platform (especially for HPC)**  
   - Activate the code for Linux platform computation in `gpdemo_crack_growth.m`.  
   - Set `gp.fitness.savepath` in `gpdemo_crack_growth_config.m`.  
   - Modify line 165 in `evalfitness_par.m` from:  
     ```matlab
     savePath = [gp.fitness.savepath, str, '\gen', str, '.mat'];
     ```  
     to:  
     ```matlab
     savePath = [gp.fitness.savepath, str, '/gen', str, '.mat'];
     ```  
   - Modify `run.sh` and `run1.sh` as per HPC platform requirements and submit jobs.  

---

## ✨Results Display (Windows Only)  
1. After the symbolic regression process  completes, locate the result folder and execute `Find_pareto_in_all_gen.m` in the `gp_show` folder.  
2. Pre-written visualization codes are stored in the `figure_of_results` subfolder. Run them directly or modify as needed.  

---

## 🗝️Extended Features  
The proposed PSR framework can not only be used to find crack propagation models, but also applies to various other fields that require unified models. It helps extract unified models from diverse experimental data.  

### Switch Dataset  
- Multiple datasets are stored in the `data` folder. Modify relevant code in `gpdemo_crack_growth_config.m` to call them.  
- Use code in the `data_make` folder to generate datasets.  After that, run `data_operation.m` to process the data for program compatibility.  

### Switch Candidate Model Space  
- Add new function/terminal nodes to the `gp_operators` folder. Select them in `gpdemo_crack_growth_config.m` to expand search space for addressing problems in different domains.  
**Note:** The fitting algorithm for model parameters may need to be specifically adjusted according to the characteristics of the candidate model. Adjust parameter fitting algorithms in `loss_cal_optimize.m` (located in the `find_parameters` folder).  

### Switch Multi-Criteria Evaluation Metrics  
- Users can modify the evaluation metrics for candidate models based on the target of different research domains. Modify evaluation metrics in `cal_loss_multidata.m` (within `find_parameters` folder).  
- Update corresponding content in `updatestats.m`.  
