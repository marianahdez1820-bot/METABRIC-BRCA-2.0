# Transcriptomic Signatures Predicting Survival and Recurrence in Estrogen Receptor-Positive Breast Cancer: A Machine Learning Analysis of the METABRIC Cohort

This repository contains the pipeline for preprocessing, feature selection, hyperparameter tuning, training, and external validation, as well as misclassification and enrichment analyses for ER-positive breast cancer data.

## Execution Order & Pipeline Structure

The folders and files are ordered in the exact sequence they should be executed, as objects created in earlier scripts are required for subsequent steps.

---

### 1. METABRIC Preprocessing

#### `1. METABRIC preprocessing/1.1 Data preproccesing.R`
This script handles the initial data pipeline:
* **Prerequisite:** Having the downloaded data corresponding to METABRIC metadata and counts.
* **Steps:** Data loading $\rightarrow$ Duplicate gene management $\rightarrow$ Metadata preprocessing.
* **Key Outputs:** Generates the core objects: `metadata.ER_POS_SURV`, `metadata.ER_POS_REC`, and `counts_data`.

---

### 2. Feature Selection

#### `2. Feature selection/2.1 Boruta.R`
Runs the Boruta feature selection algorithm for signature creation

> **Important Pipeline Notes:**
> * **Analysis Mode:** The output depends on the object loaded in **Line 9**. Use `metadata.ER_POS_SURV` for survival analysis, or switch to `metadata.ER_POS_REC` for recurrence.
> * **Gene Selection:** **Line 43 (Section 5.2)** defines the threshold for how many genes will be analyzed.
> * **Parallelization:** Adjust **Line 67 (Section 6.2)** to allocate the number of CPU threads/cores.
> * ** Number of runs:** **Line 90 (Section 7)** to modify the number of runs
> * **Tentative Genes:** **Line 117** forces a final decision on tentative genes before saving results.

---

### 3. Signature Preparation

> **Section Workflow:** This section outputs the final data objects and optimal hyperparameters used for downstream regression. **Run either 3.1 OR 3.2** depending on your desired analysis, followed by **3.3**.

#### `3. Prepare_linreg_ER/3.1 Prepare_linreg_ER_rec_global.R` (Recurrence Path)
* **Purpose:** Prepares data specifically for **recurrence analysis**.
* **Function:** Imports the gene signature from Boruta step 2.1 and utilizes `metadata.ER_POS_REC` to create `proof_genes_pt` and `ml_metadata`.

#### `3. Prepare_linreg_ER/3.2 Prepare_linreg_ER_global_surv.R` (Survival Path)
* **Purpose:** Prepares data specifically for **survival analysis**.
* **Function:** Imports the gene signature from Boruta step 2.1 and utilizes `metadata.ER_POS_SURV` to create `proof_genes_pt` and `ml_metadata`.

> **Object notes:**
> * **Signature Selection:** For both scripts, **Line 16 (Section 1.2)** is where you insert your desired signature list object.
> * **`proof_genes_pt`:** This object contains patients in rows and signature genes in columns, alongside outcome columns (*EVENT_STAT*, *EVENT_MON*, and *surv_obj*). Because it does *not* contain an intrinsic label stating whether it represents survival or recurrence, be highly careful about which preparation file you ran.
> * **`proof_genes`:** Another object with similar name and almost similar impoirtance. This object is used in many downtream analysis and forms part of the validation data preprocessing to filter only the genes that are usefull for use in the signature and to check if the validation dataset contains all the necessary genes.
> * **`ml_metadata`:** An outcome-neutral renaming of the filtered patient metadata, allowing downstream scripts to run without manual object name changes.
> * ** External Validation Loop (GSE2034 / GSE96058):** As noted in the paper, some signature genes may be missing in external datasets. During data preparation of either of these datasets, pass the generated `common_genes_meta.gse2034` or `common_genes_meta.gse96058` objects back into **Line 16** of either script 3.1 or 3.2 to re-run the pipeline with a fully compatible dataset-wide signature.
> * **Last preprocessing step:** At this point is where the data are scaled and centered for the patients and genes that will be going into the training and initial testing.

#### `3. Prepare_linreg_ER/3.3 Hyperparameter_selection.R`
* **Purpose:** Selects hyperparameters for the Cox regression model using `proof_genes_pt`.
* **Output:** Generates `best_params` (a tibble containing optimal parameters). 
* **Tip:** To inspect alternative hyperparameter options instead of automatically choosing the top C-index, change `select_best()` to `show_best()` on **Line 85 (Section 3.2)**.

---

### 4. Cox Regression

> **Section Workflow:** Subdivided into the main model (`4.1`) and downstream analyses (`4.2`). This section automatically adapts to whichever file you ran back in Section 3, so ensure your active global objects match your intended analysis (Survival vs. Recurrence).

#### `4. Cox Regression/4.1 Cox_regression_global.R`
* **Purpose:** Handles the METABRIC dataset split, model training, and initial testing.
* **Required Inputs:** `proof_genes_pt`, `ml_metadata`, and `best_params` (optional; if missing, the script recalculates the tuning grid, though this diverges from the published workflow).
* **Key Results & Outputs:** 
  * Main paper findings derive from `summary_cox` and `independent_prog`.
  * Generates downstream analysis objects `proof_genes_pt.cox` (test data with metadata),  `final_fit` (trained model), `true_cut` (cutpoint obtained from *surv_cutpoint()*) and image generation inputs (`fit_km`, `plot_roc`, `facet_labels`, `cox_p_metabric`) utilized later in folder `7. Images`.

#### `4.2 Further_analysis`
> **Section Workflow:** This section involves 3 scripts focused on identifying potential biases and causes for misclassification

##### `4. Cox Regression/4.2 Further_analysis/4.2.1 Bias analysis.R`
* **Purpose:** Identyfies multiple metrics that evaluate performance, assumptions and posible confounding factors or interactions with clinical characteristics
* **Performance evaluation:** Calculates the C-index and AUC with acompanying confidence intervals and empirical p value after 100 refolds of Monte Carlo cross validation **(Section 5)**, also calculates Brier score **(Section 2.1)**
* **Asumption evaluation:** Calculates Schonfeld and Martingales residuals **(Section 2.2)**. It also analyses the contribution of each gene to the score **(Section 3)**
* **Posible confounders:** Evaluates Wilcoxon test between patient status across distinc molecular subtypes **(Section 6.4)** and treatment modalities **(Section 6.5)**. Finally it evaluates a Cox model of interaction between score and treatment **(Section 6.5.2)**.
> **Notes:**
> If the previous files were run correctly this script shouldnt have the need to be modified except on the labels of the plots to change from survival tu recurrence and viceversa.

##### `4.2.2 Misclassification_analysis.R`
* **Purpose:** Applies the misclassification score to patients, plots the highest scored patients, obtains clinical characteristics of these patients and finally runs 2 out of 3 important steps to identify potential causal elements that led to the misclassification (multiple characteristic analysis via chi squared, fishers test and kruskal wallis; and logistic regression penalized by Firth method).
* **Inputs:** If all scripts have been run correctly the script should run smoothly since the input objects are `train_data` and `final_fit` (created on `4.1 Cox_regression_global.R`).
*  **Interations:** Only at 2 points could changes be made, those being at  line 25 **(Section 1.2.2)** to stop the plots from showing each time that for loop is performed and at line 89 **(Section 1.6.1)** to change the standard deviation that functions as threshold of selecting misclassified patients.
*  **Output objects:** The most important output object is `misclassification_diff_id ` which consist on those patients who during the for loop at the 3 time points where considered a misclassified patient at least once. This object is important for the differential expresion and GSEA analysis.

##### `4.2.3 Random_signatures.R`
* **Purpose:** Compare the gene set obtained from Boruta against 1000 random gene sets of equal size
* **Inputs:** If all scripts have been run correctly the script should run smoothly since the main entrance objects are `ml_metadata`, `counts_data` and `proof_genes`.
---

### 5. Validation

> **Section Workflow:** This is one of the more complex sections. It is suvbdivided into 3 sections corresponding to the 3 external validation sets (**GSE2034** `5.1 GSE2034`, **GSE96058** 5.2 GSE96058, **TCGA** `5.3 TCGA`). Each section is divided into 2 subsections which correspond to the data preparation and the model testing.

> * The workflow depends on survival or recurrence analysis
>   * On survival run `5.2 GSE96058` and `5.3 TCGA`. The latter with specifications given in its assigned section.
>   * On recurrence run `5.1 GSE2034` and `5.3 TCGA`. Once againthe latter with specifications given in its assigned section.

**IMPORTANT FOR ALL 3 CASES IN DATA PREPARATION:** The only external object is `proof_genes` created on folder `3. Prepare_linreg_ER/3.1 Prepare_linreg_ER_rec_global.R`. If all genes on `proof_genes` can be found on the processed database, **Section 4.3 (it is the same section in all 3 validation sets)** will output "All genes in the signature are on (Name of cohort)", if not, the same section will output how many genes are missing as well as which ones are those missing genes. In both cases an object with all the common genes between the `proof_genes` object and the validation data base is created (independently of it shares all or only a few) called `common_genes_meta.` with the dot followed by the name of the database in undercase (example `common_genes_meta.tcga`). If the output is the error stating that there are missing genes the workflow consists on returning to the utilized `3. Signature Preparation` subfolder and inserting the `common_genes_meta.` into **Line 16 (Section 1.2)** and reruning the workflow from there. If instead it outpus the message stating all genes are found, one can continue with the workflow normally.

**Output on all 3 cases of data preparation:** All preparation files output an object called `proof_genes_pt.` with the dot followed by the database name in underscore (`proof_genes_pt.tcga`, `proof_genes_pt.gse2034` and `proof_genes_pt.gse96058`). This object same as with `proof_genes_pt` from `3. Signature Preparation` consists on the patients on rows and scaled and centered genes on columns as well ad the outcome variables which are ignored to the model because of tidymodels recipe specification given at `4. Cox Regression/4.1 Cox_regression_global.R`. It also outputs metadata in differing named objects specified in each subsection

#### `5.1 GSE2034`
##### `5.1.1 GSE2034_prep`
* **Purpose:** Download, untar and read cel files followed by metadata preprocessing and count data preprocessing and finally prepares the object used for testing the signature. GSE2034 was only used for recurrence pipeline
* **Data procesing:** RMA normalization, anotation, duplicate gene management. Creation of `proof_genes_pt.gse2034` and scaling and centering of its gene counts.
* **Other output:** The other importatn output used in subsequent files is `metadata.gse.2034_er_pos`
