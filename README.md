# The effect of combined MEK and EGFRinhibition on gene expression incancer patients
The constitutive activation the MAPK pathway by KRAS-mutant tumors has complicated the de-
velopment of an effective treatment. To mitigate this resistance mechanism, simultaneous inhibition
of MEK, EGFR, and HER2 was studied in three phase I clinical trials in colorectal cancer (CRC),
non-small cell lung cancer (NSCLC), and pancreatic cancer patients. The results of these trials
were modest compared to the response in cancer cell lines and pre-clinical models. To elucidate
the differences in response, we obtained a gene signature for each drug used in the trials to deter-
mine whether each drug was inhibiting its target in the clinical setting. The gene signature was
obtained from the L1000 dataset and, by inverse transform sampling, we merged two approaches
to recover an accurate expression level of the genes: a Bayesian approach and the characteristic
direction method. To obtain a gene signature independent of proliferation we implemented a novel
approach by using a KS test and the gene expression data of drugs inhibiting targets independent
of the MAPK pathway. The performance of the pipeline was evaluated in eight independent cancer
cell line datasets and deemed valid in most. Finally, the gene signature was applied to the clinical
gene expression data, 85% and 95% patients showed a MEK and panHER inhibition signal, respec-
tively. Hence, the limited clinical efficacy must be attributed to other causes, such as a persistent
subpopulation of cells, a lack of apoptosis stimulation or immune response activation.
