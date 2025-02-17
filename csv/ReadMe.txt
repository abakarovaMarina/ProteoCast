ProteoCast data 

File Descriptions

1. Prediction Summary Files
	•	Dmel6.44PredictionsRecap.csv: summary of mutation effect predictions for all proteoforms, including FlyBase IDs, sequence length, multiple sequence alignment (MSA) properties, classification thresholds, structural data, and mutation counts.
	•	ConfidenceLocal_pLDDT.csv: provides local confidence scores along with AlphaFold2’s pLDDT scores. 

2. Polymorphism Data (DEST and DGRP Datasets)
	•	DEST_allRepresentativeProteoforms.csv / DGRP_allRepresentativeProteoforms.csv: polymorphism data mapped to representative proteoforms, including global and local confidence scores, mutational impact scores, ranks, and classification.
	•	DEST_maxFobsProteoforms.csv / DGRP_maxFobsProteoforms.csv: same as above, but filtered for proteoforms with maximum MSA diversity (F_obs).
	•	DGRP_SNPfreq.csv: SNP frequency data from the DGRP dataset, including reference and alternate alleles, allele counts, and minor allele frequencies (MAF).

3. Lethal and Hypomorphic Mutation Data
	•	Lethal_allRepresentativeProteoforms.csv / Lethal_maxFobsProteoforms.csv: lethal mutations mapped to representative proteoforms,  including global and local confidence scores, mutational impact scores, ranks, and classification.
	•	Hypomorphic_allRepresentativeProteoforms.csv / Hypomorphic_maxFobsProteoforms.csv: same as above, but for hypomorphic mutations.
	•	ProteoCast_Benchmark_Lethal_DEST2_DGRP.csv: contains the complete benchmark dataset.

4. Post-Translational Modifications (PTMs) and Short Linear Motifs (SLiMs) Data
	•	PTM.csv: list of post-translational modifications, including position, modified amino acid, and modification type.
	•	PTM_proteomeSensitivity1_4.csv: PTM annotations with sensitivity scores (local and global), structural information, and confidence scores.
	•	ELM_dataSensitivity.csv: list of short linear motifs (SLiMs) with sensitivity scores (local and global), structural information, and confidence scores.

5. Segmentation Data
	•	Segmentation_ALL_a1_2.csv to Segmentation_ALL_a2.csv: segmentation of GEMME signal into sensitivity states (0, 1,  2).

6. Protein Interaction 
	•	Interactions_GenPhys.csv: number of genetic and physical interactions for each FlyBase gene ID.

7. RNAi Experimental Data
	•	RNAi_experiments.csv: experimental details for RNA interference (RNAi) studies.
