# DNA Variation Analyzer

A Python-based bioinformatics project for identifying and summarizing DNA sequence variations by comparing sample sequences with a reference sequence.

## Project Overview

This project analyzes DNA sequences from an ATM gene dataset and compares each sample against the reference sequence **MK284930**.

The program identifies nucleotide differences, summarizes mutation patterns, and generates Excel reports and visualizations.

## Objectives

* Compare DNA sequences with a reference sequence
* Identify nucleotide substitutions
* Calculate mutations for individual samples
* Summarize mutations by position and mutation type
* Generate Excel reports for further analysis
* Create visualizations of mutation distributions

## Dataset

* **Gene:** ATM
* **Reference ID:** MK284930
* **Reference length:** 384 bp
* **Input file:** `ATM.fasta`
* **Number of records:** 47

Samples with sequence lengths different from the reference are recorded separately and excluded from direct comparison.

## Analysis Workflow

The program:

1. Reads the reference and sample sequences from the FASTA file.
2. Compares each sample with the reference sequence.
3. Detects nucleotide substitutions.
4. Records mutation positions and nucleotide changes.
5. Generates summary tables.
6. Saves results as Excel files.
7. Creates plots showing mutation patterns.

## Results

The analysis produced:

* **446 total mutation records**
* **78 unique mutation rows**
* **6 samples with length mismatches**

The generated results include:

* Master mutation table
* Mutation summary
* Sample summary
* Position summary
* Mutation type summary
* Mutation details
* Length mismatch samples

## Visualizations

The project also generates figures showing:

* Distribution of mutations
* Top samples according to mutation counts

Figures are stored in the `figures/` directory.

## Technologies Used

* Python
* Pandas
* Biopython
* Matplotlib
* FASTA sequence data
* Excel output

## Project Structure

```text
DNA-Variation-Analyzer/
│
├── ATM.fasta
├── dna_variation_analyzer.py
├── requirements.txt
├── .gitignore
│
├── figures/
│   ├── distribution_of_mutation.png
│   └── top_10_samples.png
│
└── results/
    ├── length_mismatch_samples.xlsx
    ├── master_mutation_table.xlsx
    ├── mutation_details.xlsx
    ├── mutation_summary.xlsx
    ├── mutation_type_summary.xlsx
    ├── position_summary.xlsx
    └── sample_summary.xlsx
```

## How to Run

Create and activate a Python environment, then install the required packages:

```bash
pip install -r requirements.txt
```

Run the analyzer using:

```bash
python dna_variation_analyzer.py ATM.fasta MK284930 --all
```

The analysis results will be saved in the `results/` directory and figures will be saved in the `figures/` directory.

## Skills Demonstrated

* Python programming
* FASTA sequence handling
* Sequence comparison
* Mutation identification
* Data processing with Pandas
* Biological data analysis
* Excel report generation
* Data visualization
* Basic command-line usage
* Organizing bioinformatics analysis workflows

## Future Improvements

Possible future improvements include:

* Additional sequence quality checks
* More detailed mutation classification
* Integration with NCBI sequence information
* Phylogenetic analysis
* Additional mutation visualizations
* Command-line options for more flexible analysis
