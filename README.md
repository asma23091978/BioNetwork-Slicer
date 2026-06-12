# BioNetwork-Slicer

**BioNetwork-Slicer** is a Python–Prolog workflow for the analysis of biological transition networks.  
It combines graph-based attractor analysis with dynamic slicing in BioReSolve to identify the contexts, reactions, and molecular entities associated with selected biological target states.

The project was developed for qualitative analysis of the phenylalanine–tyrosine (Phe–Tyr) catabolic pathway using Reaction Systems and Discrete Entity-Level Annotation (DEA).

---

## Features

- Target-oriented attractor analysis of biological transition networks.
- Identification of initial contexts leading to selected target attractors.
- Dynamic negative slicing through BioReSolve and SWI-Prolog.
- Extraction of positive and negative dependencies associated with target outcomes.
- Export of results as CSV and JSON files.
- Documentation of the modeling assumptions and analysis workflow.

---

## Repository structure

```text
BioNetwork-Slicer/
├── README.md
├── requirements.txt
├── .gitignore
├── slicing_script.py
├── BioReSolvePositive.pl
├── output/
│   └── hgd0_04-02.dot
├── prolog/
│   └── hgd0_04_02_reactions.pl
├── docs/
│   └── methods.tex
├── figures/
│   └── software architecture.png
└── scripts/
    └── repo_check.py
```

---

## Requirements

Before running the project, make sure you have:

- Python 3.9 or later
- SWI-Prolog installed and available from the command line
- The Python packages listed in `requirements.txt`
- The BioReSolve Prolog file `BioReSolvePositive.pl`
- The labeled transition system graph file `output/hgd0_04-02.dot`

---

## Installation

Clone the repository:

```bash
git clone https://github.com/asma23091978/BioNetwork-Slicer.git
cd BioNetwork-Slicer
```

Create a virtual environment:

```bash
python -m venv .venv
```

Activate it:

```bash
source .venv/bin/activate
```

On Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install the Python dependencies:

```bash
pip install -r requirements.txt
```

---

## Required input files

The main Python script expects the following files:

```text
BioReSolvePositive.pl
output/hgd0_04-02.dot
```

`BioReSolvePositive.pl` must be located in the root directory because the Python script loads it directly through SWI-Prolog.

The graph file must be placed inside the `output/` folder and must be named:

```text
hgd0_04-02.dot
```

---

## How to run

After installing the dependencies and placing all required files in the correct locations, run:

```bash
python slicing_script.py
```

The script performs:

1. loading of the `.dot` transition graph;
2. attractor detection;
3. target-oriented context extraction;
4. generation of `params.pl`;
5. dynamic negative slicing through BioReSolve;
6. export of CSV and JSON results.

---

## Output

The analysis generates result files under:

```text
Results/
```

Typical outputs include:

- context-to-target CSV files;
- sliced computation JSON files;
- positive and negative dependency sets;
- summary tables for downstream analysis.

Generated files such as `Results/`, `params.pl`, and temporary slicing files are excluded from version control through `.gitignore`.

---

## Methodological background

A detailed description of the modeling assumptions and analysis workflow is provided in [`docs/methods.tex`](docs/methods.tex). The method covers the construction of the Reaction Systems model for Phe--Tyr catabolism, the Discrete Entity-Level Annotation framework, generation of the labeled transition system in BioReSolve, target-oriented attractor analysis, and dynamic negative slicing.

---

## Notes on reproducibility

The repository contains the Python workflow and the Prolog reaction-system model files required for the computational analysis.  
To fully reproduce the analysis, the user must provide the complete BioReSolve Prolog engine file and the corresponding `.dot` transition graph generated from the model.

---

## Citation

If you use this repository in academic work, please cite the associated manuscript or publication describing the method and the biological model.

---

## License

A license has not been selected yet. Before making the repository public, add an appropriate open-source license, such as MIT, Apache-2.0, or GPL-3.0, depending on how you want others to use and redistribute the code.
