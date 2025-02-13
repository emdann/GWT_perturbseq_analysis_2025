# GWT_perturbseq_analysis
Analysis of genome-wide perturb-seq screen on primary T cells

## Contents

- `src` - analysis code
    - `1_preprocess/` - ingest and preprocess new experiments from cellranger outputs
- `metadata` - sample and experimental metadata, configs, gene annotations etc

## Set-up compute environment

To install required packages (including [perturbseq_tools module](https://github.com/emdann/perturbseq_tools))

```
conda env create -f environment.yml
conda activate gwt-env
```

## Contributing code

Copy directory in your local environment 
```
git clone git@github.com:emdann/GWT_perturbseq_analysis.git
```

To add new code/notebooks, add them to `src` in an appropriate directory, then run:
```
git add src/path/to/new_script.py
git commit -m 'adding script for analysis x'
git push origin master
```

To edit existing code/notebooks without creating conflicts, make edits, then run:
```
git checkout -b new-branch-name
git add .
git commit -m 'esiting script x y and z'
git push origin new-branch-name
```






