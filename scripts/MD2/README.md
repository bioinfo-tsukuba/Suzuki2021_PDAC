# Scripts to run ExtractEdges.py and DiffEdges.py for P_low_grade and P_high_grade groups

- 1. Make a pair of emFile and annFile for each group
  - `010-make_NATMI_input.R`
  - Note: Given `N` is the number of patients, you will get `N` pairs of emFile and annFile.
- 2. Run NATMI ExtractEdges.py on a pair of emFile and annFile for each group
  - `020-run_NATMI_ExtractEdges.sh`
  - Note: Given `N` is the number of patients, you will get `N` output direcotries.
- 3. Gun NATMI DiffEdges.py for a pair of P_low_grade and P_high_grade groups
  - `030-run_NATMI_DiffEdges.sh`
  - Note: Given `N` is the number of patients, you will get `N` PDF files.
