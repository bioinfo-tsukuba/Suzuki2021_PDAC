# Scripts to run ExtractEdges.py for each patient

- 1. Make a pair of emFile and annFile for each patient
  - `010-make_NATMI_input.R`
  - Note: Given `N` is the number of patients, you will get `N` pairs of emFile and annFile.
- 2. Run NATMI on a pair of emFile and annFile for each patient
  - `020-run_NATMI_each_patient.sh`
  - Note: Given `N` is the number of patients, you will get `N` output direcotries.
- 3. Get heatmap where x-axis represents ligand-receptor pairs and y-axix represents cell-type pairs.
  - `030-make_heatmap_each_patient.R`
  - Note: Given `N` is the number of patients, you will get `N` PDF files.
