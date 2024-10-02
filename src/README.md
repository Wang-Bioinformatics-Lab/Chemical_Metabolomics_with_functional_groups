# Running the codes
Make sure the correct conda enviroment located at `{project root}/bin/conda_env.yml` is installed.

```bash
conda env create -f {project root}/bin/conda_env.yml
```

# Input Format Management
To Create an MGF file with correct formating, the required fields has to exist in the mgf file. To check if the MGF file is correct run the `format_checker.py` with your mgf file:

```bash
conda activate Chemical_Metablomics
python format_checker.py {path to your mgf file}
``` 

If your mgf file does not contain the reactivity information, make sure to add that using the guideline provided at `{project root}/README.md`.

If the known compounds information are not added to the MGF file, you can use the `select_known_scans.py` and passing the knowns csv and your mgf file to get the updated mgf.
the knowns file has to be a CSV with the following column names:
`mgf_path,Name,SMILES,FGs,Scan,Adduct`

* if the known are scattered for different MGF files for different reactions, you can use `convert_knowns_to_correct_format.py` to get the merged knowns file. make sure that the knowns have the following columns: 
`Sample name,Sample,Educt ID H+,Educt ID Na+,Name,SMILE,FGs`
sample name should have the format `{folder}_{reaction}`.
Sample should be the subfolder if exists (Pool, Box) or same as folder


# Convert .ms files to mgf files
To Create MGF files from .ms files, first, activate the conda enviroment provided in the `bin` folder at the root of the project:

```bash
conda activate Chemical_Metablomics
```

Then, run the `ms_to_mgf.py` passing the directory with the ms files and the csv file containing the functional groups. An example for the functional groups file is provided in `{project root}/examples/example_fgs.csv`