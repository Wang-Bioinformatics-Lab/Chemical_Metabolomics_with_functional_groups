## Installation

You will need to have conda, mamba, and nextflow installed to run things locally. 

### GNPS2 Workflow Input information

Check the definition for the workflow input and display parameters:
https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/workflowdev/

### Checking out the deployment submodules

use the following commands from the deploy_gnps2 folder. 

You might need to checkout the module, do this by running

```
git submodule init
git submodule update --init --recursive
```

## Input Formatting
The input MGF File has to prepared with the same criteria as for the library search workflow, However, for the scans with reactivity, the `ONLINE_REACTIVITY` key needs to be added as well.
The format used to report the reactivity must be an array of dictionary type items described here:
```json
[
    {
        "compound_type": "Educt", //['Product', 'Educt']
        "reaction_name": "Hydroxylamine", // name of the reaction
        "filename_contains": "Hydroxylamine", // ['cysteine', 'AQC', 'Hydroxylamine']
        "functional_group_name": "Aldheides", // name or an indicator for the functional group smart.
        "educt_smarts": "O=[C]([#1])[C]([#1])([#1])C", //smart formula of the educt
        "reaction_smarts": "O=[C]([#1])[C]([#1])([#1])C",
        "delta_mz": 15.0109,
        "linked_ids": [
            5415
        ]
    },
    // next reactions
]
```


To perform analysis, known compounds must be provided, the known compounds need to have the following keys:

* `SMILES`: Contains thr structure of the known compound.
* `NAME (Optional)`: The name of the compound.
* `ADDUCT (Optional)`: The adduct of the compound.
* `FEATURE_ID ('Optional')`: This will be shown in the final data in case a refrence back to the original data is needed.

### Format converter
if your mgf input does not have the correct formatting, you can use the codes provided in the `src` directory to convert to correct format or check the format.

### Example Input
Check the `examples` folder for the `example_input.mgf` for an input example.

## Running the project
To run the project with the library search

```
nextflow run nf_workflow.nf --inputspectra {path to input spectra mgfs} --inputlibraries {path to input libraries} --with_library_search 1
```

If you already have the result of the library search and you want to skip that

```
nextflow run nf_workflow.nf --inputspectra {path to input spectra mgfs} --library_search_result {path to library search file} --with_library_search 0
```

