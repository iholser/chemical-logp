# Prediction of logP (Kow) values from chemical structures

## Running
To reliable run using rdkit, a docker container with a compiled version of rdkit was created. 

Start the server
```sh
docker run -p 8000:8000 ianh/chem-notebook python3 /srv/server.py
```

The static version of the front end may be opened in any browser, or can be run in react dev mode

# Development

## Data Wrangling
[Data Collection](experiments/data_collection.ipynb) was done by starting with an initial set of structures pulled from [ZINC](https://zinc.docking.org/). This dataset was supplimented with structures of simple molecules that are not usually treated as drugs, such as alkanes and, polycyclic aromatic hydrocarbons. These structures were matched with records in PubChem using [Pug REST](https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest) apis. The pubchem ids were then used to fetch experimental data where available. Data where no pubchem record could be matched were dropped as the structure information alone is not meaningful for this project. LogP values were backfilled populated by taking the experimental values where available. If not available, values were backfilled using xLogP3. Values that were still missing were backfilled with xLogP3 AA values. If no logP value was available, the record was dropped.

[Data Processing](experiments/data_wrangling.ipynb) was done to clean up the data set and generate 3D structures using [RDKit](https://rdkit.org) and create fingerprint/conformational bits using [E3FP](https://github.com/keiserlab/e3fp). RDKit structures were droped from the dataframe so that models could be developed without the dependency on RDKit or Docker.


## Model Seclection

[Model Selection](experiments/model_selection.ipynb) was done using a random 2000 sample subset of the data. Different regressors were tested and the most promising based on mean absolute error and R² values were explored further. Experiments were also done to optimize fingerprint folding. Folded fingerprints have a higher data density, but more abstraction. In this experiment it was observed that a single fold greatly improved model performance and further folding caused it to gradually decrease. For this reason, all further experiments were run with a single fold to a bit length of 1024.

## Hyperparameter Tuning

Parameter hypertuning and evalation was performed on [Support Vector Regression](), [Random Forest]() and, [XGBoost]().
Additionally, a [neural network]() was constructed and evaluated.

The best accuracy and performance was obtained using Support Vector Regression. The neural network acheived a comparable accuracy, but not as efficiently.
