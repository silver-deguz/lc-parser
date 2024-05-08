## Setup Environment
Create a python environment and install necessary dependencies.
```
python3 -m venv venv
source venv/bin/activate
python3 -m pip install requirements.txt
```

## Run script
```
python3 main.py -d <path to a single .txt file or a flat folder of .txt files> --smooth <optional: whether to smooth data> --save-plots <optional: whether to save figures to output folder on disk> 
```