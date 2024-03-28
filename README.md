# af2genomics

## Installation
```
$ git clone git@github.com:jurgjn/af2genomics.git
$ cd af2genomics
$ pip install -e .
```

Environment set up if needed (python>3.9 needed for functools.cache):
```
mamba create -n af2genomics "python>3.9" numpy scipy matplotlib seaborn requests
mamba activate af2genomics
pip install ipykernel
python -m ipykernel install --user --name=af2genomics
```
    