# af2genomics

## Installation
```
$ git clone git@github.com:jurgjn/af2genomics.git
$ cd af2genomics
$ pip install -e .
```

Set up environment:
```
mamba env create --name af2genomics-env --file envs/af2genomics.yaml
```

Alternative/minimal (python>3.9 needed for functools.cache):
```
mamba create -n af2genomics-env "python>3.9" numpy scipy matplotlib seaborn requests
mamba activate af2genomics-env
pip install ipykernel
python -m ipykernel install --user --name=af2genomics-env
```
