[View on Colab](https://colab.research.google.com/github/mayhd3/DIPC/blob/master/topinsulators.ipynb)

(Press Ctrl+F9 and scroll down.)

### Requirements:

[C2DB](https://cmr.fysik.dtu.dk/c2db/c2db.html)

#### For running cells outside of the jupyter kernel:

Ubuntu LTS:
```
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt install python3.9 python3.9-dev libxc4 libxc-dev libopenblas-base libopenblas-dev
python3.9 -m pip install --upgrade setuptools pip distlib
```

Python:
`python3 -m pip install --upgrade --user ase asr gpaw`

also Numpy, Pandas, Matplotlib, Tabulate, Pickle, and Jupyter Notebook

WSL:
[VcXsrv](https://sourceforge.net/projects/vcxsrv/files/latest/download)
