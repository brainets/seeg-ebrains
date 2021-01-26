# Python envrironement config

Required librairies : 
- NumPy
- Xarray
- Matploblib
- MNE-Python
- [Neo-python](https://github.com/NeuralEnsemble/python-neo)

I found an issue for reading ELAN's data. Hence, I made some changes to Neo and then proposed a PR to fix the issue
(see my [PR](https://github.com/NeuralEnsemble/python-neo/pull/935)). For the moment I'm still waiting for merging the modifications to Neo.
In the meantime, if you want to be able to load the files, you need to install my version of neo :

```bash
pip install git+https://github.com/EtienneCmb/python-neo
```
