# Ligand Parametrization

## Overview
This is a Python package designed to provide a simple workflow for ligand parameterization. It automates many of the 
key features encountered by users, including...

## Documentation

The online documentation is located here: https://ligandparam.readthedocs.io/en/latest/



### Developing

#### Releasing a new version

1. update `version` at `pyproject.toml`. Eg: `0.3.2`
2. Add and commit the change. Then tag the commit: `git tag 0.3.2`
3. `git push origin --tags`

gitlab's CI/CD will publish the new version automatically
