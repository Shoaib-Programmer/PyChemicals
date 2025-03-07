# PyChemicals 😁

This is a python package I created right after finishing CS50P, to help people with chemical calculations like pH, and plotting titration curves. Chempy exists, but it doesn't
handle these sort of issues, so that is why I decided to make my own.

Collab would be greatly appreciated.

## How to setup up
This package is not production ready. So installing it in development mode is recommended.

First, let's clone the repo:
```bash
git clone https://github.com/Shoaib-Programmer/PyChemicals.git
```

`cd` into PyChemicals, then install the required dependencies:
```bash
pip install -r requirements.txt
```

I am considering implementing Poetry for this project to avoid this extra step for setup, and make it production-ready.

Then install it in development mode:
```bash
pip install -e .
```

And you're all set! You can import PyChemicals into python file anywhere on your computer. Here is an example of how you can go about using it:
```python
from PyChemicals import Acid

hcl = Acid('Hydrochloric Acid')
print(hcl.ph())
```

I'm considering implementing an SQLite database to store predefined chemicals, and use IUPAC names instead of common names.
