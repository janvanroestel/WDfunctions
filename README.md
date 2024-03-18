# WDfunctions
a collection of useful functions related to white dwarf stars.

Currently contains two mass-radius functions; one by Eggleton and one by Naunenberg. In addition, there is also a function that interpolates values from a modelgrid (http://evolgroup.fcaglp.unlp.edu.ar/TRACKS/newtables.html)

# Install
Download and install manually, or use pip:
```
pip install git+https://github.com/janvanroestel/WDfunctions.git
```

# Example

```
import numpy as np
import matplotlib.pyplot as plt
import WDfunctions

# grid of masses in solar units
mrange = np.linspace(0.2,1.3,1000)

plt.plot(mrange,
  WDfunctions.WD_MTR(mrange,np.log10(15000)*np.ones_like(mrange)),)
plt.plot(mrange,
  WDfunctions.WD_MR_Eggleton(mrange),'k:')
plt.plot(mrange,
  WDfunctions.WD_MR_Nauenberg(mrange),'k--')
plt.xlabel('Radius')
plt.xlabel('Mass')
plt.show()
```
