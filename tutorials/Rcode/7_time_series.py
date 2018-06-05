import pandas as pd
import pylab as pl

import statsmodels.api as sm

pd.read_csv("HD3651_rv.dat", header=None, sep=r"\s*", names=['Day', 'Vel', 'sigVel', 'tmp'])
rv = pd.read_csv("HD3651_rv.dat", header=None, sep=r"\s*", names=['Day', 'Vel', 'sigVel', 'tmp'])

ax = rv.plot(x='Day', y='Vel')
ax.errorbar(rv['Day'], rv['Vel'], yerr=rv['sigVel'], fmt='.')
pl.text(10800, 30, 'HD 3651    P=62.257d   e=0.63', ha='right')
pl.show()

rv['Daynorm'] = rv.Day % 62.257
#(rv[,2] ~ lp( % 62.257, nn=0.5), weights=1/rv[,3]^2)

locfit_phase = sp.interpolate.UnivariateSpline(rv['Daynorm'][np.argsort(rv['Daynorm'])],
                                               rv['Vel'][np.argsort(rv['Daynorm'])],
                                               w=1.0/rv['sigVel'][np.argsort(rv['Daynorm'])]**2,
                                               s=50)

ax.plot(rv['Daynorm'][np.argsort(rv['Daynorm'])],
        spl(rv['Daynorm'][np.argsort(rv['Daynorm'])]), '-')

pl.ylim(-30,20)
     #, band='local',
pl.xlabel('JD - 2440000 mod P=62.257d')
pl.ylabel('Velocity (m/s)')
pl.show()

   
import scipy.signal as signal
ax = pl.figure().add_subplot(111)
f = np.arange(0.001,0.0501,0.001)
pgram = signal.lombscargle(rv.Day, rv.Vel, f, normalize=True)
ax.plot(f, pgram)
pl.show()
