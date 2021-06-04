# Author: adpozuelo@gmail.com
# Version: 1.0
# Date: 06/2021
# Input files: nsim_times.mc

import numpy as np

data = np.loadtxt('nsim_times.mc', dtype=np.float64, max_rows=10)

mu = np.mean(data)
sigma = np.std(data)

print('\u03BC: %e' % mu)
print('2\u03C3: %e' % (2.0 * sigma))