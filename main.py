import numpy as np
from scipy.interpolate import make_smoothing_spline

x = np.array([0.0000000,
0.0560852,
0.1121700,
0.1682560,
0.2243410,
0.2804260,
0.3365110,
0.3925970,
0.4486820,
0.5047670,
0.5608520
])

y = np.array([
2.06471,
-1.21856,
4.96948,
-1.40129,
-3.41988,
-0.605902,
2.38111,
0.653869,
4.82372,
-9.43329,
6.00037
])


spl = make_smoothing_spline(x, y)
for i in range(len(x)):
    print(spl(x[i]))