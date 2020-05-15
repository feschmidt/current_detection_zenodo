from src.tlineformulas import getlinepars
# import numpy as np
# import scipy.constants as const
# eps0 = const.epsilon_0
# c0 = const.c
# print()

# all units in SI

print('### part 1: reproducing the values from our existing device ###')

geo = {'s': 10e-6, 'w': 6e-6, 't': 20e-9, 'london': 158e-9, 'epsr': 11.9}
tlpars = getlinepars(geo)
[print(str(key)+':', val) for key, val in tlpars.items()]
print()

print('### part 2: device scaled down x10 ###')

geo = {'s': 1e-6, 'w': .6e-6, 't': 20e-9, 'london': 158e-9, 'epsr': 11.9}
tlpars = getlinepars(geo)
[print(str(key)+':', val) for key, val in tlpars.items()]
print()

print('### part 2: device scaled down x10, but thicker ###')

geo = {'s': 1e-6, 'w': .6e-6, 't': 70e-9, 'london': 158e-9, 'epsr': 11.9}
tlpars = getlinepars(geo)
[print(str(key)+':', val) for key, val in tlpars.items()]
print()

print('### part 3: playground ###')

geo = {'s': 1e-6, 'w': .6e-6, 't': 80e-9, 'london': 158e-9, 'epsr': 11.9}
tlpars = getlinepars(geo)
[print(str(key)+':', val) for key, val in tlpars.items()]
print()
