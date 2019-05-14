import numpy as np
import matplotlib.pyplot as plt

PRM_plus_SRM = np.array([[1666185, 7062635, 11891040, 3887532, 2887270, 2897970, 4504396, 27275849, 57356795, 30500369, 11178168, 6017909, 3239127, 1744225, 873830, 469820, 250365, 137331, 78608, 46804, 27645, 16640, 10566, 6385, 4184, 2595, 1644, 1084, 745, 1152], 
[25642, 108119, 178075, 62041, 48565, 45957, 58554, 378762, 808740, 390386, 197350, 128109, 83060, 51327, 30513, 18821, 11142, 7102, 4305, 2863, 1658, 1022, 627, 347, 226, 139, 94, 79, 37, 42]])

noPeaks = np.sum(PRM_plus_SRM[0])
peaks = np.sum(PRM_plus_SRM[1])
peakProb = peaks / (peaks+noPeaks)

print("PRM+SRM peakProb ", peakProb )

PRM = np.array([[21538388, 3026047, 1402430, 2454562, 1295752, 1192164, 1166074, 20207130, 68004148, 30092394, 9546249, 5788354, 3605822, 2256643, 1507760, 904472, 522053, 301978, 186097, 116874, 73402, 48187, 33528, 22922, 15240, 11086, 8107, 5621, 3837, 8860], [161190, 22952, 11279, 20277, 11703, 11025, 11263, 139400, 506245, 166161, 77856, 57210, 40901, 29231, 23722, 13966, 10730, 6555, 4508, 3016, 2154, 1448, 1026, 795, 521, 359, 268, 203, 144, 283]])

noPeaks = np.sum(PRM[0])
peaks = np.sum(PRM[1])
peakProb = peaks / (peaks+noPeaks)


print("PRM peakProb ", peakProb )

SRM = np.array([[0, 0, 3191, 2020971, 21799670, 7096354, 8098402, 50873983, 38589503, 22814084, 7991595, 5371247, 3380157, 2325063, 1498375, 1067890, 781896, 534528, 358652, 236730, 158315, 114748, 78750, 52498, 32284, 20401, 14978, 11113, 7782, 20092] ,[0, 0, 18, 12801, 166230, 46942, 58092, 340822, 258793, 127168, 65236, 55819, 43043, 35803, 25150, 21056, 19487, 13579, 9674, 7561, 5119, 4282, 3456, 2621, 1754, 1248, 824, 762, 484, 1496]])

noPeaks = np.sum(SRM[0])
peaks = np.sum(SRM[1])
peakProb = peaks / (peaks+noPeaks)


print("PRM peakProb ", peakProb )


x_axis = [-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0,
               0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,
               10.0,11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]


PRM_plus_SRM_colProbs = []

for i in range(0, PRM_plus_SRM.shape[1]):
    if(PRM_plus_SRM[0][i] + PRM_plus_SRM[1][i] != 0):
        colProb = PRM_plus_SRM[1][i] / (PRM_plus_SRM[0][i] + PRM_plus_SRM[1][i])
    else:
        colProb = 0
    print(colProb, end =", ")
    PRM_plus_SRM_colProbs.append(colProb)
print()
plt.plot(x_axis, PRM_plus_SRM_colProbs, label = 'PRM plus SRM')

PRM_colProbs = []

for i in range(0, PRM.shape[1]):
    if(PRM[0][i] + PRM[1][i] != 0):
        colProb = PRM[1][i] / (PRM[0][i] + PRM[1][i])
    else:
        colProb = 0
    print(colProb, end =", ")
    PRM_colProbs.append(colProb)
print()
plt.plot(x_axis, PRM_colProbs, label = 'PRM')

SRM_colProbs = []

for i in range(0, SRM.shape[1]):
    if(SRM[0][i] + SRM[1][i] != 0):
        colProb = SRM[1][i] / (SRM[0][i] + SRM[1][i])
    else:
        colProb = 0
    print(colProb, end =", ")
    SRM_colProbs.append(colProb)
print()
plt.plot(x_axis, SRM_colProbs, label = 'SRM')


plt.xlabel('MSGF+ Score')
plt.ylabel('Probability of spectra peak')

plt.legend()
plt.show()

#print(PRM_plus_SRM)
