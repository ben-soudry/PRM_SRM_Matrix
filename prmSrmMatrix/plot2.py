import numpy as np
import matplotlib.pyplot as plt

'''PRM_plus_SRM = np.array([[1666185, 7062635, 11891040, 3887532, 2887270, 2897970, 4504396, 27275849, 57356795, 30500369, 11178168, 6017909, 3239127, 1744225, 873830, 469820, 250365, 137331, 78608, 46804, 27645, 16640, 10566, 6385, 4184, 2595, 1644, 1084, 745, 1152], 
[25642, 108119, 178075, 62041, 48565, 45957, 58554, 378762, 808740, 390386, 197350, 128109, 83060, 51327, 30513, 18821, 11142, 7102, 4305, 2863, 1658, 1022, 627, 347, 226, 139, 94, 79, 37, 42]])
 '''
'''
PRM_plus_SRM =  np.array([[2052930, 9405310, 15903578, 5232126, 3863114, 3914868, 6553555, 36145919, 71060041, 37033511, 13700710, 7412540, 4007420, 2159105, 1078994, 577432, 304743, 165866, 94166, 55294, 32251, 19257, 12074, 7279, 4740, 2902, 1839, 1225, 827, 1282], [31057, 141751, 236270, 81521, 62960, 59917, 84244, 494205, 993209, 469614, 240194, 157021, 101618, 63209, 37809, 23294, 13791, 8643, 5091, 3316, 1933, 1167, 748, 396, 265, 154, 106, 85, 44, 45]])
'''
'''
PRM_plus_SRM = np.array([[2530447, 12463748, 21561317, 7261226, 5171758, 5249505, 9250505, 48753951, 91713036, 47312584, 17236518, 9370400, 5102645, 2744626, 1370592, 729110, 381134, 204259, 113944, 66014, 37926, 22361, 13867, 8266, 5381, 3249, 2061, 1369, 911, 1410],
[38176, 184536, 317440, 110414, 82338, 76819, 117504, 654628, 1271297, 591359, 299017, 196542, 127374, 79467, 48026, 29685, 17500, 10691, 6301, 3908, 2391, 1426, 889, 471, 311, 182, 116, 92, 48, 45]])
'''

'''PRM_plus_SRM = np.array([[1261286, 5203575, 8611392, 2884726, 2090723, 2131876, 2938574, 19890336, 45657412, 24998765, 8983517, 4809899, 2577320, 1384260, 694846, 375633, 201980, 111625, 64649, 38793, 23334, 14066, 9027, 5520, 3656, 2287, 1457, 971, 664, 1041] , 
[19746, 80025, 130235, 47075, 36258, 34730, 38513, 281099, 646903, 320625, 159130, 102437, 67091, 41217, 24396, 15215, 8914, 5730, 3586, 2465, 1429, 856, 529, 289, 192, 120, 74, 72, 31, 40]])'''


'''PRM_plus_SRM = np.array([[1261286, 5203575, 8611392, 2884726, 2090723, 2131876, 2938574, 19890336, 45657412, 24998765, 8983517, 4809899, 2577320, 1384260, 694846, 375633, 201980, 111625, 64649, 38793, 23334, 14066, 9027, 5520, 3656, 2287, 1457, 971, 664, 1041], [19746, 80025, 130235, 47075, 36258, 34730, 38513, 281099, 646903, 320625, 159130, 102437, 67091, 41217, 24396, 15215, 8914, 5730, 3586, 2465, 1429, 856, 529, 289, 192, 120, 74, 72, 31, 40]])'''

'''PRM_plus_SRM = np.array ([[512354, 1903910, 2959825, 1049872, 746318, 816896, 886075, 6764140, 20951673, 12865063, 4214089, 2214902, 1166424, 619730, 310814, 171267, 94846, 53636, 32728, 20376, 12922, 7945, 5345, 3357, 2251, 1423, 915, 642, 446, 706], [8261, 30002, 47477, 17913, 13485, 13628, 11598, 102195, 301457, 164954, 75748, 47868, 31805, 19430, 11239, 7051, 4123, 2752, 1892, 1364, 805, 445, 271, 155, 103, 67, 31, 28, 16, 25]])'''


'''PRM_plus_SRM = np.array([[293030, 849432, 1417789, 503884, 365551, 390147, 372801, 3107801, 11593991, 7810377, 2363112, 1228183, 637294, 335544, 167435, 93096, 52872, 30062, 19026, 12144, 7954, 5008, 3388, 2142, 1486, 944, 596, 438, 305, 488], [4697, 13838, 23436, 8732, 6609, 6700, 5471, 48716, 167840, 100440, 43424, 27236, 18126, 11238, 6488, 4003, 2324, 1498, 1051, 693, 412, 242, 150, 92, 62, 31, 18, 18, 12, 14]])'''


'''PRM_plus_SRM = np.array([[21538388, 3026047, 1405621, 4475533, 23095422, 8288518, 9264476, 71081113, 106593651, 52906478, 17537844, 11159601, 6985979, 4581706, 3006135, 1972362, 1303949, 836506, 544749, 353604, 231717, 162935, 112278, 75420, 47524, 31487, 23085, 16734, 11619, 28952], [161190, 22952, 11297, 33078, 177933, 57967, 69355, 480222, 765038, 293329, 143092, 113029, 83944, 65034, 48872, 35022, 30217, 20134, 14182, 10577, 7273, 5730, 4482, 3416, 2275, 1607, 1092, 965, 628, 1779]])'''


'''PRM_plus_SRM = np.array([[0, 0, 3191, 2020971, 21799670, 7096354, 8098402, 50873983, 38589503, 22814084, 7991595, 5371247, 3380157, 2325063, 1498375, 1067890, 781896, 534528, 358652, 236730, 158315, 114748, 78750, 52498, 32284, 20401, 14978, 11113, 7782, 20092], [0, 0, 18, 12801, 166230, 46942, 58092, 340822, 258793, 127168, 65236, 55819, 43043, 35803, 25150, 21056, 19487, 13579, 9674, 7561, 5119, 4282, 3456, 2621, 1754, 1248, 824, 762, 484, 1496]])
'''


'''PRM_plus_SRM = np.array([[21538388, 3026047, 1402430, 2454562, 1295752, 1192164, 1166074, 20207130, 68004148, 30092394, 9546249, 5788354, 3605822, 2256643, 1507760, 904472, 522053, 301978, 186097, 116874, 73402, 48187, 33528, 22922, 15240, 11086, 8107, 5621, 3837, 8860], [ 

161190, 22952, 11279, 20277, 11703, 11025, 11263, 139400, 506245, 166161, 77856, 57210, 40901, 29231, 23722, 13966, 10730, 6555, 4508, 3016, 2154, 1448, 1026, 795, 521, 359, 268, 203, 144, 283]])'''

PRM_plus_SRM = np.array([[2888535, 308066, 164197, 408339, 2774819, 915015, 882148, 11425580, 20393652, 14056078, 3423855, 2240256, 1348678, 890283, 575073, 377261, 256919, 162409, 110553, 68840, 46303, 34286, 24289, 17142, 10545, 7558, 5693, 4472, 3226, 8635], [23133, 2796, 1507, 3335, 23214, 7419, 7511, 79073, 154029, 79255, 28735, 24028, 16737, 14283, 10826, 7601, 7048, 4425, 3307, 2409, 1632, 1362, 1068, 749, 457, 326, 248, 153, 139, 352]])


PRM_plus_SRM = np.array([[812550, 96991, 49888, 105650, 802164, 231618, 257227, 4296270, 7439948, 6290470, 1269699, 839892, 492904, 331835, 213195, 138975, 94304, 59767, 41069, 24634, 16483, 12292, 9084, 6510, 3993, 2863, 2174, 1712, 1218, 3425, 6857, 877, 448, 884, 6706, 2110, 2279, 27973, 56975, 37866, 11146, 10238, 6753, 6593, 4709, 3208, 2908, 1971, 1516, 1015, 676, 560, 419, 254, 163, 109, 78, 54, 53, 114]])

noPeaks = np.sum(PRM_plus_SRM[0])
peaks = np.sum(PRM_plus_SRM[1])
peakProb = peaks / (peaks+noPeaks)

print("PRM+SRM peakProb ", peakProb )

x_axis = [-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0,
               0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,
               10.0,11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]


ProbNoPeak = []



for i in range(0, PRM_plus_SRM.shape[1]):
    colProb = PRM_plus_SRM[0][i] / noPeaks; 
    ProbNoPeak.append(colProb  ) 
    print(colProb, end =", ")
print()

plt.plot(x_axis, ProbNoPeak, label = 'Score Probability given no peak')

ProbPeak = []

for i in range(0, PRM_plus_SRM.shape[1]):
    colProb = PRM_plus_SRM[1][i] / peaks; 
    ProbPeak.append(colProb) 
    print(colProb, end =", ")
print()


plt.plot(x_axis, ProbPeak, label = 'Score Probability given peak')


PRM_colProbs = []
plt.xlabel('MSGF+ Score')
plt.ylabel('Probability of Score')

plt.legend()
plt.show()

#print(PRM_plus_SRM)