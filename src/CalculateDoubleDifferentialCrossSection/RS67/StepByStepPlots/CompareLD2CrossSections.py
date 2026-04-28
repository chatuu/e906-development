import uproot
import numpy as np
import ROOT
import sys
import os
import math

# Set ROOT batch mode to avoid opening windows
ROOT.gROOT.SetBatch(True)
# Remove the statistics box
ROOT.gStyle.SetOptStat(0)
# Set color palette
ROOT.gStyle.SetPalette(ROOT.kBird)

# ==========================================
# 0. Hardcoded Shivangi LD2 Data
# ==========================================
# Format: "xF_bin": [ [mass_min, mass_max, xsec, stat_err, sys_err_low, sys_err_high], ... ]
SHIVANGI_DATA = {
    "0.0-0.05": [
        [4.2, 4.5, 0.000000, 0.000000, 0.000000, 0.000000],
        [4.5, 4.8, 6.981630, 21.029500, 5.645100, 6.479820],
        [4.8, 5.1, 8.469470, 4.128100, 1.881600, 0.967758],
        [5.1, 5.4, 3.796700, 2.412110, 0.919370, 0.713182],
        [5.4, 5.7, 2.979740, 1.094640, 0.564971, 0.326210],
        [5.7, 6.0, 2.172180, 0.637578, 0.321968, 0.165598],
        [6.0, 6.3, 1.878820, 0.577268, 0.300901, 0.169372],
        [6.3, 6.6, 1.743900, 0.362662, 0.241713, 0.128480],
        [6.6, 6.9, 0.669427, 0.207319, 0.114798, 0.045106],
        [6.9, 7.2, 0.832782, 0.254222, 0.124520, 0.058104],
        [7.2, 7.5, 0.750034, 0.274599, 0.097649, 0.058388],
        [7.5, 7.8, 0.117855, 0.126358, 0.023396, 0.008710],
        [7.8, 8.1, -0.006599, 0.006625, 0.009332, 0.009332],
        [8.1, 8.4, 0.000000, 0.000000, 0.000000, 0.000000],
        [8.4, 8.7, 0.000000, 0.000000, 0.000000, 0.000000]
    ],
    "0.05-0.1": [
        [4.2, 4.5, 73.641700, 79.316400, 8.799050, 6.227190],
        [4.5, 4.8, 3.443820, 10.813400, 3.972070, 3.589540],
        [4.8, 5.1, 7.500030, 2.386270, 1.818010, 0.978490],
        [5.1, 5.4, 6.530890, 1.253780, 1.064250, 0.456672],
        [5.4, 5.7, 5.396090, 0.839860, 0.788197, 0.374682],
        [5.7, 6.0, 3.408610, 0.661027, 0.469841, 0.270061],
        [6.0, 6.3, 1.706820, 0.367789, 0.282820, 0.099183],
        [6.3, 6.6, 1.864610, 0.311974, 0.250431, 0.141118],
        [6.6, 6.9, 1.114220, 0.233166, 0.143045, 0.087882],
        [6.9, 7.2, 0.333038, 0.276485, 0.064875, 0.046357],
        [7.2, 7.5, 0.325292, 0.150757, 0.046568, 0.023367],
        [7.5, 7.8, 0.381444, 0.201633, 0.064572, 0.025667],
        [7.8, 8.1, 0.217261, 0.166724, 0.046110, 0.017452],
        [8.1, 8.4, 0.000000, 0.000000, 0.000000, 0.000000],
        [8.4, 8.7, 0.000000, 0.000000, 0.000000, 0.000000]
    ],
    "0.1-0.15": [
        [4.2, 4.5, 88.741100, 51.468200, 12.368200, 8.388610],
        [4.5, 4.8, 8.148940, 4.118380, 2.947340, 2.352450],
        [4.8, 5.1, 9.401260, 1.977410, 1.782850, 1.054710],
        [5.1, 5.4, 7.410540, 1.105090, 1.178680, 0.581015],
        [5.4, 5.7, 3.729020, 0.711054, 0.616052, 0.352300],
        [5.7, 6.0, 2.359510, 0.530481, 0.410399, 0.186288],
        [6.0, 6.3, 2.512210, 0.342960, 0.388923, 0.162203],
        [6.3, 6.6, 1.223720, 0.281131, 0.182256, 0.082056],
        [6.6, 6.9, 1.688640, 0.272056, 0.226370, 0.128013],
        [6.9, 7.2, 0.775320, 0.190602, 0.106397, 0.057597],
        [7.2, 7.5, 0.414149, 0.151742, 0.058246, 0.030156],
        [7.5, 7.8, -0.008720, 0.004389, 0.012332, 0.012332],
        [7.8, 8.1, 0.079242, 0.085332, 0.016350, 0.006130],
        [8.1, 8.4, 0.283884, 0.208679, 0.042447, 0.019807],
        [8.4, 8.7, 0.000000, 0.000000, 0.000000, 0.000000]
    ],
    "0.15-0.2": [
        [4.2, 4.5, 11.442100, 4.807990, 3.479870, 1.528250],
        [4.5, 4.8, 14.256800, 2.548230, 2.601890, 1.290640],
        [4.8, 5.1, 10.803100, 1.262550, 1.745920, 0.666489],
        [5.1, 5.4, 7.895290, 0.861893, 1.243790, 0.499019],
        [5.4, 5.7, 4.885780, 0.615830, 0.730210, 0.345509],
        [5.7, 6.0, 3.150490, 0.432706, 0.497097, 0.201732],
        [6.0, 6.3, 2.066170, 0.346122, 0.303426, 0.140223],
        [6.3, 6.6, 1.795420, 0.266748, 0.251647, 0.131087],
        [6.6, 6.9, 1.143970, 0.203887, 0.155552, 0.085652],
        [6.9, 7.2, 0.551947, 0.181500, 0.090457, 0.037186],
        [7.2, 7.5, 0.202207, 0.169812, 0.030842, 0.013954],
        [7.5, 7.8, 0.451037, 0.163141, 0.056857, 0.036198],
        [7.8, 8.1, 0.078792, 0.081407, 0.011720, 0.005515],
        [8.1, 8.4, 0.104492, 0.104701, 0.011610, 0.009499],
        [8.4, 8.7, 0.000000, 0.000000, 0.000000, 0.000000]
    ],
    "0.2-0.25": [
        [4.2, 4.5, 15.763800, 4.800060, 3.121440, 1.535270],
        [4.5, 4.8, 14.613000, 2.040810, 2.336120, 1.506720],
        [4.8, 5.1, 9.928910, 1.149570, 1.618620, 0.789370],
        [5.1, 5.4, 7.453230, 0.760436, 1.107360, 0.536569],
        [5.4, 5.7, 5.027120, 0.516373, 0.735864, 0.345279],
        [5.7, 6.0, 3.057280, 0.409108, 0.461974, 0.215135],
        [6.0, 6.3, 2.616430, 0.303814, 0.352523, 0.197466],
        [6.3, 6.6, 1.548610, 0.216608, 0.231501, 0.108063],
        [6.6, 6.9, 0.596618, 0.144817, 0.090971, 0.041178],
        [6.9, 7.2, 0.622141, 0.141942, 0.089412, 0.044563],
        [7.2, 7.5, 0.338635, 0.155707, 0.038128, 0.030380],
        [7.5, 7.8, 0.180351, 0.097171, 0.036607, 0.013671],
        [7.8, 8.1, 0.139775, 0.100467, 0.017640, 0.011205],
        [8.1, 8.4, 0.000000, 0.000000, 0.000000, 0.000000],
        [8.4, 8.7, -0.004113, 0.004123, 0.005817, 0.005817]
    ],
    "0.25-0.3": [
        [4.2, 4.5, 16.000800, 2.265740, 2.947040, 1.552680],
        [4.5, 4.8, 12.833900, 1.424580, 2.230010, 1.049460],
        [4.8, 5.1, 9.069020, 0.860092, 1.380850, 0.649653],
        [5.1, 5.4, 6.380560, 0.576810, 0.961338, 0.425091],
        [5.4, 5.7, 4.601670, 0.438541, 0.671629, 0.313771],
        [5.7, 6.0, 3.337990, 0.289372, 0.496995, 0.233488],
        [6.0, 6.3, 2.302090, 0.258824, 0.333007, 0.159940],
        [6.3, 6.6, 1.418580, 0.203797, 0.198848, 0.103565],
        [6.6, 6.9, 0.674957, 0.199071, 0.109121, 0.041675],
        [6.9, 7.2, 0.755658, 0.146970, 0.107505, 0.054536],
        [7.2, 7.5, 0.585872, 0.138340, 0.074548, 0.046601],
        [7.5, 7.8, 0.238431, 0.099016, 0.029485, 0.019493],
        [7.8, 8.1, 0.056041, 0.056105, 0.006227, 0.005095],
        [8.1, 8.4, 0.000000, 0.000000, 0.000000, 0.000000],
        [8.4, 8.7, 0.000000, 0.000000, 0.000000, 0.000000]
    ],
    "0.3-0.35": [
        [4.2, 4.5, 18.414400, 1.903580, 2.990350, 1.441910],
        [4.5, 4.8, 13.109200, 1.149120, 1.973460, 0.929509],
        [4.8, 5.1, 9.353010, 0.758214, 1.421380, 0.641271],
        [5.1, 5.4, 5.786790, 0.471191, 0.877552, 0.376938],
        [5.4, 5.7, 4.024580, 0.382738, 0.608252, 0.263445],
        [5.7, 6.0, 3.268790, 0.314715, 0.453884, 0.237468],
        [6.0, 6.3, 2.153780, 0.234741, 0.313102, 0.148997],
        [6.3, 6.6, 1.593330, 0.202644, 0.225343, 0.115513],
        [6.6, 6.9, 1.059040, 0.161117, 0.158719, 0.073789],
        [6.9, 7.2, 0.473227, 0.133559, 0.066389, 0.034526],
        [7.2, 7.5, 0.317486, 0.132835, 0.037267, 0.027324],
        [7.5, 7.8, 0.275755, 0.149961, 0.034674, 0.022184],
        [7.8, 8.1, 0.125477, 0.076587, 0.022342, 0.008542],
        [8.1, 8.4, 0.182039, 0.109284, 0.028058, 0.012498],
        [8.4, 8.7, 0.097122, 0.097291, 0.010791, 0.008829]
    ],
    "0.35-0.4": [
        [4.2, 4.5, 16.114900, 1.567490, 2.704230, 1.383150],
        [4.5, 4.8, 11.222300, 0.888990, 1.757090, 0.773069],
        [4.8, 5.1, 9.554130, 0.676913, 1.363180, 0.673738],
        [5.1, 5.4, 6.461790, 0.445358, 0.949548, 0.443987],
        [5.4, 5.7, 3.587320, 0.314927, 0.499415, 0.258065],
        [5.7, 6.0, 2.736480, 0.254792, 0.400135, 0.189523],
        [6.0, 6.3, 1.846420, 0.200720, 0.263223, 0.133052],
        [6.3, 6.6, 1.259450, 0.184697, 0.179919, 0.090613],
        [6.6, 6.9, 0.968037, 0.141874, 0.137057, 0.070122],
        [6.9, 7.2, 0.461695, 0.100452, 0.067696, 0.032612],
        [7.2, 7.5, 0.513777, 0.121965, 0.070488, 0.038176],
        [7.5, 7.8, 0.368536, 0.111049, 0.056487, 0.025370],
        [7.8, 8.1, 0.087097, 0.064399, 0.014359, 0.005864],
        [8.1, 8.4, 0.113105, 0.081454, 0.014686, 0.008827],
        [8.4, 8.7, -0.002512, 0.002515, 0.003552, 0.003552]
    ],
    "0.4-0.45": [
        [4.2, 4.5, 13.315500, 1.100600, 2.090480, 1.021850],
        [4.5, 4.8, 10.638900, 0.704015, 1.609040, 0.695027],
        [4.8, 5.1, 7.983320, 0.547969, 1.206930, 0.521560],
        [5.1, 5.4, 5.246620, 0.372582, 0.755911, 0.366467],
        [5.4, 5.7, 3.814650, 0.312167, 0.552637, 0.263393],
        [5.7, 6.0, 2.774400, 0.264613, 0.377381, 0.204742],
        [6.0, 6.3, 1.817150, 0.208020, 0.278669, 0.119891],
        [6.3, 6.6, 1.377800, 0.165198, 0.201448, 0.097507],
        [6.6, 6.9, 0.754537, 0.141176, 0.104935, 0.055437],
        [6.9, 7.2, 0.353849, 0.123757, 0.045349, 0.027955],
        [7.2, 7.5, 0.359050, 0.120580, 0.047786, 0.027395],
        [7.5, 7.8, 0.115598, 0.059850, 0.017312, 0.008058],
        [7.8, 8.1, 0.031045, 0.114122, 0.002441, 0.005360],
        [8.1, 8.4, -0.003027, 0.002144, 0.004281, 0.004281],
        [8.4, 8.7, 0.070176, 0.070262, 0.007797, 0.006380]
    ],
    "0.45-0.5": [
        [4.2, 4.5, 12.590200, 0.925469, 1.853700, 0.895362],
        [4.5, 4.8, 9.653780, 0.648818, 1.441630, 0.649382],
        [4.8, 5.1, 6.547850, 0.415911, 0.967275, 0.442238],
        [5.1, 5.4, 4.799760, 0.343633, 0.662839, 0.348554],
        [5.4, 5.7, 2.963270, 0.261030, 0.431269, 0.202562],
        [5.7, 6.0, 2.345380, 0.201434, 0.330750, 0.170415],
        [6.0, 6.3, 1.640920, 0.176018, 0.221838, 0.123477],
        [6.3, 6.6, 1.188460, 0.146243, 0.168498, 0.085996],
        [6.6, 6.9, 0.547846, 0.111213, 0.073965, 0.041273],
        [6.9, 7.2, 0.575059, 0.108145, 0.084583, 0.040536],
        [7.2, 7.5, 0.307998, 0.086352, 0.049329, 0.020848],
        [7.5, 7.8, 0.143668, 0.067031, 0.022894, 0.009739],
        [7.8, 8.1, 0.136669, 0.070267, 0.019165, 0.009975],
        [8.1, 8.4, 0.139864, 0.082652, 0.018886, 0.010536],
        [8.4, 8.7, 0.000000, 0.000000, 0.000000, 0.000000]
    ],
    "0.5-0.55": [
        [4.2, 4.5, 10.757000, 0.739267, 1.677310, 0.681624],
        [4.5, 4.8, 8.061030, 0.529713, 1.124740, 0.580037],
        [4.8, 5.1, 5.132000, 0.331645, 0.730985, 0.365221],
        [5.1, 5.4, 3.516390, 0.268716, 0.523406, 0.235650],
        [5.4, 5.7, 2.375990, 0.215388, 0.348266, 0.162761],
        [5.7, 6.0, 1.670830, 0.183768, 0.232404, 0.120767],
        [6.0, 6.3, 1.395110, 0.157982, 0.188091, 0.105231],
        [6.3, 6.6, 0.718617, 0.114704, 0.107841, 0.050032],
        [6.6, 6.9, 0.767950, 0.116245, 0.097774, 0.061049],
        [6.9, 7.2, 0.433715, 0.093191, 0.057164, 0.033384],
        [7.2, 7.5, 0.337953, 0.089629, 0.044600, 0.025983],
        [7.5, 7.8, 0.124972, 0.059034, 0.022142, 0.008491],
        [7.8, 8.1, 0.172587, 0.078953, 0.022949, 0.013179],
        [8.1, 8.4, 0.039928, 0.042496, 0.007630, 0.002843],
        [8.4, 8.7, 0.111459, 0.080350, 0.014637, 0.008608]
    ],
    "0.55-0.6": [
        [4.2, 4.5, 8.154150, 0.559218, 1.207230, 0.553150],
        [4.5, 4.8, 5.102610, 0.375346, 0.799865, 0.318574],
        [4.8, 5.1, 4.468730, 0.303053, 0.616293, 0.325795],
        [5.1, 5.4, 2.564620, 0.220668, 0.389129, 0.167328],
        [5.4, 5.7, 2.075650, 0.171809, 0.292195, 0.151026],
        [5.7, 6.0, 1.102600, 0.156854, 0.158541, 0.076590],
        [6.0, 6.3, 0.922566, 0.143240, 0.124812, 0.069378],
        [6.3, 6.6, 0.598144, 0.101774, 0.085581, 0.042985],
        [6.6, 6.9, 0.483784, 0.092035, 0.068780, 0.034933],
        [6.9, 7.2, 0.350973, 0.081189, 0.048919, 0.025740],
        [7.2, 7.5, 0.148791, 0.058872, 0.024190, 0.010037],
        [7.5, 7.8, 0.169350, 0.065182, 0.021479, 0.013511],
        [7.8, 8.1, 0.061771, 0.044406, 0.007953, 0.004859],
        [8.1, 8.4, 0.036234, 0.037559, 0.005590, 0.002487],
        [8.4, 8.7, -0.003250, 0.002303, 0.004596, 0.004596]
    ],
    "0.6-0.65": [
        [4.2, 4.5, 6.552570, 0.440831, 0.938395, 0.456657],
        [4.5, 4.8, 4.745350, 0.327720, 0.664606, 0.339056],
        [4.8, 5.1, 2.587740, 0.204594, 0.380556, 0.177361],
        [5.1, 5.4, 2.025360, 0.169657, 0.286132, 0.146958],
        [5.4, 5.7, 1.330700, 0.132961, 0.190636, 0.095539],
        [5.7, 6.0, 0.940483, 0.122042, 0.135802, 0.067137],
        [6.0, 6.3, 0.775924, 0.114262, 0.107003, 0.057406],
        [6.3, 6.6, 0.511816, 0.085808, 0.073361, 0.036732],
        [6.6, 6.9, 0.362222, 0.087886, 0.051803, 0.026039],
        [6.9, 7.2, 0.132603, 0.048958, 0.021151, 0.008986],
        [7.2, 7.5, 0.124017, 0.052139, 0.017546, 0.008988],
        [7.5, 7.8, 0.046598, 0.035172, 0.009080, 0.003378],
        [7.8, 8.1, 0.057957, 0.041671, 0.007474, 0.004552],
        [8.1, 8.4, 0.108811, 0.062976, 0.012090, 0.009892],
        [8.4, 8.7, 0.050495, 0.052314, 0.007745, 0.003475]
    ],
    "0.65-0.7": [
        [4.2, 4.5, 4.336390, 0.340513, 0.607913, 0.313180],
        [4.5, 4.8, 2.633900, 0.214143, 0.389845, 0.177071],
        [4.8, 5.1, 1.913380, 0.176877, 0.268825, 0.136434],
        [5.1, 5.4, 1.188750, 0.126864, 0.185984, 0.081094],
        [5.4, 5.7, 1.051250, 0.112181, 0.156003, 0.073683],
        [5.7, 6.0, 0.678980, 0.099614, 0.092745, 0.050639],
        [6.0, 6.3, 0.371705, 0.083415, 0.057021, 0.025578],
        [6.3, 6.6, 0.229882, 0.071712, 0.035420, 0.015786],
        [6.6, 6.9, 0.198536, 0.057036, 0.028795, 0.014129],
        [6.9, 7.2, 0.058161, 0.033387, 0.017759, 0.008613],
        [7.2, 7.5, 0.081935, 0.042332, 0.012099, 0.005761],
        [7.5, 7.8, 0.138990, 0.058498, 0.019885, 0.009989],
        [7.8, 8.1, 0.031426, 0.031445, 0.003492, 0.002857],
        [8.1, 8.4, 0.049122, 0.049163, 0.005458, 0.004466],
        [8.4, 8.7, 0.000000, 0.000000, 0.000000, 0.000000]
    ],
    "0.7-0.75": [
        [4.2, 4.5, 2.523080, 0.228897, 0.353579, 0.180163],
        [4.5, 4.8, 1.637320, 0.170258, 0.248101, 0.106626],
        [4.8, 5.1, 1.190920, 0.134728, 0.172774, 0.082625],
        [5.1, 5.4, 0.799661, 0.093516, 0.112390, 0.058258],
        [5.4, 5.7, 0.607620, 0.084819, 0.085535, 0.044211],
        [5.7, 6.0, 0.436197, 0.072154, 0.061328, 0.031769],
        [6.0, 6.3, 0.212258, 0.052910, 0.035930, 0.014282],
        [6.3, 6.6, 0.111894, 0.054338, 0.022534, 0.008403],
        [6.6, 6.9, 0.072894, 0.054740, 0.010138, 0.005355],
        [6.9, 7.2, 0.049863, 0.029177, 0.006165, 0.004077],
        [7.2, 7.5, 0.068319, 0.040261, 0.009099, 0.005209],
        [7.5, 7.8, 0.028912, 0.029960, 0.004452, 0.001986],
        [7.8, 8.1, 0.063265, 0.044849, 0.007029, 0.005751],
        [8.1, 8.4, -0.000354, 0.000355, 0.000500, 0.000500],
        [8.4, 8.7, 0.000000, 0.000000, 0.000000, 0.000000]
    ],
    "0.75-0.8": [
        [4.2, 4.5, 1.442870, 0.150774, 0.208504, 0.100840],
        [4.5, 4.8, 0.810564, 0.121057, 0.110761, 0.059513],
        [4.8, 5.1, 0.642336, 0.097489, 0.088905, 0.046501],
        [5.1, 5.4, 0.595968, 0.081022, 0.081301, 0.044497],
        [5.4, 5.7, 0.378078, 0.070492, 0.050674, 0.028666],
        [5.7, 6.0, 0.234967, 0.052291, 0.034916, 0.016455],
        [6.0, 6.3, 0.167463, 0.046776, 0.026090, 0.011443],
        [6.3, 6.6, 0.139475, 0.045053, 0.017992, 0.010951],
        [6.6, 6.9, 0.111844, 0.046802, 0.015186, 0.008385],
        [6.9, 7.2, 0.035692, 0.026669, 0.006407, 0.002438],
        [7.2, 7.5, -0.000886, 0.000629, 0.001253, 0.001253],
        [7.5, 7.8, 0.002265, 0.001173, 0.000252, 0.000206],
        [7.8, 8.1, 0.000000, 0.000000, 0.000000, 0.000000],
        [8.1, 8.4, 0.000000, 0.000000, 0.000000, 0.000000],
        [8.4, 8.7, 0.000000, 0.000000, 0.000000, 0.000000]
    ]
}


# ==========================================
# 1. Define Physics Constants
# ==========================================
PROTONS_ON_TARGET_LH2 = 1.57319e+17
PROTONS_ON_TARGET_LD2 = 7.51113e+16
PROTONS_ON_TARGET_FLASK = 3.57904e+16

LH2_TARGET_DENSITY_MOL_CM2 = 3.597
AVOGADRO_CONSTANT = 6.022e23
PROTONS_PER_NUCLEON_LH2 = 1.008
XF_BIN_WIDTH = 0.05
THD_THH_RATIO = 0.1084/3.5966 # Ratio of target thicknesses T_HD / T_HH

# Beam Attenuation
val_exp = -(50.8 * 0.0708) / 52.0
BEAM_ATTENUATION = (52.0 / (0.0708 * 50.8)) * (1.0 - math.exp(val_exp))

# Global Normalization Constants
GLOBAL_CONSTANT_LH2 = (PROTONS_PER_NUCLEON_LH2 * 1e33) / (
LH2_TARGET_DENSITY_MOL_CM2 * AVOGADRO_CONSTANT * PROTONS_ON_TARGET_LH2 * BEAM_ATTENUATION * XF_BIN_WIDTH
)

GLOBAL_CONSTANT_LD2 = (PROTONS_PER_NUCLEON_LH2 * 1e33) / (
LH2_TARGET_DENSITY_MOL_CM2 * AVOGADRO_CONSTANT * PROTONS_ON_TARGET_LD2 * BEAM_ATTENUATION * XF_BIN_WIDTH
)

# Flask Normalization Factors
FLASK_NORM_LH2 = PROTONS_ON_TARGET_LH2 / PROTONS_ON_TARGET_FLASK
FLASK_NORM_LD2 = PROTONS_ON_TARGET_LD2 / PROTONS_ON_TARGET_FLASK

print("=== Physics Constants ===")
print(f"PoT LH2: {PROTONS_ON_TARGET_LH2:.4e}")
print(f"PoT LD2: {PROTONS_ON_TARGET_LD2:.4e}")
print(f"PoT Flask: {PROTONS_ON_TARGET_FLASK:.4e}")
print(f"Flask Norm LH2: {FLASK_NORM_LH2:.4f}")
print(f"Flask Norm LD2: {FLASK_NORM_LD2:.4f}")
print(f"Global Constant LH2: {GLOBAL_CONSTANT_LH2:.4e}")
print(f"Global Constant LD2: {GLOBAL_CONSTANT_LD2:.4e}")
print("=========================")


# ==========================================
# 2. Define the User's Cut Function
# ==========================================
def e906_chuck_cuts(tree, cut=4.2, beam_offset=1.6):
    events = tree.arrays(library="np")
    class EventNamespace:
        def __init__(self, data):
            self.__dict__.update(data)
    e = EventNamespace(events)

    dimuon_cut = (
    (np.abs(e.dx) < 0.25) & (np.abs(e.dy - beam_offset) < 0.22) &
    (e.dz < -5.) & (e.dz > -280.) &
    (np.abs(e.dpx) < 1.8) & (np.abs(e.dpy) < 2.0) &
    (e.dpx**2 + e.dpy**2 < 5.) &
    (e.dpz < 116.) & (e.dpz > 38.) &
    (e.mass > cut) & (e.mass < 8.8) &
    (e.dx**2 + (e.dy - beam_offset)**2 < 0.06) &
    (e.xF < 0.95) & (e.xF > -0.1) &
    (e.xT > 0.05) & (e.xT <= 0.58) &
    (np.abs(e.costh) < 0.5) &
    (np.abs(e.trackSeparation) < 270.) &
    (e.chisq_dimuon < 18)
    )

    track1_cut = (
    (e.chisq1_target < 15.) & (e.pz1_st1 > 9.) & (e.pz1_st1 < 75.) &
    (e.nHits1 > 13) & (e.x1_t**2 + (e.y1_t - beam_offset)**2 < 320.) &
    (e.x1_d**2 + (e.y1_d - beam_offset)**2 < 1100.) &
    (e.x1_d**2 + (e.y1_d - beam_offset)**2 > 16.) &
    (e.chisq1_target < 1.5 * e.chisq1_upstream) &
    (e.chisq1_target < 1.5 * e.chisq1_dump) &
    (e.z1_v < -5.) & (e.z1_v > -320.) &
    (e.chisq1/(e.nHits1 - 5) < 12) & ((e.y1_st1)/(e.y1_st3 ) < 1.) &
    (np.abs(np.abs(e.px1_st1 - e.px1_st3) - 0.416) < 0.008) &
    (np.abs(e.py1_st1 - e.py1_st3) < 0.008) &
    (np.abs(e.pz1_st1 - e.pz1_st3) < 0.08) &
    ((e.y1_st1) * (e.y1_st3) > 0.) & (np.abs(e.py1_st1) > 0.02)
    )

    track2_cut = (
    (e.chisq2_target < 15.) & (e.pz2_st1 > 9.) & (e.pz2_st1 < 75.) &
    (e.nHits2 > 13) & (e.x2_t**2 + (e.y2_t - beam_offset)**2 < 320.) &
    (e.x2_d**2 + (e.y2_d - beam_offset)**2 < 1100.) &
    (e.x2_d**2 + (e.y2_d - beam_offset)**2 > 16.) &
    (e.chisq2_target < 1.5 * e.chisq2_upstream) &
    (e.chisq2_target < 1.5 * e.chisq2_dump) &
    (e.z2_v < -5.) & (e.z2_v > -320.) &
    (e.chisq2/(e.nHits2 - 5) < 12) & ((e.y2_st1 )/(e.y2_st3) < 1.) &
    (np.abs(np.abs(e.px2_st1 - e.px2_st3) - 0.416) < 0.008) &
    (np.abs(e.py2_st1 - e.py2_st3) < 0.008) &
    (np.abs(e.pz2_st1 - e.pz2_st3) < 0.08) &
    ((e.y2_st1) * (e.y2_st3) > 0.) & (np.abs(e.py2_st1) > 0.02)
    )

    tracks_cut = (
    (np.abs(e.chisq1_target + e.chisq2_target - e.chisq_dimuon) < 2.) &
    ((e.y1_st3) * (e.y2_st3) < 0.) &
    (e.nHits1 + e.nHits2 > 29) & (e.nHits1St1 + e.nHits2St1 > 8) &
    (np.abs(e.x1_st1 + e.x2_st1) < 42)
    )

    occ_cut = (
    (e.D1 < 400) & (e.D2 < 400) & (e.D3 < 400) &
    (e.D1 + e.D2 + e.D3 < 1000)
    )

    total_cut_mask = (track1_cut & track2_cut & tracks_cut & dimuon_cut & occ_cut)

    filtered_events = {}
    for key, val in events.items():
        filtered_events[key] = val[total_cut_mask]
    return filtered_events


# ==========================================
# 3. Setup Binning and Histograms
# ==========================================
mass_bins_np = np.array([4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7], dtype=float)
xf_bins_np = np.round(np.arange(0.0, 0.85, 0.05), 2)

def create_histograms(file_label):
    hists = {}
    def make_th2(name, title):
        h = ROOT.TH2D(f"{name}_{file_label}", f"{title} ({file_label});Mass [GeV];x_{{F}}",
        len(mass_bins_np)-1, mass_bins_np,
        len(xf_bins_np)-1, xf_bins_np)
        h.Sumw2()
        h.SetStats(0)
        h.GetXaxis().CenterTitle()
        h.GetYaxis().CenterTitle()
        h.GetXaxis().SetTitleOffset(1.2)
        h.GetYaxis().SetTitleOffset(1.2)
        return h

    hists["Y_total"] = make_th2("Y_total", "Total Yield (result)")
    hists["Y_mix"] = make_th2("Y_mix", "Mix Yield (result_mix)")
    hists["E_total_reco"] = make_th2("E_total_reco", "Avg Reco Eff (Total)")
    hists["E_mix_reco"] = make_th2("E_mix_reco", "Avg Reco Eff (Mix)")
    hists["E_total_hodo"] = make_th2("E_total_hodo", "Avg Hodo Eff (Total)")
    hists["E_mix_hodo"] = make_th2("E_mix_hodo", "Avg Hodo Eff (Mix)")
    hists["E_total_final"] = make_th2("E_total_final", "Avg Final Eff (Total)")
    hists["E_mix_final"] = make_th2("E_mix_final", "Avg Final Eff (Mix)")
    hists["E_final_signal"] = make_th2("E_final_signal", "Avg Signal Efficiency")
    hists["Y_corrected"] = make_th2("Y_corrected", "Corrected Yield (Total Error)")
    hists["Y_corrected_stat"] = make_th2("Y_corrected_stat", "Corrected Yield (Stat Error)")
    hists["Y_corrected_sys"] = make_th2("Y_corrected_sys", "Corrected Yield (Sys Error)")

    return hists


# ==========================================
# 4. Calculation Helpers
# ==========================================
def get_weighted_mean_and_error(eff_arr, err_arr):
    if len(eff_arr) == 0: return 0.0, 0.0
    mean = np.mean(eff_arr)
    err_on_mean = np.sqrt(np.sum(err_arr**2)) / len(eff_arr)
    return mean, err_on_mean

def add_latex_to_bin(hist, x_center, y_center, value, error):
    val_f = float(value)
    err_f = float(error)
    if np.isnan(val_f) or np.isinf(val_f): val_f = 0.0
    if np.isnan(err_f) or np.isinf(err_f): err_f = 0.0
    latex_text = f"#splitline{{{val_f:.3f}}}{{#pm {err_f:.3f}}}"
    l = ROOT.TLatex(x_center, y_center, latex_text)
    l.SetTextSize(0.015)
    l.SetTextAlign(22)
    l.SetTextColor(ROOT.kBlack)
    hist.GetListOfFunctions().Add(l)


# ==========================================
# 5. Basic Processing
# ==========================================
def process_file_data(filename, file_label, output_file):
    print(f"\nProcessing {filename} as {file_label}...")
    try: f = uproot.open(filename)
    except Exception as e:
        print(f"Error opening {filename}: {e}")
        return {}

    output_file.cd()
    hists = create_histograms(file_label)
    try:
        data_total = e906_chuck_cuts(f["result"])
        data_mix = e906_chuck_cuts(f["result_mix"])
    except KeyError as e:
        print(f" Missing TTree in file: {e}")
        return {}

    for i_m in range(len(mass_bins_np) - 1):
        m_low = mass_bins_np[i_m]
        m_high = mass_bins_np[i_m+1]
        m_center = (m_low + m_high) / 2.0
        for i_x in range(len(xf_bins_np) - 1):
            x_low = xf_bins_np[i_x]
            x_high = xf_bins_np[i_x+1]
            x_center = (x_low + x_high) / 2.0
            mask_tot = (data_total["mass"] >= m_low) & (data_total["mass"] < m_high) & \
            (data_total["xF"] >= x_low) & (data_total["xF"] < x_high)
            mask_mix = (data_mix["mass"] >= m_low) & (data_mix["mass"] < m_high) & \
            (data_mix["xF"] >= x_low) & (data_mix["xF"] < x_high)
            N_tot = np.sum(mask_tot)
            N_mix = np.sum(mask_mix)
            err_N_tot = np.sqrt(N_tot)
            err_N_mix = np.sqrt(N_mix)

            eff_reco_tot = data_total["recoeff"][mask_tot]
            err_reco_tot = data_total["recoeff_error"][mask_tot]
            eff_hodo_tot = data_total["hodoeff"][mask_tot]
            err_hodo_tot = data_total["hodoeff_error"][mask_tot]
            eff_final_tot = eff_reco_tot * eff_hodo_tot
            err_final_tot = np.sqrt( (eff_hodo_tot * err_reco_tot)**2 + (eff_reco_tot * err_hodo_tot)**2 )

            eff_reco_mix = data_mix["recoeff"][mask_mix]
            err_reco_mix = data_mix["recoeff_error"][mask_mix]
            eff_hodo_mix = data_mix["hodoeff"][mask_mix]
            err_hodo_mix = data_mix["hodoeff_error"][mask_mix]
            eff_final_mix = eff_reco_mix * eff_hodo_mix
            err_final_mix = np.sqrt( (eff_hodo_mix * err_reco_mix)**2 + (eff_reco_mix * err_hodo_mix)**2 )

            mean_reco_tot, e_mean_reco_tot = get_weighted_mean_and_error(eff_reco_tot, err_reco_tot)
            mean_hodo_tot, e_mean_hodo_tot = get_weighted_mean_and_error(eff_hodo_tot, err_hodo_tot)
            mean_final_tot, e_mean_final_tot = get_weighted_mean_and_error(eff_final_tot, err_final_tot)
            mean_reco_mix, e_mean_reco_mix = get_weighted_mean_and_error(eff_reco_mix, err_reco_mix)
            mean_hodo_mix, e_mean_hodo_mix = get_weighted_mean_and_error(eff_hodo_mix, err_hodo_mix)
            mean_final_mix, e_mean_final_mix = get_weighted_mean_and_error(eff_final_mix, err_final_mix)

            # Initialize variables BEFORE the conditional block to prevent UnboundLocalError
            diff_yield = N_tot - N_mix
            val_sig_eff = 0.0
            err_sig_eff = 0.0
            val_corr_yield = 0.0
            err_corr_stat = 0.0
            err_corr_sys = 0.0
            err_corr_total = 0.0
            
            if diff_yield > 0:
                numerator = (N_tot * mean_final_tot) - (N_mix * mean_final_mix)
                val_sig_eff = numerator / diff_yield
                term1 = (N_tot * e_mean_final_tot)**2
                term2 = (N_mix * e_mean_final_mix)**2
                err_sig_eff = (1.0 / diff_yield) * np.sqrt(term1 + term2)

            if val_sig_eff > 0 and diff_yield > 0:
                val_corr_yield = diff_yield / val_sig_eff
                err_corr_stat = np.sqrt(N_tot + N_mix) / val_sig_eff
                err_corr_sys = val_corr_yield * (err_sig_eff / val_sig_eff)
                err_corr_total = np.sqrt(err_corr_stat**2 + err_corr_sys**2)

            def fill_bin(key, val, err):
                hists[key].SetBinContent(i_m + 1, i_x + 1, val)
                hists[key].SetBinError(i_m + 1, i_x + 1, err)
                if "stat" not in key and "sys" not in key and "Centroid" not in key:
                    add_latex_to_bin(hists[key], m_center, x_center, val, err)

            fill_bin("Y_total", N_tot, err_N_tot)
            fill_bin("Y_mix", N_mix, err_N_mix)
            fill_bin("E_total_reco", mean_reco_tot, e_mean_reco_tot)
            fill_bin("E_mix_reco", mean_reco_mix, e_mean_reco_mix)
            fill_bin("E_total_hodo", mean_hodo_tot, e_mean_hodo_tot)
            fill_bin("E_mix_hodo", mean_hodo_mix, e_mean_hodo_mix)
            fill_bin("E_total_final", mean_final_tot, e_mean_final_tot)
            fill_bin("E_mix_final", mean_final_mix, e_mean_final_mix)
            fill_bin("E_final_signal", val_sig_eff, err_sig_eff)
            fill_bin("Y_corrected", val_corr_yield, err_corr_total)
            fill_bin("Y_corrected_stat", val_corr_yield, err_corr_stat)
            fill_bin("Y_corrected_sys", val_corr_yield, err_corr_sys)

    print(f" Saving PDFs for {file_label}...")
    for name, hist in hists.items():
        if "stat" in name or "sys" in name or "Centroid" in name:
            continue
        c_name = f"c_{name}_{file_label}"
        c = ROOT.TCanvas(c_name, c_name, 1200, 900)
        c.SetRightMargin(0.15)
        c.SetLeftMargin(0.12)
        c.SetBottomMargin(0.12)
        c.SetTickx(1)
        c.SetTicky(1)
        local_min = sys.float_info.max
        local_max = -sys.float_info.max
        has_data = False
        for ix in range(1, hist.GetNbinsX() + 1):
            for iy in range(1, hist.GetNbinsY() + 1):
                val = hist.GetBinContent(ix, iy)
                if val != 0.0:
                    has_data = True
                    if val < local_min: local_min = val
                    if val > local_max: local_max = val
        if has_data:
            hist.SetMinimum(local_min)
            hist.SetMaximum(local_max)
        else:
            hist.SetMinimum(0.0)
            hist.SetMaximum(1.0)
        hist.Draw("COLZ")
        c.SaveAs(f"{name}_{file_label}.pdf")
        c.Close()

    return hists


# ==========================================
# 6. Subtracted Yield Calculations
# ==========================================
def generate_subtracted_plot(hists_target, hists_flask, flask_norm, target_label, output_file):
    print(f"\nGenerating Subtracted Plot: Y_corrected_{target_label} - Y_corrected_Flask...")
    if "Y_corrected_stat" not in hists_target or "Y_corrected_stat" not in hists_flask:
        print("Error: Y_corrected_stat histogram missing.")
        return None

    output_file.cd()
    h_target_stat = hists_target["Y_corrected_stat"]
    h_target_sys = hists_target["Y_corrected_sys"]
    h_flask_stat = hists_flask["Y_corrected_stat"]
    h_flask_sys = hists_flask["Y_corrected_sys"]
    name = f"Y_corrected_Subtracted_{target_label}"
    title = f"Corrected Yield ({target_label} - Flask)"
    h_sub = ROOT.TH2D(name, f"{title};Mass [GeV];x_{{F}}",
    len(mass_bins_np)-1, mass_bins_np,
    len(xf_bins_np)-1, xf_bins_np)
    h_sub.Sumw2()
    h_sub.SetStats(0)
    h_sub.GetXaxis().CenterTitle()
    h_sub.GetYaxis().CenterTitle()
    h_sub.GetXaxis().SetTitleOffset(1.2)
    h_sub.GetYaxis().SetTitleOffset(1.2)

    h_sub_stat = h_sub.Clone(f"{name}_stat")
    h_sub_sys = h_sub.Clone(f"{name}_sys")

    for i_m in range(len(mass_bins_np) - 1):
        m_center = (mass_bins_np[i_m] + mass_bins_np[i_m+1]) / 2.0
        bin_x = i_m + 1
        for i_x in range(len(xf_bins_np) - 1):
            x_center = (xf_bins_np[i_x] + xf_bins_np[i_x+1]) / 2.0
            bin_y = i_x + 1
            y_target = h_target_stat.GetBinContent(bin_x, bin_y)
            e_target_stat = h_target_stat.GetBinError(bin_x, bin_y)
            e_target_sys = h_target_sys.GetBinError(bin_x, bin_y)
            y_flask = h_flask_stat.GetBinContent(bin_x, bin_y)
            e_flask_stat = h_flask_stat.GetBinError(bin_x, bin_y)
            e_flask_sys = h_flask_sys.GetBinError(bin_x, bin_y)
            val_sub = y_target - (flask_norm * y_flask)
            err_sub_stat = np.sqrt(e_target_stat**2 + (flask_norm * e_flask_stat)**2)
            err_sub_sys = np.sqrt(e_target_sys**2 + (flask_norm * e_flask_sys)**2)
            err_sub_total = np.sqrt(err_sub_stat**2 + err_sub_sys**2)
            h_sub.SetBinContent(bin_x, bin_y, val_sub)
            h_sub.SetBinError(bin_x, bin_y, err_sub_total)
            add_latex_to_bin(h_sub, m_center, x_center, val_sub, err_sub_total)

            h_sub_stat.SetBinContent(bin_x, bin_y, val_sub)
            h_sub_stat.SetBinError(bin_x, bin_y, err_sub_stat)

            h_sub_sys.SetBinContent(bin_x, bin_y, val_sub)
            h_sub_sys.SetBinError(bin_x, bin_y, err_sub_sys)

    c = ROOT.TCanvas(f"c_{name}", f"c_{name}", 1200, 900)
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.SetTickx(1)
    c.SetTicky(1)
    local_min = sys.float_info.max
    local_max = -sys.float_info.max
    has_data = False
    for ix in range(1, h_sub.GetNbinsX() + 1):
        for iy in range(1, h_sub.GetNbinsY() + 1):
            val = h_sub.GetBinContent(ix, iy)
            if val != 0.0:
                has_data = True
                if val < local_min: local_min = val
                if val > local_max: local_max = val
    if has_data:
        h_sub.SetMinimum(local_min)
        h_sub.SetMaximum(local_max)
    else:
        h_sub.SetMinimum(0.0)
        h_sub.SetMaximum(1.0)
    h_sub.Draw("COLZ")
    c.SaveAs(f"{name}.pdf")
    c.Close()
    print(f" Subtracted plot saved for {target_label}.")
    return {"stat": h_sub_stat, "sys": h_sub_sys}

def generate_pd_subtracted_plot(hists_ld2, hists_lh2, hists_flask, output_file):
    print("\nGenerating Subtracted Plot for LD2 pd formula...")
    output_file.cd()

    h_ld2_stat = hists_ld2["Y_corrected_stat"]
    h_ld2_sys = hists_ld2["Y_corrected_sys"]
    h_lh2_stat = hists_lh2["Y_corrected_stat"]
    h_lh2_sys = hists_lh2["Y_corrected_sys"]
    h_flask_stat = hists_flask["Y_corrected_stat"]
    h_flask_sys = hists_flask["Y_corrected_sys"]

    c_lh2 = THD_THH_RATIO * (PROTONS_ON_TARGET_LD2 / PROTONS_ON_TARGET_LH2)
    c_flask_sub = FLASK_NORM_LD2 - (c_lh2 * FLASK_NORM_LH2)

    name = "Y_corrected_Subtracted_LD2"
    title = "Corrected Yield (LD2 - LH2 - Flask)"
    h_pd = ROOT.TH2D(name, f"{title};Mass [GeV];x_{{F}}",
    len(mass_bins_np)-1, mass_bins_np,
    len(xf_bins_np)-1, xf_bins_np)
    h_pd.Sumw2()
    h_pd.SetStats(0)
    h_pd.GetXaxis().CenterTitle()
    h_pd.GetYaxis().CenterTitle()
    h_pd.GetXaxis().SetTitleOffset(1.2)
    h_pd.GetYaxis().SetTitleOffset(1.2)

    h_pd_stat = h_pd.Clone(f"{name}_stat")
    h_pd_sys = h_pd.Clone(f"{name}_sys")

    for i_m in range(len(mass_bins_np) - 1):
        m_center = (mass_bins_np[i_m] + mass_bins_np[i_m+1]) / 2.0
        bin_x = i_m + 1
        for i_x in range(len(xf_bins_np) - 1):
            x_center = (xf_bins_np[i_x] + xf_bins_np[i_x+1]) / 2.0
            bin_y = i_x + 1
            y_ld2 = h_ld2_stat.GetBinContent(bin_x, bin_y)
            y_lh2 = h_lh2_stat.GetBinContent(bin_x, bin_y)
            y_flask = h_flask_stat.GetBinContent(bin_x, bin_y)

            e_ld2_stat = h_ld2_stat.GetBinError(bin_x, bin_y)
            e_lh2_stat = h_lh2_stat.GetBinError(bin_x, bin_y)
            e_flask_stat = h_flask_stat.GetBinError(bin_x, bin_y)
            e_ld2_sys = h_ld2_sys.GetBinError(bin_x, bin_y)
            e_lh2_sys = h_lh2_sys.GetBinError(bin_x, bin_y)
            e_flask_sys = h_flask_sys.GetBinError(bin_x, bin_y)
            val_pd = y_ld2 - (c_lh2 * y_lh2) - (c_flask_sub * y_flask)
            err_pd_stat = np.sqrt(e_ld2_stat**2 + (c_lh2 * e_lh2_stat)**2 + (c_flask_sub * e_flask_stat)**2)
            err_pd_sys = np.sqrt(e_ld2_sys**2 + (c_lh2 * e_lh2_sys)**2 + (c_flask_sub * e_flask_sys)**2)
            err_pd_total = np.sqrt(err_pd_stat**2 + err_pd_sys**2)
            h_pd.SetBinContent(bin_x, bin_y, val_pd)
            h_pd.SetBinError(bin_x, bin_y, err_pd_total)
            add_latex_to_bin(h_pd, m_center, x_center, val_pd, err_pd_total)

            h_pd_stat.SetBinContent(bin_x, bin_y, val_pd)
            h_pd_stat.SetBinError(bin_x, bin_y, err_pd_stat)

            h_pd_sys.SetBinContent(bin_x, bin_y, val_pd)
            h_pd_sys.SetBinError(bin_x, bin_y, err_pd_sys)

    c = ROOT.TCanvas(f"c_{name}", f"c_{name}", 1200, 900)
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.SetTickx(1)
    c.SetTicky(1)
    local_min = sys.float_info.max
    local_max = -sys.float_info.max
    has_data = False
    for ix in range(1, h_pd.GetNbinsX() + 1):
        for iy in range(1, h_pd.GetNbinsY() + 1):
            val = h_pd.GetBinContent(ix, iy)
            if val != 0.0:
                has_data = True
                if val < local_min: local_min = val
                if val > local_max: local_max = val
    if has_data:
        h_pd.SetMinimum(local_min)
        h_pd.SetMaximum(local_max)
    else:
        h_pd.SetMinimum(0.0)
        h_pd.SetMaximum(1.0)
    h_pd.Draw("COLZ")
    c.SaveAs(f"{name}.pdf")
    c.Close()
    print(" Subtracted plot saved for LD2 pd equation.")

    return {"stat": h_pd_stat, "sys": h_pd_sys}


# ==========================================
# 7. Cross Section Calculation & Plotting
# ==========================================
def calculate_and_plot_cross_section(h_sub_dict, target_label, global_constant, output_file):
    print(f"\nCalculating and Plotting Double Differential Cross-Sections for {target_label}...")
    acc_path = "acceptance_mass_xF.root"
    psip_path = "All_PsiP_Contaminations.root"
    # Dynamically assign theoretical files based on the target
    if target_label == "LD2":
        theory_ct18_path = "CT18_xFnew_d_1sigma.root"
        theory_nnpdf_path = "NNPDF40_xFnew_d.root"
    else:
        theory_ct18_path = "CT18_xFnew_p_1sigma.root"
        theory_nnpdf_path = "NNPDF40_xFnew_p.root"
    if not os.path.exists(acc_path):
        raise FileNotFoundError(f"CRITICAL ERROR: Acceptance file '{acc_path}' not found!")

    # Open Files
    try:
        acc_file = ROOT.TFile.Open(acc_path)
        if not acc_file or acc_file.IsZombie():
            raise OSError(f"Acceptance file corrupted.")
        f_ct18 = None
        if os.path.exists(theory_ct18_path):
            f_ct18 = ROOT.TFile.Open(theory_ct18_path)
        else:
            print(f"Warning: Theoretical file '{theory_ct18_path}' not found.")
        f_nnpdf = None
        if os.path.exists(theory_nnpdf_path):
            f_nnpdf = ROOT.TFile.Open(theory_nnpdf_path)
        else:
            print(f"Warning: Theoretical file '{theory_nnpdf_path}' not found.")
        f_psip = None
        if os.path.exists(psip_path):
            f_psip = ROOT.TFile.Open(psip_path)
        else:
            print(f"Warning: PsiP Contamination file '{psip_path}' not found. Systematic will remain 0.")
    except Exception as e:
        print(f"Error opening files: {e}")
        sys.exit(1)

    output_file.cd()
    h_sub_stat = h_sub_dict["stat"]
    h_sub_sys = h_sub_dict["sys"]

    n_xf_bins = len(xf_bins_np) - 1
    n_mass_bins = len(mass_bins_np) - 1
    # Initialize PsiP Contamination Table Content
    latex_psip_table_content = r"""\begin{longtable}{|c|c|c|c|c|}
\caption{$\psi'$ Contamination Table for %s} \label{tab:psip_contamination_%s} \\
\hline
\textbf{xF bin} & \textbf{Mass bin (GeV)} & \textbf{$\psi'$ contamination} & \textbf{$\sigma \pm \delta \sigma^{stat.} \pm \delta \sigma^{syst.}$} & \textbf{$\delta\sigma_{\psi'}^{\rm syst.}$ (nb GeV$^2$)} \\
\hline
\endfirsthead

\multicolumn{5}{c}%%
{{\bfseries \tablename\ \thetable{} -- continued from previous page}} \\
\hline
\textbf{xF bin} & \textbf{Mass bin (GeV)} & \textbf{$\psi'$ contamination} & \textbf{$\sigma \pm \delta \sigma^{stat.} \pm \delta \sigma^{syst.}$} & \textbf{$\delta\sigma_{\psi'}^{\rm syst.}$ (nb GeV$^2$)} \\
\hline
\endhead

\hline \multicolumn{5}{|r|}{{Continued on next page}} \\ \hline
\endfoot

\hline
\endlastfoot
""" % (target_label, target_label)


    for i_x in range(n_xf_bins):
        xf_bin_index = i_x
        root_bin_y = i_x + 1
        xf_min = xf_bins_np[i_x]
        xf_max = xf_bins_np[i_x+1]
        
        # Determine the key for Shivangi's data dict based on xF
        xf_key = f"{xf_min:.1f}-{xf_max:.2f}".replace(".0-", ".00-")
        if xf_min == 0.0: xf_key = "0.0-0.05"
        elif xf_min == 0.05: xf_key = "0.05-0.1"
        elif xf_min == 0.1: xf_key = "0.1-0.15"
        elif xf_min == 0.15: xf_key = "0.15-0.2"
        elif xf_min == 0.2: xf_key = "0.2-0.25"
        elif xf_min == 0.25: xf_key = "0.25-0.3"
        elif xf_min == 0.3: xf_key = "0.3-0.35"
        elif xf_min == 0.35: xf_key = "0.35-0.4"
        elif xf_min == 0.4: xf_key = "0.4-0.45"
        elif xf_min == 0.45: xf_key = "0.45-0.5"
        elif xf_min == 0.5: xf_key = "0.5-0.55"
        elif xf_min == 0.55: xf_key = "0.55-0.6"
        elif xf_min == 0.6: xf_key = "0.6-0.65"
        elif xf_min == 0.65: xf_key = "0.65-0.7"
        elif xf_min == 0.7: xf_key = "0.7-0.75"
        elif xf_min == 0.75: xf_key = "0.75-0.8"

        # Determine the Acceptance Histogram (Fallback to LH2 if target-specific one doesn't exist)
        acc_hist_name = f"h_ratio_{target_label}_xF_bin{xf_bin_index}"
        h_acc = acc_file.Get(acc_hist_name)
        if not h_acc:
            fallback_name = f"h_ratio_LH2_xF_bin{xf_bin_index}"
            h_acc = acc_file.Get(fallback_name)
            if h_acc:
                print(f" Note: Using fallback acceptance histogram '{fallback_name}'.")

        if not h_acc:
            print(f"Warning: Acceptance histogram for xF bin {xf_bin_index} not found. Skipping bin.")
            continue
        # Load the PsiP/DY ratio histogram for this xF bin
        h_ratio_psip = None
        if f_psip:
            h_ratio_psip = f_psip.Get(f"hRatio_PsiP_DY_xF_{xf_bin_index}")
        print(f"\n--- {target_label} xF Bin {xf_bin_index}: {xf_min:.2f} < xF < {xf_max:.2f} ---")

        # --- Data Graph (User's Script, Red) ---
        g_xsec = ROOT.TGraphErrors()
        g_xsec.SetName(f"g_xsec_{target_label}_xF_{xf_bin_index}")
        g_xsec.SetTitle(f";Mass [GeV];M^{{3}} d^{{2}}\\sigma / dM dx_{{F}} [nb GeV^{{2}}]")
        
        # --- Systematic Error Band Graph (User's Script, Red) ---
        g_sys = ROOT.TGraphErrors()
        g_sys.SetName(f"g_sys_{target_label}_xF_{xf_bin_index}")
        
        # --- Shivangi's Data Graph (Blue) ---
        g_xsec_shivangi = ROOT.TGraphErrors()
        g_xsec_shivangi.SetName(f"g_xsec_shivangi_{target_label}_xF_{xf_bin_index}")
        
        # --- Shivangi's Systematic Error Band Graph (Blue) ---
        g_sys_shivangi = ROOT.TGraphErrors()
        g_sys_shivangi.SetName(f"g_sys_shivangi_{target_label}_xF_{xf_bin_index}")

        point_idx = 0
        point_idx_shivangi = 0

        # Variables to track dynamic min/max for the Y-axis
        y_min_data = sys.float_info.max
        y_max_data = -sys.float_info.max
        
        for i_m in range(n_mass_bins):
            root_bin_x = i_m + 1
            mass_min = mass_bins_np[i_m]
            mass_max = mass_bins_np[i_m+1]
            geometric_center = (mass_min + mass_max) / 2.0
            mass_width = mass_max - mass_min
            
            # --- 1. User's Calculations (Red) ---
            Y_sub = h_sub_stat.GetBinContent(root_bin_x, root_bin_y)
            Y_sub_stat_err = h_sub_stat.GetBinError(root_bin_x, root_bin_y)
            Y_sub_sys_err = h_sub_sys.GetBinError(root_bin_x, root_bin_y)
            
            if Y_sub > 0:
                psip_ratio = 0.0
                psip_ratio_err = 0.0
                if h_ratio_psip:
                    ratio_bin = h_ratio_psip.FindBin(geometric_center)
                    psip_ratio = h_ratio_psip.GetBinContent(ratio_bin)
                    psip_ratio_err = h_ratio_psip.GetBinError(ratio_bin)

                if psip_ratio > 1.0:
                    psip_ratio = 0.0
                    psip_ratio_err = 0.0

                acc_bin = h_acc.FindBin(geometric_center)
                acceptance = h_acc.GetBinContent(acc_bin)
                acceptance_err = h_acc.GetBinError(acc_bin)
                print(f" Mass Centroid {geometric_center:.3f} GeV: Acceptance = {acceptance:.5f} +/- {acceptance_err:.5f}")
                
                if acceptance > 0:
                    numerator = global_constant * Y_sub
                    denominator = mass_width * acceptance
                    xsec = numerator / denominator
                    
                    stat_unc = (Y_sub_stat_err / Y_sub) * xsec
                    sys_acc = (acceptance_err / acceptance) * xsec
                    sys_tot_yield = (Y_sub_sys_err / Y_sub) * xsec
                    sys_psip_cont = psip_ratio * xsec if Y_sub > 0 else 0.0
                    total_systematic_unc = np.sqrt(sys_acc**2 + sys_tot_yield**2 + sys_psip_cont**2)
                    
                    scaling_factor = geometric_center**3
                    scaled_xsec = xsec * scaling_factor
                    scaled_stat_err = stat_unc * scaling_factor
                    scaled_sys_err = total_systematic_unc * scaling_factor
                    scaled_sys_psip = sys_psip_cont * scaling_factor
                    
                    if psip_ratio > 0.0:
                        s_xf = f"[{xf_min:.2f}, {xf_max:.2f})"
                        s_mass = f"[{mass_min:.2f}, {mass_max:.2f})"
                        s_ratio = f"{psip_ratio:.4f} $\\pm$ {psip_ratio_err:.4f}"
                        s_sigma_psip_col = f"{scaled_xsec:.4f} $\\pm$ {scaled_stat_err:.4f} $\\pm$ {scaled_sys_err:.4f}"
                        s_sigma_psip = f"{scaled_sys_psip:.4f}"
                        latex_psip_table_content += f"{s_xf} & {s_mass} & {s_ratio} & {s_sigma_psip_col} & {s_sigma_psip} \\\\ \n"
                        latex_psip_table_content += r"\hline" + "\n"
                        
                    # Plot Red Data
                    g_sys.SetPoint(point_idx, geometric_center, scaled_xsec)
                    g_sys.SetPointError(point_idx, 0.0, scaled_sys_err) # <-- CHANGED TO 0.0 X ERROR
                    
                    g_xsec.SetPoint(point_idx, geometric_center, scaled_xsec)
                    g_xsec.SetPointError(point_idx, 0.0, scaled_stat_err)
                    
                    max_err_for_range = max(scaled_stat_err, scaled_sys_err)
                    y_high = scaled_xsec + max_err_for_range
                    y_low = scaled_xsec - max_err_for_range
                    if y_low <= 0: y_low = scaled_xsec * 0.5
                    if y_high > y_max_data: y_max_data = y_high
                    if y_low < y_min_data: y_min_data = y_low
                    
                    point_idx += 1

            # --- 2. Shivangi's Hardcoded Data (Blue) ---
            if target_label == "LD2" and xf_key in SHIVANGI_DATA:
                # Find matching mass bin in dictionary
                for row in SHIVANGI_DATA[xf_key]:
                    # row format: [mass_min, mass_max, xsec, stat_err, sys_err_low, sys_err_high]
                    # Give tiny tolerance for float matching
                    if abs(row[0] - mass_min) < 0.01 and abs(row[1] - mass_max) < 0.01:
                        shiv_xsec = row[2]
                        if shiv_xsec != 0.0:
                            # Divide by 10 since the table header says (x10^-1)
                            shiv_xsec = shiv_xsec / 10.0
                            shiv_stat = row[3] / 10.0
                            # Averaging asymmetric sys errs for ROOT TGraphErrors symmetric error band
                            shiv_sys = ((row[4] / 10.0) + (row[5] / 10.0)) / 2.0 
                            
                            # Use Geometric center just like the new script logic
                            g_sys_shivangi.SetPoint(point_idx_shivangi, geometric_center, shiv_xsec)
                            g_sys_shivangi.SetPointError(point_idx_shivangi, 0.0, shiv_sys) # <-- CHANGED TO 0.0 X ERROR
                            
                            g_xsec_shivangi.SetPoint(point_idx_shivangi, geometric_center, shiv_xsec)
                            g_xsec_shivangi.SetPointError(point_idx_shivangi, 0.0, shiv_stat)
                            
                            max_err_for_range = max(shiv_stat, shiv_sys)
                            y_high = shiv_xsec + max_err_for_range
                            y_low = shiv_xsec - max_err_for_range
                            if y_low <= 0: y_low = shiv_xsec * 0.5
                            if y_high > y_max_data: y_max_data = y_high
                            if y_low < y_min_data: y_min_data = y_low
                            
                            point_idx_shivangi += 1
                        break


        # --- Plotting ---
        if g_xsec.GetN() > 0 or g_xsec_shivangi.GetN() > 0:
            c_xsec = ROOT.TCanvas(f"c_xsec_{target_label}_{xf_bin_index}", f"Cross Section {target_label} xF {xf_bin_index}", 800, 600)
            c_xsec.SetLogy()
            c_xsec.SetTickx(1)
            c_xsec.SetTicky(1)
            mg = ROOT.TMultiGraph()
            mg.SetTitle(f";Mass [GeV];M^{{3}} d^{{2}}\\sigma / dM dx_{{F}} [nb GeV^{{2}}]")
            
            if y_min_data < y_max_data:
                mg.SetMinimum(y_min_data * 0.2)
                mg.SetMaximum(y_max_data * 5.0)
            else:
                mg.SetMinimum(1e-3)
                mg.SetMaximum(3.0)
                
            leg = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)
            leg.SetBorderSize(0)
            
            gr_name = f"gr_xFbin{xf_bin_index}"
            
            # Theory Lines
            if f_ct18:
                g_ct18 = f_ct18.Get(gr_name)
                if g_ct18:
                    g_ct18.SetLineColor(ROOT.kGreen + 2)
                    g_ct18.SetFillColorAlpha(ROOT.kGreen - 5, 0.5)
                    g_ct18.SetFillStyle(3002)
                    #mg.Add(g_ct18, "L3")
                    #leg.AddEntry(g_ct18, "CT18 NLO", "lf")
            if f_nnpdf:
                g_nnpdf = f_nnpdf.Get(gr_name)
                if g_nnpdf:
                    g_nnpdf.SetLineColor(ROOT.kGray + 2)
                    g_nnpdf.SetFillColorAlpha(ROOT.kGray + 1, 0.5)
                    g_nnpdf.SetFillStyle(3002)
                    #mg.Add(g_nnpdf, "L3")
                    #leg.AddEntry(g_nnpdf, "NNPDF4.0 NLO", "lf")

            # Shivangi's Blue Data
            if target_label == "LD2" and g_xsec_shivangi.GetN() > 0:
                # --- CHANGED: Thick inner systematic error bar ---
                g_sys_shivangi.SetMarkerSize(0)
                g_sys_shivangi.SetLineColor(ROOT.kBlue)
                g_sys_shivangi.SetLineWidth(4) # Thick Line Width
                mg.Add(g_sys_shivangi, "Z") # "Z" prevents horizontal caps from drawing
                
                # Standard statistical error bar and marker
                g_xsec_shivangi.SetMarkerStyle(21) # Square marker
                g_xsec_shivangi.SetMarkerColor(ROOT.kBlue)
                g_xsec_shivangi.SetLineColor(ROOT.kBlue)
                g_xsec_shivangi.SetLineWidth(1)
                mg.Add(g_xsec_shivangi, "P")
                
                leg.AddEntry(g_xsec_shivangi, f"Shivangi's Data ({target_label})", "lep")
                leg.AddEntry(g_sys_shivangi, "Shivangi Syst.", "l") # "l" matches the line for legend

            # User's Red Data
            if g_xsec.GetN() > 0:
                # --- CHANGED: Thick inner systematic error bar ---
                g_sys.SetMarkerSize(0)
                g_sys.SetLineColor(ROOT.kRed)
                g_sys.SetLineWidth(4) # Thick Line Width
                mg.Add(g_sys, "Z") # "Z" prevents horizontal caps from drawing
                
                # Standard statistical error bar and marker
                g_xsec.SetMarkerStyle(20) # Circle marker
                g_xsec.SetMarkerColor(ROOT.kRed)
                g_xsec.SetLineColor(ROOT.kRed)
                g_xsec.SetLineWidth(1)
                mg.Add(g_xsec, "P")
                
                leg.AddEntry(g_xsec, f"Current Data ({target_label})", "lep")
                leg.AddEntry(g_sys, "Current Syst.", "l") # "l" matches the line for legend

            mg.Draw("A")
            mg.GetXaxis().CenterTitle()
            mg.GetYaxis().CenterTitle()
            mg.GetXaxis().SetLimits(4.2, 8.7)
            leg.Draw()

            # Annotations
            target_prefix = "pp" if target_label == "LH2" else "pd" if target_label == "LD2" else target_label
            internal_title_text = f"{target_prefix} cross-section within {xf_min:.2f} #leq x_{{F}} < {xf_max:.2f}"
            internal_title = ROOT.TLatex()
            internal_title.SetNDC(True)
            internal_title.SetTextFont(42)
            internal_title.SetTextSize(0.04)
            internal_title.SetTextAlign(13)
            internal_title.DrawLatex(0.14, 0.86, internal_title_text)

            prelim = ROOT.TLatex()
            prelim.SetNDC(True)
            prelim.SetTextColor(ROOT.kBlue)
            prelim.SetTextAlign(33)
            prelim.SetTextSize(0.05)
            #prelim.DrawLatex(0.82, 0.6, "Preliminary")

            lumi_note = ROOT.TLatex()
            lumi_note.SetNDC(False)
            lumi_note.SetTextColor(ROOT.kBlack)
            lumi_note.SetTextAlign(13)
            lumi_note.SetTextSize(0.025)
            if y_min_data < y_max_data:
                y_min_actual = y_min_data * 0.2
            else:
                y_min_actual = 1e-3
            if y_min_actual <= 0:
                y_min_actual = 1e-3
            lumi_y_pos = y_min_actual * 1.5
            #lumi_note.DrawLatex(4.3, lumi_y_pos, "10% global uncertainty due to the integrated luminosity is not included in the error bands")
            
            # Save Output
            plot_name = f"Comparison_{target_label}_xF_{xf_min:.2f}_{xf_max:.2f}.pdf"
            c_xsec.SaveAs(plot_name)
            g_xsec.Write(f"g_xsec_{target_label}_{xf_bin_index}")
            g_sys.Write(f"g_sys_{target_label}_{xf_bin_index}")
            
            if target_label == "LD2":
                g_xsec_shivangi.Write(f"g_xsec_shivangi_{target_label}_{xf_bin_index}")
                g_sys_shivangi.Write(f"g_sys_shivangi_{target_label}_{xf_bin_index}")
                
            c_xsec.Close()
            print(f" Saved plot: {plot_name}")
        else:
            print(f" No valid points for xF bin {xf_bin_index}.")

    latex_psip_table_content += r"\end{longtable}" + "\n"
    psip_table_filename = f"Table_PsiP_Contamination_{target_label}.tex"
    with open(psip_table_filename, "w") as f:
        f.write(latex_psip_table_content)
    print(f"\nSaved PsiP contamination LaTeX table to {psip_table_filename}")

    acc_file.Close()
    if f_ct18: f_ct18.Close()
    if f_nnpdf: f_nnpdf.Close()
    if f_psip: f_psip.Close()


# ==========================================
# Main Execution
# ==========================================
def main():
    output_filename = "Processed_Kinematic_Hists.root"
    print(f"Creating ROOT output file: {output_filename}")
    out_file = ROOT.TFile(output_filename, "RECREATE")

    hists_lh2 = process_file_data("../../../HodoEfficiency/RS67/LH2/merged_RS67_3089_LH2_recoeff_hodoeff.root", "LH2", out_file)
    hists_flask = process_file_data("../../../HodoEfficiency/RS67/EmptyFlask/merged_RS67_3089_Flask_recoeff_hodoeff.root", "Flask", out_file)
    hists_ld2 = process_file_data("../../../HodoEfficiency/RS67/LD2/merged_RS67_3089_LD2_recoeff_hodoeff.root", "LD2", out_file)

    # 1. Process standard LH2 cross-sections
    if hists_lh2 and hists_flask:
        h_sub_dict_lh2 = generate_subtracted_plot(hists_lh2, hists_flask, FLASK_NORM_LH2, "LH2", out_file)
        if h_sub_dict_lh2:
            calculate_and_plot_cross_section(h_sub_dict_lh2, "LH2", GLOBAL_CONSTANT_LH2, out_file)
        else:
            print("\nSkipping LH2 subtraction plotting because LH2 or Flask histograms were not generated.")

    # 2. Process LD2 using the user's specific pd extraction equation
    if hists_ld2 and hists_lh2 and hists_flask:
        h_sub_dict_pd = generate_pd_subtracted_plot(hists_ld2, hists_lh2, hists_flask, out_file)
        if h_sub_dict_pd:
            calculate_and_plot_cross_section(h_sub_dict_pd, "LD2", GLOBAL_CONSTANT_LD2, out_file)
        else:
            print("\nSkipping LD2 pd extraction plotting because LD2, LH2, or Flask histograms are missing.")

    # ---> FLUSH EVERYTHING TO DISK <---
    out_file.Write()

    out_file.Close()
    print("\nAll histograms, tables, and cross-section plots generated successfully.")

if __name__ == "__main__":
    main()