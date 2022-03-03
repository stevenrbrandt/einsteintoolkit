import numpy as np

eta = 0.24000000000000002

Domain = [0., 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02, 0.022, 0.024, 0.026, 0.028, 0.03,
          0.032, 0.034, 0.036, 0.038, 0.04, 0.042, 0.044, 0.046, 0.048,0.05, 0.052, 0.054, 0.056, 0.058, 0.06, 0.062,
          0.064, 0.066, 0.068, 0.07, 0.072, 0.074, 0.076, 0.078, 0.08, 0.082, 0.084, 0.086, 0.088, 0.09, 0.092, 0.094,
          0.096, 0.098, 0.1, 0.102, 0.104, 0.106, 0.108, 0.11, 0.112, 0.114, 0.116, 0.118, 0.12, 0.122, 0.124, 0.126,
          0.128, 0.13, 0.132, 0.134, 0.136, 0.138, 0.14, 0.142, 0.144, 0.146, 0.148, 0.15, 0.152, 0.154, 0.156, 0.158,
          0.16, 0.162, 0.164, 0.166, 0.168, 0.17, 0.172, 0.174, 0.176, 0.178, 0.18, 0.182, 0.184, 0.186, 0.188, 0.19,
          0.192, 0.194, 0.196, 0.198, 0.2, 0.202, 0.204, 0.206, 0.208, 0.21, 0.212, 0.214, 0.216, 0.218, 0.22, 0.222,
          0.224, 0.226, 0.228, 0.23, 0.232, 0.234, 0.236, 0.238, 0.24, 0.242, 0.244, 0.246, 0.248, 0.25]

a1Range = [-12.9867168, -12.3922341, -11.8826575, -11.5044206, -11.2010798, -10.9291101, -10.6867579, -10.4603116,
           -10.2705936, -10.1096886, -9.94972042, -9.7993101, -9.67017397, -9.54568872, -9.42387643, -9.30374853,
           -9.19233175, -9.08443818, -8.97983162, -8.87983878, -8.78334655, -8.68874044, -8.5958085, -8.50609937,
           -8.42033859, -8.3368237, -8.25318119, -8.17059891, -8.0893624, -8.00915281, -7.92946769, -7.85089647,
           -7.77416021, -7.69874159, -7.6230207, -7.54778069, -7.47506574, -7.40420786, -7.33414965, -7.26472564,
           -7.19583939, -7.12778746, -7.06152986, -6.9966294, -6.93284618, -6.86298297, -6.78959864, -6.73156524,
           -6.67432199, -6.61818773, -6.56085433, -6.50488932, -6.45206856, -6.40421195, -6.35784803, -6.3138088,
           -6.26993274, -6.22465506, -6.18289146, -6.14263195, -6.10558528, -6.07015742, -6.03568352, -6.00473068,
           -5.97536374, -5.9488801, -5.92520739, -5.89986979, -5.87411946, -5.85174778, -5.83384542, -5.82778792,
           -5.82798078, -5.8256649, -5.82049035, -5.81770024, -5.81990472, -5.82666469, -5.83577969, -5.84678194,
           -5.8593275, -5.87380902, -5.88835704, -5.90227283, -5.91692664, -5.93375623, -5.95112447, -5.97214313,
           -5.99595282, -6.01998938, -6.0451582, -6.07037506, -6.09500227, -6.12087467, -6.1457118, -6.16991637,
           -6.19524357, -6.21896485, -6.24302387, -6.26778731, -6.29335024, -6.32209441, -6.35041468, -6.37940448,
           -6.40967178, -6.44669514, -6.4822263, -6.51034108, -6.53379526, -6.5510243, -6.5791139, -6.61538596,
           -6.65893348, -6.70365878, -6.7474139, -6.79400218, -6.8495302, -6.91475643, -6.98767056, -7.06220586,
           -7.16459696, -7.29175945, -7.43577701, -7.60975023, -7.86628294, -8.249368]

a2Range = [103.047993, 97.488689, 92.7138991, 89.1259949, 86.2216369, 83.5996035, 81.2588265, 79.0565575, 77.1901259,
           75.5962053, 74.0126613, 72.5193589, 71.2189287, 69.960102, 68.7278401, 67.5137535, 66.3813558, 65.2852124,
           64.2251046, 63.2086791, 62.2258845, 61.2625581, 60.3162743, 59.4001895, 58.5221931, 57.666497, 56.8109787,
           55.9670566, 55.1376064, 54.3195185, 53.5077914, 52.7084763, 51.9289985, 51.1641196, 50.3976552, 49.6369195,
           48.9014988, 48.184488, 47.4747413, 46.7739132, 46.0794666, 45.3892075, 44.7162594, 44.0576325, 43.4108107,
           42.7414687, 42.0328661, 41.4379676, 40.8519633, 40.2763931, 39.6864942, 39.1075193, 38.5558168, 38.0487152,
           37.5534784, 37.0855881, 36.6187593, 36.1354218, 35.6891039, 35.2559943, 34.8506441, 34.4544071, 34.0606759,
           33.7102444, 33.3728535, 33.0425553, 32.7368071, 32.424474, 32.1183409, 31.8446322, 31.6055987, 31.4663075,
           31.3800944, 31.2792375, 31.1629634, 31.0686905, 31.0129002, 30.9902402, 30.9910664, 31.0100087, 31.0447735,
           31.096186, 31.1545669, 31.2147854, 31.2846909, 31.3747464, 31.4774613, 31.6133923, 31.7749697, 31.9441796,
           32.1267434, 32.3145467, 32.5026478, 32.7052324, 32.9044074, 33.103595, 33.3182241, 33.5238512, 33.7360635,
           33.9573682, 34.1872127, 34.4454878, 34.7018803, 34.9661177, 35.2447179, 35.576539, 35.8979802, 36.1672599,
           36.4055123, 36.592871, 36.8680094, 37.2087221, 37.6092315, 38.0164539, 38.4141907, 38.8351392, 39.3247457,
           39.8822802, 40.4926366, 41.1147783, 41.9315173, 42.9304965, 44.0528309, 45.3955045, 47.3274085, 50.0808317]

a3Range = [-112.067972, -105.543106, -99.9455353, -95.7264397, -92.3043493, -89.2092924, -86.4492519, -83.8470286,
           -81.6361633, -79.7483937, -77.874168, -76.1068934, -74.5632033, -73.0674934, -71.6039815, -70.1633203,
           -68.8189591, -67.519153, -66.2647734, -65.0620239, -63.8992811, -62.7604276, -61.6423525, -60.5597465,
           -59.5224349, -58.5119795, -57.5024912, -56.5074902, -55.53044, -54.5676408, -53.6130924, -52.6741387,
           -51.759812, -50.8637645, -49.9665188, -49.0765995, -48.2171201, -47.379661, -46.5505356, -45.7333476,
           -44.9242629, -44.118089, -43.3322908, -42.5639012, -41.8098449, -41.043258, -40.2287372, -39.5338046,
           -38.8499135, -38.1781038, -37.4888781, -36.8113298, -36.1641156, -35.5671876, -34.9832785, -34.4326665,
           -33.8826254, -33.3121583, -32.7856761, -32.2743216, -31.7933477, -31.3201055, -30.8470575, -30.4267354,
           -30.0203142, -29.6142881, -29.2352514, -28.8518367, -28.478952, -28.1429549, -27.8435553, -27.6538903,
           -27.5225161, -27.3767334, -27.2165867, -27.0811087, -26.986459, -26.9258851, -26.8924026, -26.8792128,
           -26.8840645, -26.9066881, -26.939216, -26.976347, -27.0248923, -27.0956124, -27.1838684, -27.3095255,
           -27.4637058, -27.6286706, -27.8094514, -27.9976633, -28.1882304, -28.3957603, -28.6018695, -28.8099936,
           -29.0365609, -29.2555204, -29.4831745, -29.7220365, -29.9712378, -30.2520714, -30.5319354, -30.8213449,
           -31.1276025, -31.4916256, -31.8453217, -32.1448723, -32.4128347, -32.6260766, -32.9352624, -33.3159744,
           -33.762117, -34.2149428, -34.6570779, -35.1244967, -35.6653064, -36.2763704, -36.9410027, -37.61746,
           -38.493775, -39.5607747, -40.7566336, -42.1829659, -44.2164331, -47.0600574]

b1Range = [-0.0409516519, -0.109149589, -0.18280716, -0.243877638, -0.29427273, -0.337607837, -0.376514835, -0.408102051,
           -0.430340476, -0.447661929, -0.455696226, -0.455620099, -0.451829661, -0.448183565, -0.444485773,
           -0.441429709, -0.43796189, -0.435013771, -0.43317134, -0.431090718, -0.428804204, -0.426754544, -0.4248826,
           -0.422471012, -0.419498947, -0.416391867, -0.413725074, -0.411268992, -0.408996729, -0.407008137,
           -0.405444696, -0.404208778, -0.403206029, -0.402514117, -0.402371401, -0.402655408, -0.403094271,
           -0.403632835, -0.404112538, -0.405565386, -0.407338113, -0.407607441, -0.408032297, -0.409123955,
           -0.410936785, -0.415814407, -0.417806894, -0.420418936, -0.423769675, -0.427415314, -0.431003883,
           -0.434312853, -0.437542567, -0.440912596, -0.445155727, -0.450792563, -0.456134467, -0.461308333,
           -0.467704275, -0.475422557, -0.482471914, -0.489063765, -0.496174659, -0.503981188, -0.513338867,
           -0.525302557, -0.536971352, -0.546967233, -0.555720559, -0.565412016, -0.576313878, -0.586599947,
           -0.596427341, -0.604860719, -0.611970131, -0.619519476, -0.628821459, -0.639872749, -0.65090668,
           -0.661411175, -0.671571431, -0.682134421, -0.692764125, -0.703573421, -0.715112958, -0.727369423,
           -0.738621336, -0.749657644, -0.761066387, -0.77191596, -0.782873523, -0.793717364, -0.804225625, -0.815026186,
           -0.825216179, -0.835238792, -0.845945128, -0.855787288, -0.865303002, -0.874495405, -0.883003789,
           -0.891163975, -0.898663309, -0.905882454, -0.91320713, -0.919704848, -0.926061344, -0.933262667, -0.940928324,
           -0.947831369, -0.953153169, -0.95693697, -0.959330857, -0.960185331, -0.959843982, -0.958686944, -0.956619418,
           -0.952843922, -0.946447305, -0.939771599, -0.928224314, -0.911651306, -0.891719944, -0.867854188, -0.82502669,
           -0.73315571]

b2Range = [-0.00166722842, 0.222250337, 0.448252485, 0.63206119, 0.781799547, 0.910224883, 1.02375301, 1.1157851,
           1.17923771, 1.22659824, 1.24771432, 1.24548805, 1.2320439, 1.2192957, 1.2064752, 1.19588522, 1.18435194,
           1.17453472, 1.16812456, 1.16113472, 1.15387353, 1.14747797, 1.14166252, 1.13468889, 1.12673472, 1.11873428,
           1.11191795, 1.10593347, 1.10080563, 1.09667658, 1.09370713, 1.09190117, 1.09129942, 1.09194833, 1.0939851,
           1.09725207, 1.1015065, 1.10645069, 1.11139159, 1.11928716, 1.12839242, 1.13375028, 1.14006609, 1.14868099,
           1.15958421, 1.18126715, 1.19236439, 1.20564444, 1.22112152, 1.23765157, 1.25462197, 1.27123058, 1.28793848,
           1.30524443, 1.32551845, 1.34960096, 1.37341634, 1.39766424, 1.4249208, 1.45597978, 1.48516877, 1.51368441,
           1.5446322, 1.5767945, 1.61477856, 1.66514258, 1.71502079, 1.75728168, 1.79328757, 1.83222422, 1.87523717,
           1.91514764, 1.95277569, 1.98471429, 2.01131332, 2.03946121, 2.07428419, 2.11572946, 2.15704477, 2.19631812,
           2.23416378, 2.27326935, 2.3122869, 2.35158306, 2.39309063, 2.43665605, 2.47607718, 2.51394563, 2.55248145,
           2.58860308, 2.62460385, 2.65980703, 2.6935071, 2.72765112, 2.75939003, 2.79018574, 2.82260365, 2.85186527,
           2.87952241, 2.90554772, 2.92874698, 2.94977991, 2.96834961, 2.98551006, 3.0023529, 3.01499219, 3.02720026,
           3.04304093, 3.06075921, 3.07707208, 3.08528399, 3.08627843, 3.08073147, 3.06990182, 3.05552041, 3.03783353,
           3.0153212, 2.98607118, 2.94810984, 2.90941344, 2.85250477, 2.77660316, 2.68817051, 2.58115094, 2.4104722,
           2.11713253]

a1 = np.interp(eta, Domain, a1Range)
a2 = np.interp(eta, Domain, a2Range)
a3 = np.interp(eta, Domain, a3Range)
b1 = -np.interp(eta, Domain, b1Range)
b2 = -np.interp(eta, Domain, b2Range)
