import FWCore.ParameterSet.Config as cms

from L1Trigger.L1TCalorimeter.caloParams_cfi import caloParamsSource
import L1Trigger.L1TCalorimeter.caloParams_cfi
caloStage2Params = L1Trigger.L1TCalorimeter.caloParams_cfi.caloParams.clone()

# towers
caloStage2Params.towerLsbH        = cms.double(0.5)
caloStage2Params.towerLsbE        = cms.double(0.5)
caloStage2Params.towerLsbSum      = cms.double(0.5)
caloStage2Params.towerNBitsH      = cms.int32(8)
caloStage2Params.towerNBitsE      = cms.int32(8)
caloStage2Params.towerNBitsSum    = cms.int32(9)
caloStage2Params.towerNBitsRatio  = cms.int32(3)
caloStage2Params.towerEncoding    = cms.bool(True)

# regions
caloStage2Params.regionLsb        = cms.double(0.5)
caloStage2Params.regionPUSType    = cms.string("None")
caloStage2Params.regionPUSParams  = cms.vdouble()

# EG
caloStage2Params.egLsb                      = cms.double(0.5)
caloStage2Params.egSeedThreshold            = cms.double(2.)
caloStage2Params.egNeighbourThreshold       = cms.double(1.)
caloStage2Params.egHcalThreshold            = cms.double(0.)
caloStage2Params.egTrimmingLUTFile          = cms.FileInPath("L1Trigger/L1TCalorimeter/data/egTrimmingLUT_10_v16.01.19.txt")
caloStage2Params.egMaxHcalEt                = cms.double(0.)
caloStage2Params.egMaxPtHOverE          = cms.double(128.)
caloStage2Params.egMaxHOverELUTFile         = cms.FileInPath("L1Trigger/L1TCalorimeter/data/HoverEIdentification_0.995_v15.12.23.txt")
caloStage2Params.egCompressShapesLUTFile    = cms.FileInPath("L1Trigger/L1TCalorimeter/data/egCompressLUT_v4.txt")
caloStage2Params.egShapeIdType              = cms.string("compressed")
caloStage2Params.egShapeIdVersion           = cms.uint32(0)
caloStage2Params.egShapeIdLUTFile           = cms.FileInPath("L1Trigger/L1TCalorimeter/data/shapeIdentification_adapt0.99_compressedieta_compressedE_compressedshape_v15.12.08.txt")#Not used any more in the current emulator version, merged with calibration LUT
caloStage2Params.egPUSType                  = cms.string("None")
caloStage2Params.egIsolationType            = cms.string("compressed")
#caloStage2Params.egIsoLUTFile               = cms.FileInPath("L1Trigger/L1TCalorimeter/data/IsoIdentification_adapt_extrap_v16.07.29.txt")
caloStage2Params.egIsoLUTFile               = cms.FileInPath("L1Trigger/L1TCalorimeter/data/IsoIdentification_adapt_extrap_FW_v16.08.08.txt") # new SK Sep '16
caloStage2Params.egIsoAreaNrTowersEta       = cms.uint32(2)
caloStage2Params.egIsoAreaNrTowersPhi       = cms.uint32(4)
caloStage2Params.egIsoVetoNrTowersPhi       = cms.uint32(2)
#caloStage2Params.egIsoPUEstTowerGranularity = cms.uint32(1)
#caloStage2Params.egIsoMaxEtaAbsForTowerSum  = cms.uint32(4)
#caloStage2Params.egIsoMaxEtaAbsForIsoSum    = cms.uint32(27)
caloStage2Params.egPUSParams                = cms.vdouble(1,4,32) #Isolation window in firmware goes up to abs(ieta)=32 for now
caloStage2Params.egCalibrationType          = cms.string("compressed")
caloStage2Params.egCalibrationVersion       = cms.uint32(0)
#caloStage2Params.egCalibrationLUTFile       = cms.FileInPath("L1Trigger/L1TCalorimeter/data/corrections_Trimming10_compressedieta_compressedE_compressedshape_v16.03.14_shapeIdentification_adapt0.99_compressedieta_compressedE_compressedshape_v15.12.08.txt")
# 2018 EG calib
caloStage2Params.egCalibrationLUTFile       = cms.FileInPath("L1Trigger/L1TCalorimeter/data/corrections_Trimming10_compressedieta_compressedE_compressedshape_PANTELIS_v2_NEW_CALIBRATIONS_withShape_v17.04.04.txt")


# Tau
caloStage2Params.tauLsb                        = cms.double(0.5)
caloStage2Params.tauSeedThreshold              = cms.double(0.)
caloStage2Params.tauNeighbourThreshold         = cms.double(0.)
caloStage2Params.tauIsoAreaNrTowersEta         = cms.uint32(2)
caloStage2Params.tauIsoAreaNrTowersPhi         = cms.uint32(4)
caloStage2Params.tauIsoVetoNrTowersPhi         = cms.uint32(2)
caloStage2Params.tauPUSType                    = cms.string("None")
caloStage2Params.tauIsoLUTFile                 = cms.FileInPath("L1Trigger/L1TCalorimeter/data/Tau_Iso_LUT_Option_22_NewLayer1Calibration_SK1616_noCompressionBlock_FW_v6.2.0.txt")
caloStage2Params.tauIsoLUTFile2                = cms.FileInPath("L1Trigger/L1TCalorimeter/data/Tau_Iso_LUT_Option_22_NewLayer1Calibration_SK1616_noCompressionBlock_FW_v6.2.0.txt")
caloStage2Params.tauCalibrationLUTFile         = cms.FileInPath("L1Trigger/L1TCalorimeter/data/Tau_Calibration_LUT_NewLayer1Calibration_SK1616_FW_v11.0.0.txt")
caloStage2Params.tauCompressLUTFile            = cms.FileInPath("L1Trigger/L1TCalorimeter/data/tauCompressAllLUT_12bit_v3.txt")
caloStage2Params.tauPUSParams                  = cms.vdouble(1,4,32)

# jets
caloStage2Params.jetLsb                = cms.double(0.5)
caloStage2Params.jetSeedThreshold      = cms.double(4.0)
caloStage2Params.jetNeighbourThreshold = cms.double(0.)
caloStage2Params.jetPUSType            = cms.string("ChunkyDonut")

# Calibration options
# function6PtParams22EtaBins or None
#caloStage2Params.jetCalibrationType    = cms.string("None")
#caloStage2Params.jetCalibrationType = cms.string("function8PtParams22EtaBins")
caloStage2Params.jetCalibrationType = cms.string("LUT")

#Vector with 6 parameters for eta bin, from low eta to high
# 1,0,1,0,1,1 gives no correction
# must be in this form as may require > 255 arguments

# Or vector with 8 parameters, which caps correction value below given pT
# as 6 parameters, but last two are max correction and L1 pT below which cap is applied, respectively

jetCalibParamsVector = cms.vdouble()
jetCalibParamsVector.extend([
        1,0,1,0,1,1,1.36123039014,1024,
        1,0,1,0,1,1,1.37830172245,1024,
        1,0,1,0,1,1,1.37157036457,1024,
        1,0,1,0,1,1,1.42460009989,1024,
        10.1179757811,-697.422255848,55.9767511168,599.040770412,0.00930772659892,-21.9921521313,1.77585386314,24.1202894336,
        12.2578170485,-736.96846599,45.3225355911,848.976802835,0.00946235693865,-21.7970133915,2.04623980351,19.6049149791,
        14.0198255047,-769.175319944,38.687351315,1072.9785137,0.00951954709279,-21.6277409602,2.08021511285,22.265051562,
        14.119589176,-766.199501821,38.7767169666,1059.63374337,0.00952979125289,-21.6477483043,2.05901166216,23.8125466978,
        13.7594864391,-761.860391454,39.9060363401,1019.30588542,0.00952105483129,-21.6814176696,2.03808638982,22.2127275989,
        10.2635352836,-466.890522023,32.5408463829,2429.03382746,0.0111274121697,-22.0890253377,2.04880080215,22.5083699943,
        5.46086027683,-150.888778124,18.3292242153,16968.6469599,0.0147496053457,-22.4089831889,2.08107691501,22.4129703515,
        5.46086027683,-150.888778124,18.3292242153,16968.6469599,0.0147496053457,-22.4089831889,2.08107691501,22.4129703515,
        10.2635352836,-466.890522023,32.5408463829,2429.03382746,0.0111274121697,-22.0890253377,2.04880080215,22.5083699943,
        13.7594864391,-761.860391454,39.9060363401,1019.30588542,0.00952105483129,-21.6814176696,2.03808638982,22.2127275989,
        14.119589176,-766.199501821,38.7767169666,1059.63374337,0.00952979125289,-21.6477483043,2.05901166216,23.8125466978,
        14.0198255047,-769.175319944,38.687351315,1072.9785137,0.00951954709279,-21.6277409602,2.08021511285,22.265051562,
        12.2578170485,-736.96846599,45.3225355911,848.976802835,0.00946235693865,-21.7970133915,2.04623980351,19.6049149791,
        10.1179757811,-697.422255848,55.9767511168,599.040770412,0.00930772659892,-21.9921521313,1.77585386314,24.1202894336,
        1,0,1,0,1,1,1.42460009989,1024,
        1,0,1,0,1,1,1.37157036457,1024,
        1,0,1,0,1,1,1.37830172245,1024,
        1,0,1,0,1,1,1.36123039014,1024
])
caloStage2Params.jetCalibrationParams  = jetCalibParamsVector 

#caloStage2Params.jetCompressPtLUTFile     = cms.FileInPath("L1Trigger/L1TCalorimeter/data/lut_pt_compress.txt")
#caloStage2Params.jetCompressEtaLUTFile    = cms.FileInPath("L1Trigger/L1TCalorimeter/data/lut_30to40_hfHighPt_experiment2_changeLimits_eta.txt")
#caloStage2Params.jetCalibrationLUTFile    = cms.FileInPath("L1Trigger/L1TCalorimeter/data/lut_30to40_hfHighPt_experiment2_changeLimits_add_mult.txt")
caloStage2Params.jetCompressPtLUTFile     = cms.FileInPath("L1Trigger/L1TCalorimeter/data/lut_pt_compress_2017v1.txt")
caloStage2Params.jetCompressEtaLUTFile    = cms.FileInPath("L1Trigger/L1TCalorimeter/data/lut_eta_compress_2017v1.txt")
caloStage2Params.jetCalibrationLUTFile    = cms.FileInPath("L1Trigger/L1TCalorimeter/data/lut_calib_2018v1.txt")



# sums: 0=ET, 1=HT, 2=MET, 3=MHT
caloStage2Params.etSumLsb                = cms.double(0.5)
caloStage2Params.etSumEtaMin             = cms.vint32(1, 1, 1, 1, 1)
caloStage2Params.etSumEtaMax             = cms.vint32(28,  28,  28,  28, 27)
caloStage2Params.etSumEtThreshold        = cms.vdouble(0.,  30.,  0.,  30., 0.)

caloStage2Params.etSumXCalibrationLUTFile         = cms.FileInPath("L1Trigger/L1TCalorimeter/data/lut_etSumPUS_dummy.txt")
caloStage2Params.etSumYCalibrationLUTFile         = cms.FileInPath("L1Trigger/L1TCalorimeter/data/lut_etSumPUS_dummy.txt")
caloStage2Params.etSumEttCalibrationLUTFile       = cms.FileInPath("L1Trigger/L1TCalorimeter/data/lut_etSumPUS_dummy.txt")
caloStage2Params.etSumEcalSumCalibrationLUTFile   = cms.FileInPath("L1Trigger/L1TCalorimeter/data/lut_etSumPUS_dummy.txt")


# Layer 1 LUT specification
#
# Et-dependent scale factors
# ECal/HCal scale factors will be a 9*28 array:
#   28 eta scale factors (1-28)
#   in 9 ET bins (10, 15, 20, 25, 30, 35, 40, 45, Max)
#  So, index = etBin*28+ieta




# 2018 Layer 1 SF
caloStage2Params.layer1ECalScaleETBins = cms.vint32([6, 9, 12, 15, 20, 25, 30, 35, 40, 45, 55, 70, 256])
caloStage2Params.layer1ECalScaleFactors = cms.vdouble([
        1.128436, 1.102229, 1.128385, 1.127897, 1.142444, 1.115476, 1.104283, 1.124583, 1.115929, 1.115196, 1.130342, 1.127173, 1.130640, 1.125474, 1.126652, 1.143535, 1.148905, 1.309035, 1.156021, 1.292685, 1.314302, 1.327634, 1.341229, 1.364885, 1.411117, 1.432419, 1.288526, 1.082139, 
        1.078545, 1.072734, 1.075464, 1.081920, 1.078434, 1.072281, 1.079780, 1.082043, 1.094741, 1.074544, 1.082784, 1.084089, 1.086375, 1.099718, 1.092858, 1.092855, 1.105166, 1.256155, 1.126301, 1.215671, 1.226302, 1.268900, 1.281721, 1.310629, 1.356976, 1.386428, 1.220159, 1.066925, 
        1.052366, 1.053986, 1.055250, 1.051033, 1.055017, 1.062249, 1.059624, 1.065355, 1.062623, 1.054089, 1.060477, 1.074504, 1.075570, 1.078549, 1.071588, 1.080279, 1.078463, 1.211087, 1.103915, 1.186517, 1.194161, 1.234868, 1.250080, 1.274639, 1.327394, 1.362218, 1.161404, 1.062366, 
        1.044640, 1.043507, 1.046185, 1.042067, 1.042425, 1.044121, 1.050677, 1.051604, 1.046070, 1.040140, 1.052732, 1.055652, 1.057201, 1.062982, 1.059512, 1.054542, 1.063873, 1.189094, 1.091948, 1.165298, 1.177338, 1.213632, 1.223587, 1.259376, 1.312025, 1.330172, 1.160220, 1.059058, 
        1.032947, 1.033877, 1.036016, 1.036056, 1.037819, 1.036489, 1.040341, 1.035373, 1.042736, 1.030510, 1.039291, 1.043943, 1.051946, 1.049653, 1.045154, 1.048874, 1.043392, 1.146608, 1.083743, 1.161479, 1.164940, 1.197187, 1.229915, 1.238886, 1.289410, 1.344620, 1.078591, 1.051894, 
        1.025813, 1.028301, 1.026054, 1.032050, 1.029899, 1.032383, 1.033763, 1.034211, 1.033892, 1.023902, 1.034960, 1.039866, 1.039984, 1.042478, 1.041047, 1.044143, 1.038748, 1.146814, 1.069148, 1.134356, 1.147952, 1.175102, 1.202532, 1.234549, 1.285897, 1.280056, 1.055845, 1.050155, 
        1.025370, 1.024465, 1.023378, 1.024989, 1.026322, 1.025140, 1.026122, 1.028451, 1.029161, 1.020083, 1.031555, 1.032971, 1.036222, 1.042410, 1.038053, 1.036796, 1.037195, 1.123576, 1.071556, 1.129229, 1.129561, 1.170449, 1.190240, 1.218357, 1.270482, 1.302586, 1.047321, 1.049100, 
        1.018591, 1.019825, 1.020823, 1.019265, 1.021761, 1.021521, 1.024053, 1.024121, 1.024979, 1.015315, 1.026035, 1.028734, 1.030409, 1.031414, 1.030694, 1.033450, 1.035642, 1.103688, 1.066969, 1.117955, 1.135950, 1.163170, 1.180714, 1.228736, 1.254963, 1.307361, 1.047123, 1.047264, 
        1.017483, 1.016714, 1.018925, 1.017087, 1.020438, 1.018852, 1.020796, 1.022534, 1.023495, 1.013378, 1.024097, 1.026067, 1.029037, 1.030731, 1.028759, 1.032480, 1.034680, 1.101491, 1.069770, 1.110644, 1.129222, 1.147881, 1.176695, 1.219110, 1.253033, 1.308691, 1.040706, 1.046607, 
        1.015432, 1.014445, 1.016057, 1.014908, 1.019115, 1.016567, 1.020411, 1.019852, 1.020255, 1.010779, 1.023433, 1.023674, 1.027479, 1.027385, 1.027332, 1.027537, 1.029061, 1.091079, 1.063278, 1.108876, 1.122727, 1.171282, 1.172058, 1.211259, 1.245839, 1.303968, 1.033863, 1.047743, 
        1.014370, 1.013304, 1.013397, 1.014261, 1.013673, 1.013183, 1.018534, 1.016581, 1.017015, 1.008220, 1.019515, 1.021560, 1.024502, 1.025611, 1.025905, 1.025863, 1.027252, 1.085230, 1.063040, 1.112256, 1.116617, 1.140393, 1.159214, 1.191434, 1.240601, 1.268525, 1.033247, 1.042853, 
        1.010174, 1.009843, 1.011520, 1.011041, 1.012957, 1.009075, 1.013178, 1.013301, 1.015033, 1.005133, 1.017533, 1.018564, 1.020319, 1.022634, 1.022429, 1.022338, 1.025613, 1.077639, 1.057895, 1.107098, 1.111157, 1.136106, 1.161737, 1.179259, 1.232736, 1.290141, 1.018941, 1.014733, 
        1.000302, 1.007651, 1.000751, 1.007791, 1.008949, 1.005394, 1.009599, 1.010180, 1.010865, 1.001827, 1.012447, 1.015231, 1.019545, 1.020611, 1.022404, 1.019032, 1.023113, 1.065127, 1.054688, 1.102754, 1.106151, 1.125574, 1.134480, 1.180965, 1.231939, 1.277289, 1.018941, 1.014733
    ])

"""
caloStage2Params.layer1ECalScaleETBins = cms.vint32([10, 15, 20, 25, 30, 35, 40, 45, 256])
caloStage2Params.layer1ECalScaleFactors = cms.vdouble([
    1.1847, 1.16759, 1.17779, 1.19955, 1.21125, 1.214, 1.21503, 1.22515, 1.24151, 1.27836, 1.30292, 1.33526, 1.42338, 1.4931, 1.49597, 1.50405, 1.52785, 1.81552, 1.59856, 1.75692, 1.76496, 1.77562, 1.69527, 1.66827, 1.61861, 1.56645, 1.56645, 1.56645,
    1.1351, 1.12589, 1.12834, 1.13725, 1.14408, 1.1494, 1.14296, 1.14852, 1.1578, 1.17634, 1.18038, 1.19386, 1.23758, 1.27605, 1.27818, 1.28195, 1.34881, 1.71053, 1.37338, 1.52571, 1.54801, 1.53316, 1.4397, 1.40497, 1.37743, 1.33914, 1.33914, 1.33914,
    1.18043, 1.17823, 1.1751, 1.17608, 1.19152, 1.196, 1.20125, 1.2068, 1.22584, 1.22476, 1.22395, 1.22302, 1.25137, 1.28097, 1.29871, 1.2862, 1.33489, 1.60937, 1.28365, 1.41367, 1.42521, 1.42041, 1.36784, 1.34922, 1.32754, 1.29825, 1.29825, 1.29825,
    1.11664, 1.11852, 1.11861, 1.12367, 1.12405, 1.14814, 1.14304, 1.15337, 1.16607, 1.18698, 1.17048, 1.17463, 1.2185, 1.23842, 1.23214, 1.24744, 1.30047, 1.47152, 1.22868, 1.33121, 1.34841, 1.35178, 1.30048, 1.28537, 1.27012, 1.24159, 1.24159, 1.24159,
    1.08422, 1.08146, 1.08706, 1.08906, 1.08636, 1.10092, 1.10363, 1.11102, 1.1186, 1.13301, 1.12369, 1.14377, 1.16477, 1.17801, 1.18782, 1.17168, 1.24593, 1.36835, 1.20252, 1.28349, 1.29828, 1.30328, 1.26848, 1.25817, 1.2464, 1.22259, 1.22259, 1.22259,
    1.07444, 1.06774, 1.06883, 1.0707, 1.07881, 1.08859, 1.08285, 1.08747, 1.09736, 1.10678, 1.10008, 1.10717, 1.12858, 1.15383, 1.15826, 1.14855, 1.19911, 1.32567, 1.17553, 1.25976, 1.27926, 1.28459, 1.24524, 1.23706, 1.22597, 1.20006, 1.20006, 1.20006,
    1.06224, 1.05968, 1.05767, 1.06254, 1.06729, 1.0691, 1.07125, 1.07312, 1.08124, 1.08966, 1.08695, 1.08826, 1.10611, 1.13115, 1.12641, 1.13093, 1.17074, 1.28958, 1.16217, 1.22844, 1.24812, 1.25352, 1.22065, 1.21287, 1.20544, 1.18344, 1.18344, 1.18344,
    1.03589, 1.03224, 1.03229, 1.03623, 1.03979, 1.04403, 1.04574, 1.049, 1.04821, 1.06183, 1.0588, 1.06655, 1.08582, 1.10289, 1.10052, 1.10506, 1.143, 1.27373, 1.1459, 1.2156, 1.23455, 1.23968, 1.20753, 1.20127, 1.19629, 1.16809, 1.16809, 1.16809,
    1.03456, 1.02955, 1.03079, 1.03509, 1.03949, 1.0437, 1.04236, 1.04486, 1.0517, 1.05864, 1.05516, 1.06167, 1.07738, 1.0985, 1.09317, 1.09559, 1.13557, 1.26076, 1.14118, 1.20545, 1.22137, 1.22802, 1.19936, 1.19676, 1.19088, 1.16709, 1.16709, 1.16709,
    ])
caloStage2Params.layer1HCalScaleETBins = cms.vint32([10, 15, 20, 25, 30, 35, 40, 45, 256])
caloStage2Params.layer1HCalScaleFactors = cms.vdouble([
    1.511112, 1.519900, 1.499483, 1.488560, 1.528111, 1.475114, 1.476616, 1.514163, 1.515306, 1.542464, 1.511663, 1.593745, 1.493667, 1.485315, 1.419925, 1.349169, 1.312518, 1.423302, 1.478461, 1.525868, 1.525868, 1.525868, 1.525868, 1.525868, 1.525868, 1.525868, 1.525868, 1.525868,
    1.383350, 1.365700, 1.368470, 1.354610, 1.348480, 1.329720, 1.272250, 1.301710, 1.322210, 1.360860, 1.333850, 1.392200, 1.403060, 1.394870, 1.322050, 1.244570, 1.206910, 1.321870, 1.344160, 1.403270, 1.403270, 1.403270, 1.403270, 1.403270, 1.403270, 1.403270, 1.403270, 1.403270,
    1.245690, 1.238320, 1.245420, 1.234830, 1.243730, 1.249790, 1.179450, 1.213620, 1.219030, 1.252130, 1.209560, 1.250710, 1.280490, 1.262800, 1.254060, 1.186810, 1.127830, 1.260000, 1.275140, 1.305850, 1.305850, 1.305850, 1.305850, 1.305850, 1.305850, 1.305850, 1.305850, 1.305850,
    1.189940, 1.189120, 1.177120, 1.179690, 1.185510, 1.150590, 1.151830, 1.167860, 1.154310, 1.163190, 1.161700, 1.136100, 1.161870, 1.195050, 1.153910, 1.117900, 1.106750, 1.208120, 1.160020, 1.232800, 1.232800, 1.232800, 1.232800, 1.232800, 1.232800, 1.232800, 1.232800, 1.232800,
    1.122540, 1.129520, 1.125080, 1.115150, 1.118250, 1.096190, 1.108170, 1.087490, 1.109750, 1.099780, 1.081000, 1.050610, 1.078270, 1.079460, 1.047740, 1.041400, 1.041750, 1.116880, 1.097730, 1.125780, 1.125780, 1.125780, 1.125780, 1.125780, 1.125780, 1.125780, 1.125780, 1.125780,
    1.110470, 1.117340, 1.115980, 1.088490, 1.088260, 1.078230, 1.062720, 1.054690, 1.053270, 1.086640, 1.050620, 1.038470, 1.046440, 1.059130, 1.012240, 1.039030, 1.036040, 1.088460, 1.078880, 1.090600, 1.090600, 1.090600, 1.090600, 1.090600, 1.090600, 1.090600, 1.090600, 1.090600,
    1.115970, 1.111010, 1.113170, 1.079390, 1.076850, 1.063730, 1.039300, 1.049910, 1.040100, 1.025820, 1.015830, 1.015850, 1.010810, 1.014210, 0.980321, 1.023580, 1.045990, 1.073220, 1.057750, 1.059850, 1.059850, 1.059850, 1.059850, 1.059850, 1.059850, 1.059850, 1.059850, 1.059850,
    1.061180, 1.059770, 1.071210, 1.064420, 1.065340, 1.043070, 1.041400, 1.022680, 1.017410, 1.017690, 1.005610, 1.006360, 0.999420, 0.990866, 0.986723, 0.989036, 0.995116, 1.045620, 1.024330, 1.040660, 1.040660, 1.040660, 1.040660, 1.040660, 1.040660, 1.040660, 1.040660, 1.040660,
    1.083150, 1.067090, 1.083180, 1.061010, 1.075640, 1.051640, 1.038760, 1.042670, 1.010910, 1.011580, 1.006560, 0.984468, 0.986642, 0.985799, 0.968133, 1.000290, 1.011210, 1.046690, 1.016670, 1.020470, 1.020470, 1.020470, 1.020470, 1.020470, 1.020470, 1.020470, 1.020470, 1.020470,
    ])"""

caloStage2Params.layer1HCalScaleETBins = cms.vint32([6, 9, 12, 15, 20, 25, 30, 35, 40, 45, 55, 70, 256])
caloStage2Params.layer1HCalScaleFactors = cms.vdouble([
        1.691347, 1.704095, 1.729441, 1.735242, 1.726367, 1.780424, 1.794996, 1.815904, 1.817388, 1.894632, 1.932656, 1.957527, 1.970890, 2.005818, 2.041546, 2.042775, 1.989288, 1.594904, 
        1.659821, 1.676038, 1.495936, 1.505035, 1.512590, 1.511470, 1.494893, 1.378435, 1.430994, 1.500227, 1.531796, 1.547539, 1.559295, 1.561478, 1.568922, 1.601485, 1.616591, 1.620739, 
        1.642884, 1.678420, 1.692987, 1.728681, 1.728957, 1.766650, 1.782739, 1.782875, 1.751371, 1.431918, 1.487225, 1.483881, 1.336485, 1.349895, 1.363924, 1.375728, 1.377818, 1.310078, 
        1.334588, 1.399686, 1.465418, 1.462800, 1.475840, 1.474735, 1.474407, 1.506928, 1.526279, 1.524000, 1.532718, 1.583398, 1.608380, 1.623528, 1.619634, 1.646501, 1.667856, 1.674628, 
        1.635381, 1.350235, 1.394938, 1.383940, 1.244552, 1.256971, 1.261180, 1.282746, 1.279512, 1.221092, 1.241831, 1.351526, 1.390201, 1.404198, 1.416259, 1.404045, 1.418265, 1.437914, 
        1.450857, 1.463511, 1.462653, 1.501891, 1.518896, 1.548252, 1.545831, 1.565901, 1.574314, 1.575115, 1.557629, 1.301893, 1.326949, 1.312526, 1.197573, 1.210304, 1.222283, 1.239081, 
        1.240673, 1.185591, 1.207651, 1.275166, 1.314260, 1.335228, 1.340603, 1.323027, 1.324793, 1.347954, 1.349916, 1.363145, 1.359628, 1.402624, 1.416518, 1.457202, 1.461053, 1.484090, 
        1.500787, 1.498450, 1.471731, 1.215732, 1.253565, 1.243598, 1.157168, 1.164428, 1.175435, 1.189310, 1.192682, 1.142038, 1.162810, 1.230426, 1.262901, 1.265380, 1.274364, 1.276111, 
        1.282349, 1.291748, 1.305521, 1.301818, 1.305124, 1.336506, 1.345742, 1.357458, 1.370139, 1.381995, 1.394554, 1.388952, 1.363805, 1.166810, 1.204780, 1.193913, 1.118331, 1.124657, 
        1.136138, 1.148564, 1.147392, 1.085564, 1.109949, 1.184837, 1.221607, 1.219692, 1.235950, 1.230444, 1.234908, 1.245100, 1.256813, 1.252608, 1.263569, 1.284188, 1.300083, 1.309901, 
        1.312849, 1.335500, 1.339967, 1.328269, 1.309282, 1.128239, 1.173002, 1.163030, 1.077388, 1.087037, 1.085620, 1.099773, 1.097418, 1.047416, 1.080447, 1.135984, 1.186335, 1.189457, 
        1.186903, 1.191054, 1.192951, 1.218812, 1.222226, 1.220196, 1.221331, 1.264243, 1.284869, 1.277098, 1.263366, 1.276293, 1.291829, 1.275918, 1.248086, 1.095700, 1.143874, 1.132783, 
        1.054939, 1.055922, 1.055405, 1.058330, 1.062463, 1.012972, 1.028538, 1.089975, 1.155949, 1.153120, 1.157186, 1.163320, 1.157607, 1.174722, 1.181157, 1.179473, 1.186948, 1.192614, 
        1.207973, 1.215075, 1.252322, 1.231549, 1.241483, 1.224214, 1.207592, 1.069829, 1.112551, 1.107158, 1.025349, 1.026181, 1.028466, 1.035129, 1.030918, 0.977843, 1.004295, 1.075236, 
        1.122942, 1.124839, 1.130900, 1.139241, 1.134602, 1.141732, 1.154381, 1.154366, 1.162207, 1.167863, 1.182334, 1.189497, 1.179567, 1.185553, 1.205978, 1.188532, 1.154839, 1.058371, 
        1.096597, 1.086545, 0.997724, 1.000690, 1.005683, 1.009107, 1.006028, 0.962736, 0.974019, 1.035748, 1.094997, 1.098600, 1.101567, 1.102895, 1.106445, 1.113255, 1.114956, 1.118930, 
        1.128154, 1.135288, 1.145308, 1.151612, 1.142554, 1.153640, 1.154025, 1.138100, 1.127446, 1.034945, 1.069153, 1.062188, 0.977909, 0.972598, 0.972539, 0.978454, 0.975065, 0.941113, 
        0.948722, 1.004971, 1.055020, 1.054883, 1.059317, 1.061911, 1.062005, 1.066707, 1.074156, 1.064278, 1.072810, 1.076579, 1.084072, 1.091055, 1.090640, 1.086634, 1.095179, 1.075771, 
        1.051884, 1.005930, 1.033331, 1.024734, 0.943637, 0.941986, 0.937779, 0.943865, 0.928477, 0.902234, 0.908232, 0.960607, 1.005841, 1.011405, 1.012527, 1.015557, 1.014508, 1.020877, 
        1.019076, 1.015173, 1.015651, 1.019594, 1.026845, 1.024959, 1.025915, 1.029455, 1.017985, 1.016933, 0.989723, 0.977768, 0.993744, 0.985200, 0.907247, 0.903328, 0.912164, 0.898908, 
        0.886431, 0.851162, 0.863541, 0.890523
    ])


# HF 1x1 scale factors will be a 5*12 array:
#  12 eta scale factors (30-41)
#  in 5 REAL ET bins (5, 20, 30, 50, Max)
#  So, index = etBin*12+ietaHF
caloStage2Params.layer1HFScaleETBins = cms.vint32([5, 20, 30, 50, 256])
caloStage2Params.layer1HFScaleFactors = cms.vdouble([
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 
    1.767080, 1.767080, 1.755186, 1.769951, 1.763527, 1.791043, 1.898787, 1.982235, 2.071074, 2.193011, 2.356886, 2.403384, 
    2.170477, 2.170477, 2.123540, 2.019866, 1.907698, 1.963179, 1.989122, 2.035251, 2.184642, 2.436399, 2.810884, 2.923750, 
    1.943941, 1.943941, 1.899826, 1.813950, 1.714978, 1.736184, 1.785928, 1.834211, 1.944230, 2.153565, 2.720887, 2.749795, 
    1.679984, 1.679984, 1.669753, 1.601871, 1.547276, 1.577805, 1.611497, 1.670184, 1.775022, 1.937061, 2.488311, 2.618629, 
    ])
