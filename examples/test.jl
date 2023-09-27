using PyPlot

using AssetForJwstNirspecFs

hpz = 1e25 

λz2,z2 = AssetForJwstNirspecFs.extract_spectrum("AZ84/CAL/jw01231002001_04102_0000", 
                     "_nrs1_cal.fits",
                     3,
                     hpz)

λz4,z4 = AssetForJwstNirspecFs.extract_spectrum("AZ84/CAL/jw01231002001_04104_0000", 
                     "_nrs1_cal.fits",
                     3,
                     hpz)

λz6,z6 = AssetForJwstNirspecFs.extract_spectrum("AZ84/CAL/jw01231002001_04106_0000", 
                     "_nrs1_cal.fits",
                     3,
                     hpz)

λz8,z8 = AssetForJwstNirspecFs.extract_spectrum("AZ84/CAL/jw01231002001_04108_0000", 
                     "_nrs1_cal.fits",
                     3,
                     hpz)


plot(λz2,z2)
plot(λz4,z4)
plot(λz6,z6)
plot(λz8,z8)
