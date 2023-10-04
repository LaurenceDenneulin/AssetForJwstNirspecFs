using PyPlot

using AssetForJwstNirspecFs

hpz = 1e25 #regularization hyperparameter

nexp = 3 #number of expozures

λ,z = AssetForJwstNirspecFs.extract_spectrum("pathtodata/beginingofthefilename", 
                     "endofthefilename.fits",
                     nexp,
                     hpz)

plot(λ,z)

