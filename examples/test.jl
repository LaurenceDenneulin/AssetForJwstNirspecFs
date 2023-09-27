using EasyFITS
using PyPlot
using DelimitedFiles
using InterpolationKernels
using LinearInterpolators
using Statistics
using LinearAlgebra
using LazyAlgebra
using PointSpreadFunctions
using OptimPackNextGen
using InverseProblem
import OptimPackNextGen.Powell.Bobyqa
import OptimPackNextGen.Brent

using ASSET
using AssetForJwstNirspecFs

for hpz in [1e22 5e22 1e23 5e23 1e24 5e24 1e25 5e25 1e26 5e26]
hpb=1e27
Rz=hpz*tikhonov()
Rb=hpb*tikhonov()

    data1, wgt1, bpm_map1, lambda1, pos1 = load_data("AZ84/CAL/jw01231002001_04108_00001_nrs1_cal.fits");
    data2, wgt2, bpm_map2, lambda2, pos2 = load_data("AZ84/CAL/jw01231002001_04108_00002_nrs1_cal.fits");
    data3, wgt3, bpm_map3, lambda3, pos3 = load_data("AZ84/CAL/jw01231002001_04108_00003_nrs1_cal.fits");


data=cat(data1,data2, data3, dims=3);
wgt=cat(wgt1 , wgt2 , wgt3 , dims=3);
wgt[isnan.(wgt)] .= 0.;
wgt[wgt .== Inf] .= 0.;
    
bpm_map=cat(bpm_map1, bpm_map2, bpm_map3, dims=3);
lambda=cat(lambda1, lambda2, lambda3, dims=3);
pos= cat(pos1,pos2,pos3, dims=1);

ρ = cat(geometric_calibration(data[:,:,1], bpm_map[:,:,1], pos1, order=3)[1],
        geometric_calibration(data[:,:,2], bpm_map[:,:,2], pos2, order=3)[1],
        geometric_calibration(data[:,:,3], bpm_map[:,:,3], pos3, order=3)[1], dims=3);
m,n =size(data)

αinit = 0.1*maximum(lambda)^(-2)
βinit = 0.25
ρinit, Nz = slit2cam(data, pos)
Nz = 400#le-lb
z=zeros(Nz);
b=zeros(m,n)


λmin=minimum(lambda[lambda.!=0])
λmax=maximum(lambda[lambda.!=0])
λz=range(λmin, length=Nz, stop=λmax);
#λz=range(0.54, step=1e-2, stop=5.39)
#Nz = λz.len
z=zeros(Nz);



F = SparseInterpolator(ker, lambda, λz)
parinit = [ αinit, βinit]
centinit = copy(ρinit)
h=ASSET.chromwmwGaussianPSF(parinit)
b=ASSET.BkgMdl(zeros(m,n), Rb)
display(h)
D=ASSET.CalibratedData(data, wgt, ρ, lambda)

       z,hf,cf =ASSET.extract_spectrum!(z, F, centinit, h, D ,Rz,b; extract_kwds=(psf_params_bnds=[(0.0, 0.1);(0., 1.)], psf_center_bnds=[(ρinit[1] -3, ρinit[1] +3),(ρinit[2] -3, ρinit[2] +3),(ρinit[3] -1, ρinit[3] +1)] ),max_iter=4)
         
     
display(hf)    
writedlm("jw01231002001_04108_global"*"_$hpz"*".txt", cat(λz,z, dims=2))

for k=1:3
    cent=[centinit[k]]
    Fk = SparseInterpolator(ker, lambda[:,:,k][:,:,:], λz)
    Dk = ASSET.CalibratedData(data[:,:,k][:,:,:] .- b.b, wgt[:,:,k][:,:,:], ρ[:,:,k][:,:,:], lambda[:,:,k][:,:,:])
    z,hf,cf =extract_spectrum!(z, Fk, cent, hf, Dk ,Rz; auto_calib=Val(false));

    writedlm("jw01231002001_04108_0000$k"*"_$hpz"*".txt", cat(λz,z, dims=2))
end



    data1, wgt1, bpm_map1, lambda1, pos1 = load_data("AZ84/CAL/jw01231002001_04102_00001_nrs1_cal.fits");
    data2, wgt2, bpm_map2, lambda2, pos2 = load_data("AZ84/CAL/jw01231002001_04102_00002_nrs1_cal.fits");
    data3, wgt3, bpm_map3, lambda3, pos3 = load_data("AZ84/CAL/jw01231002001_04102_00003_nrs1_cal.fits");


data=cat(data1,data2, data3, dims=3);
wgt=cat(wgt1 , wgt2 , wgt3 , dims=3);
wgt[isnan.(wgt)] .= 0.;
wgt[wgt .== Inf] .= 0.;
    
bpm_map=cat(bpm_map1, bpm_map2, bpm_map3, dims=3);
lambda=cat(lambda1, lambda2, lambda3, dims=3);
pos= cat(pos1,pos2,pos3, dims=1);

ρ = cat(geometric_calibration(data[:,:,1], bpm_map[:,:,1], pos1; order=2)[1],
        geometric_calibration(data[:,:,2], bpm_map[:,:,2], pos2; order=2)[1],
        geometric_calibration(data[:,:,3], bpm_map[:,:,3], pos3; order=2)[1], dims=3);
m,n =size(data)

αinit = 0.1*maximum(lambda)^(-2)
βinit = 0.25
ρinit, Nz= slit2cam(data, pos)
Nz = 1000#le-lb
z=zeros(Nz);
b=zeros(m,n)

λmin=minimum(lambda[lambda.!=0])
λmax=maximum(lambda[lambda.!=0])
λz=range(λmin, length=Nz, stop=λmax);
#λz=range(0.96, step=0.625e-3, stop=1.84)
#Nz = λz.len
#z=zeros(Nz);

F = SparseInterpolator(ker, lambda, λz)
parinit = [ αinit, βinit]
centinit = copy(ρinit)
#h=hf
h=ASSET.chromwmwGaussianPSF(parinit)
b=ASSET.BkgMdl(zeros(m,n), Rb)
D=ASSET.CalibratedData(data, wgt, ρ, lambda)

         z,hf,cf =ASSET.extract_spectrum!(z, F, centinit, h, D ,Rz,b; extract_kwds=(psf_params_bnds=[(0.0, 0.1);(0.1, 0.5)], psf_center_bnds=[(ρinit[1] -1, ρinit[1] +1),(ρinit[2] -1, ρinit[2] +1),(ρinit[3] -1, ρinit[3] +1)] ))
    
    
         #z,hf,cf =ASSET.extract_spectrum!(z, F, centinit, hf, D ,Rz,b;auto_calib=Val(false))
    
    
    writedlm("jw01231002001_04102_global"*"_$hpz"*".txt", cat(λz,z, dims=2))

display(hf)
for k=1:3
    cent=[centinit[k]]
    Fk = SparseInterpolator(ker, lambda[:,:,k][:,:,:], λz)
    Dk = ASSET.CalibratedData(data[:,:,k][:,:,:] .- b.b, wgt[:,:,k][:,:,:], ρ[:,:,k][:,:,:], lambda[:,:,k][:,:,:])
    z,hf,cf =extract_spectrum!(z, Fk, cent, hf, Dk ,Rz; auto_calib=Val(false));

    writedlm("jw01231002001_04102_0000$k"*"_$hpz"*".txt", cat(λz,z, dims=2))
end




    data1, wgt1, bpm_map1, lambda1, pos1 = load_data("AZ84/CAL/jw01231002001_04104_00001_nrs1_cal.fits");
    data2, wgt2, bpm_map2, lambda2, pos2 = load_data("AZ84/CAL/jw01231002001_04104_00002_nrs1_cal.fits");
    data3, wgt3, bpm_map3, lambda3, pos3 = load_data("AZ84/CAL/jw01231002001_04104_00003_nrs1_cal.fits");


data=cat(data1,data2, data3, dims=3);
wgt=cat(wgt1 , wgt2 , wgt3 , dims=3);
wgt[isnan.(wgt)] .= 0.;
wgt[wgt .== Inf] .= 0.;
    
bpm_map=cat(bpm_map1, bpm_map2, bpm_map3, dims=3);
lambda=cat(lambda1, lambda2, lambda3, dims=3);
pos= cat(pos1,pos2,pos3, dims=1);

ρ = cat(geometric_calibration(data[:,:,1], bpm_map[:,:,1], pos1; order=2)[1],
        geometric_calibration(data[:,:,2], bpm_map[:,:,2], pos2; order=2)[1],
        geometric_calibration(data[:,:,3], bpm_map[:,:,3], pos3; order=2)[1], dims=3);
m,n =size(data)

αinit = 0.1*maximum(lambda)^(-2)
βinit = 0.25
ρinit,Nz = slit2cam(data, pos)
Nz = 1000#le-lb
z=zeros(Nz);
b=zeros(m,n)

λmin=minimum(lambda[lambda.!=0])
λmax=maximum(lambda[lambda.!=0])
λz=range(λmin, length=Nz, stop=λmax);
#λz= range(1.64, step=0.125e-2, stop=3.07)
#Nz = λz.len
z=zeros(Nz);


F = SparseInterpolator(ker, lambda, λz)
parinit = [ αinit, βinit]
centinit = copy(ρinit)
#h=hf
h=ASSET.chromwmwGaussianPSF(parinit)
b=ASSET.BkgMdl(zeros(m,n), Rb)
D=ASSET.CalibratedData(data, wgt, ρ, lambda)

         z,hf,cf =ASSET.extract_spectrum!(z, F, centinit, h, D ,Rz,b; extract_kwds=(psf_params_bnds=[(0.0, 0.1);(0.1, 1.)], psf_center_bnds=[(ρinit[1] -3, ρinit[1] +3),(ρinit[2] -3, ρinit[2] +3),(ρinit[3] -1, ρinit[3] +1)] ))
    
    
         #  z,hf,cf =ASSET.extract_spectrum!(z, F, centinit, hf, D ,Rz,b;auto_calib=Val(false))
    
    
    writedlm("jw01231002001_04104_global"*"_$hpz"*".txt", cat(λz,z, dims=2))
display(hf)
for k=1:3
    cent=[centinit[k]]
    Fk = SparseInterpolator(ker, lambda[:,:,k][:,:,:], λz)
    Dk = ASSET.CalibratedData(data[:,:,k][:,:,:] .- b.b, wgt[:,:,k][:,:,:], ρ[:,:,k][:,:,:], lambda[:,:,k][:,:,:])
    z,hf,cf =extract_spectrum!(z, Fk, cent, hf, Dk ,Rz; auto_calib=Val(false));

    writedlm("jw01231002001_04104_0000$k"*"_$hpz"*".txt", cat(λz,z, dims=2))
end




    data1, wgt1, bpm_map1, lambda1, pos1 = load_data("AZ84/CAL/jw01231002001_04106_00001_nrs1_cal.fits");
    data2, wgt2, bpm_map2, lambda2, pos2 = load_data("AZ84/CAL/jw01231002001_04106_00002_nrs1_cal.fits");
    data3, wgt3, bpm_map3, lambda3, pos3 = load_data("AZ84/CAL/jw01231002001_04106_00003_nrs1_cal.fits");


data=cat(data1,data2, data3, dims=3);
wgt=cat(wgt1 , wgt2 , wgt3 , dims=3);
wgt[isnan.(wgt)] .= 0.;
wgt[wgt .== Inf] .= 0.;
    
bpm_map=cat(bpm_map1, bpm_map2, bpm_map3, dims=3);
lambda=cat(lambda1, lambda2, lambda3, dims=3);
pos= cat(pos1,pos2,pos3, dims=1);

ρ = cat(geometric_calibration(data[:,:,1], bpm_map[:,:,1], pos1; order=2)[1],
        geometric_calibration(data[:,:,2], bpm_map[:,:,2], pos2; order=2)[1],
        geometric_calibration(data[:,:,3], bpm_map[:,:,3], pos3; order=2)[1], dims=3);
m,n =size(data)

αinit = 0.1*maximum(lambda)^(-2)
βinit = 0.25
ρinit,Nz = slit2cam(data, pos)
Nz = 1000#Nz
z=zeros(Nz);
b=zeros(m,n)

λmin=minimum(lambda[lambda.!=0])
λmax=maximum(lambda[lambda.!=0])
λz=range(λmin, length=Nz, stop=λmax);
#λz=range(2.85, step=0.125e-2, stop=5.1)
#Nz = λz.len
#z=zeros(Nz);

F = SparseInterpolator(ker, lambda, λz)
parinit = [ αinit, βinit]
centinit = copy(ρinit)
#h=hf
h=ASSET.chromwmwGaussianPSF(parinit)
b=ASSET.BkgMdl(zeros(m,n), Rb)
D=ASSET.CalibratedData(data, wgt, ρ, lambda)


         z,hf,cf =ASSET.extract_spectrum!(z, F, centinit, h, D ,Rz,b; extract_kwds=(psf_params_bnds=[(0.0, 0.1);(0.1, 0.5)], psf_center_bnds=[(ρinit[1] -1, ρinit[1] +1),(ρinit[2] -1, ρinit[2] +1),(ρinit[3] -1, ρinit[3] +1)] ))

       #z,hf,cf =ASSET.extract_spectrum!(z, F, centinit, hf, D ,Rz,b;auto_calib=Val(false))
    

    writedlm("jw01231002001_04106_global"*"_$hpz"*".txt", cat(λz,z, dims=2))

display(hf)
for k=1:3
    cent=[centinit[k]]
    Fk = SparseInterpolator(ker, lambda[:,:,k][:,:,:], λz)
    Dk = ASSET.CalibratedData(data[:,:,k][:,:,:] .- b.b, wgt[:,:,k][:,:,:], ρ[:,:,k][:,:,:], lambda[:,:,k][:,:,:])
    z,hf,cf =extract_spectrum!(z, Fk, cent, hf, Dk ,Rz; auto_calib=Val(false));

    writedlm("jw01231002001_04106_0000$k"*"_$hpz"*".txt", cat(λz,z, dims=2))
end

end




