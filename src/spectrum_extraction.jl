

"""
    z=extract_spectrum()
    
    performs the spectral extraction from inegrations files. 
    
    -'filename_beg' is the first part of the String before the integration number;
    -'filename_end' is the last part of the String after the integration number;
    -'nfiles' is the amount of integrations,
    -'hpz' is the hyperparameter for the spectrum regularization
"""
function extract_spectrum(filename_beg::AbstractString,
                             filename_end::AbstractString, 
                             nfiles::Integer, 
                             hpz::AbstractFloat;
                             ker = CatmullRomSpline(Float64, Flat), 
                             hpb::AbstractFloat=0.,
                             sky_sub=true, 
                             order::Integer=3,
                             max_iter = 10,
                             psf_params_bnds=[(0.0, 0.1);(0., 1.)],
                             center_bnds=1.)

    data, wgt, lambda, ρ, cent = load_data(filename_beg,filename_end, nfiles; order=order, sky_sub=sky_sub)
    
    Rz=hpz*tikhonov()
    Rb=hpb*tikhonov()

    m,n,=size(data)
    αinit = 0.1*maximum(lambda)^(-2)
    βinit = 0.25
    ρinit =cent[1]
    Nz= n

    λmin=minimum(lambda[lambda.!=0])
    λmax=maximum(lambda[lambda.!=0])
    λz=range(λmin, length=Nz, stop=λmax);
    F = SparseInterpolator(ker, lambda, λz)
    
    z=zeros(Nz); 
    parinit = [ αinit, βinit]
    centinit = copy(ρinit)
    h=ASSET.chromwmwGaussianPSF(parinit)
    if sky_sub
        b=undef
    else
        b=ASSET.BkgMdl(zeros(m,n), Rb)
    end
    D=ASSET.CalibratedData(data, wgt, ρ, lambda)

    psf_center_bnds=[(ρinit[k] - center_bnds, ρinit[k] + center_bnds) for k=1:length(ρinit)]

       z,hf,cf =ASSET.extract_spectrum!(z, F, centinit, h, D ,Rz,b; extract_kwds=(psf_params_bnds=psf_params_bnds, psf_center_bnds=psf_center_bnds),max_iter=max_iter)

    return λz,z
end 
