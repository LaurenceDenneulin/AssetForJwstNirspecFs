"""
    z=extract_spectrum()
    
    performs the spectral extraction from inegrations files. 
    
    -'dir' is the path to the directory where the '_cal' data are stored (one directory per filter).
  
    -'hpz' is the hyperparameter for the spectrum regularization
"""
function extract_spectrum(dir::AbstractString,
                             hpz::AbstractFloat;
                             ker = CatmullRomSpline(Float64, Flat), 
                             hpb::AbstractFloat=0.,
                             sky_sub=true, 
                             order::Integer=3,
                             max_iter = 2,
                             psf_params_bnds=[(0.0, 0.1);(0., 1.)],
                             center_bnds=1.,
                             save=false)

    if save
        mkpath("save")
    end
    
    data, wgt, lambda, ρ, cent = load_data(dir; order=order, sky_sub=sky_sub, save=save)
    
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

    if save
        filenames=readdir(dir)
        numpos=findfirst("_nrs",filenames[1])[1]-1
        extpos=findfirst(".",filenames[1])[1]
        filename_beg=filenames[1][1:numpos-1]
        filename_end=filenames[1][numpos+1:extpos-1]   
        writedlm("save/"*filename_beg*"X"*filename_end[1:end-5]*"spectrum_hp=$hpz"*".txt", [λz z])
    end
    return λz,z
end 
