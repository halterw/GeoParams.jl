function Txx_rot,Tyy_rot,Txy_rot, η_ve, anifacve = LocalViscoElasticTrialStress_Newton(Exx_rot, Eyy_rot, Exy_rot, Txx_o_rot, Tyy_o_rot, Txy_o_rot, ηv, ηe, anifacv, aniface, tol, nitmax, noisy)

    Exx_ve    = Exx_rot .+ Txx_o_rot./(2*ηe)
    Eyy_ve    = Eyy_rot .+ Tyy_o_rot./(2*ηe)
    Exy_ve    = Exy_rot .+ Txy_o_rot./(2*ηe).*aniface 
    Ezz_ve    = -Exx_ve .-Eyy_ve
    # Construct initial guess for viscosity
    Eii         = 0.5*(Exx_rot.^2 .+ Eyy_rot.^2) .+ Exy_rot.^2 ./anifacv 
    ηve_n      = 1/(1/ηv .+ 1/ηe)                 # VE principal viscosity
    ηve_s      = 1/(anifacv./ηv .+ aniface./ηe)   # VE shear viscosity
    fxx0       = 0.
    fxy0       = 0.
    for iter=1:nitmax
        
        Txx       = 2*ηve_n .* Exx_ve
        Tyy       = 2*ηve_n .* Eyy_ve
        Txy       = 2*ηve_s .* Exy_ve
        Tzz       = -Txx .- Tyy 

        Y2_v      = 0.5*(Txx.^2 .+ Tyy.^2 .+ Tzz.^2) .+ Txy.^2 .*anifacv 

        # Residuals
        fxx       = 1 - ηve_n./ηe         .- ηve_n./ηv               
        fxy       = 1 - ηve_s./(ηe./aniface) .- ηve_s./(ηv./anifacv) 
        Wxx       = 2*(Exx_ve.*Txx .+ Eyy_ve.*Tyy .+ Ezz_ve.*Tzz)
        Wxy       = 2*Exy_ve.*Txy
        ieta_v    = 1 ./ ηv

        if iter==1 fxx0 = fxx  end
        if iter==1 fxy0 = fxy  end
        @printf("Newton ani VE It. %02d: fxx = %2.2e --- fxy = %2.2e\n", iter,  fxx, fxy)
        if (abs(fxx/fxx0)<tol || abs(fxx)<tol) && (abs(fxy/fxy0)<tol || abs(fxy)<tol)
            break
        end

        # Jacobian, check: AnisotropicVE_William_v2.ipynb
        dfxxde_n = Wxx .* ηve_n        .* ieta_v .* (1 - npwl) ./ (2 * Y2_v) - ieta_v - 1 ./ ηe;
        dfxxde_s = Wxy .* anifacv      .* ηve_n  .* ieta_v .* (1 - npwl) ./ Y2_v;
        dfxyde_n = Wxx .* anifacv      .* ηve_s  .* ieta_v .* (1 - npwl) ./ (2 * Y2_v);
        dfxyde_s = Wxy .* anifacv .^ 2 .* ηve_s  .* ieta_v .* (1 - npwl) ./ Y2_v - aniface ./ ηe - anifacv .* ieta_v;

        f      = [fxx; fxy]
        J      = [dfxxde_n dfxxde_s; dfxyde_n dfxyde_s]
        dx     = -J\f
        ηve_n += dx[1]
        ηve_s += dx[2]
    end

    ηve  = ηve_n
    anifacve = ηve_n/ηve_s
    
    # A posteriori stress evaluation
    Txx_rot =  2.0 * Exx_ve .* ηve
    Tyy_rot =  2.0 * Eyy_ve .* ηve
    Txy_rot =  2.0 * Exy_ve .* ηve ./ anifacve
    return Txx_rot,Tyy_rot,Txy_rot, η_ve, anifacve

end