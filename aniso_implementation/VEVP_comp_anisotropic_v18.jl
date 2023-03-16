using Plots, Printf
import Statistics:mean
include("./fcts_VE_trial.jl")
# MOTIVATED BY COMBINED RHEOLOGY:
# In previous version localiteration used one variable: ηpwl
# Adding creep mechanisms would result in adding new isolated viscosities (ηcst, ηlin...)
# To avoid this, I try to formulate the local iteration in terms of ηve_normal and ηve_shear
visu       = true
plasticity = true
showMD7    = true

function main() 

    # Parameters
    R         = 0*8.314
    εbg       = 1.0
    Δt        = .001
    nt        = 10
    My        = 1e6*365.25*24*3600
    PS        = 1
    
    npwl      = 3.0             # POWER LAW EXPONENT
    anifac_v  = 1.4             # VISCOUS anisotropy factor anifac_v = η_n/η_s
    anifac_e  = 4.0             # ELASTIC ----------------- anifac_e = G_n/G_s
    anifac_p  = 3.0             # PLASTIC ----------------- anifac_p = τxx_ve/τxy_ve
    a1 = 1.0; a2 = 1.0; a3 = 1.0 # friction angle anisotropy
    aniso_ang = LinRange( 0.0, π/2, 21 ) # Director angle (w.r. to horizontal), see Mühlhaus et al. 2002
    # aniso_ang = 0.0
    layer_ang = aniso_ang .- π/2 # Layer angle (w.r. to horizontal): in 2D the layers are orthogonal to the director
    
    # Numerics
    noisy   = 4
    tol     = 1.0e-12
    nitmax  = 20
    
    # Ambient conditions
    # Total strain rate
    εxx  = PS*εbg
    εyy  = -εxx*1.0 # Introduce a bit of compaction
    εzz  = -εxx-εyy
    εxy  = (1.0-PS)*εbg
    # Volumetric and deviatoric strain rates
    div  = εxx + εyy
    εxxd = εxx - 1.0/3.0*div
    εyyd = εyy - 1.0/3.0*div
    ezzd = -εxxd - εyyd
    P_ini = 0.0
    T     = 1.
    
    # Power law flow params
    τbg    = 1.0
    Bpwl   = 2^npwl*εbg/τbg^(npwl)
    τ_chk  = 2*Bpwl^(-1.0/npwl)*εbg^(1.0/npwl)
    ηpwl   = Bpwl^(-1.0/npwl)*εbg^(1.0/npwl-1)
    τ_chk  = 2*ηpwl*εbg
    if abs(τbg - τ_chk)/τbg > 1e-6 error("Power-law breaks down") end
    Cpwl   = (2*Bpwl^(-1/npwl))^(-npwl)
    G      = 1e1
    K      = 2.0e1
    ηe     = G*Δt
    C      = 0.5
    fric   = 0*30*π/180.0
    dil    = 0*10*π/180.0
    ηvp    = 0*1e18
    
    # Storage/visualisation
    τii_Ray    = zeros(nt, length(aniso_ang))
    τxx_ve_Ray = zeros(nt, length(aniso_ang))
    τyy_ve_Ray = zeros(nt, length(aniso_ang))
    τxy_ve_Ray = zeros(nt, length(aniso_ang))
    p_Ray      = zeros(nt, length(aniso_ang))
    τii_MD7    = zeros(nt, length(aniso_ang))
    τxx_ve_MD7 = zeros(nt, length(aniso_ang))
    τyy_ve_MD7 = zeros(nt, length(aniso_ang))
    τxy_ve_MD7 = zeros(nt, length(aniso_ang))
    p_MD7      = zeros(nt, length(aniso_ang))
    
    # Loop over all orientations
    for i in eachindex(aniso_ang)
    
        # Loop over time and model V-E-VP stress build-up 
        τxx_ve = 0.; τyy_ve = 0.; τxy_ve = 0.; P = P_ini;

        for it=1:nt
            if noisy>0 @printf("****** Time step %03d - layer angle = %02f ******\n", it, layer_ang[i]) end
            
            # Swap old deviatoric stresses
            τxx_ve0 = τxx_ve; τyy_ve0 = τyy_ve; τxy_ve0 = τxy_ve; P0 = P

            # Pressure update
            P = P0 - K*Δt*div
                
            ###############
    
            # Transformation matrix: towards principal plane
            Q            = [cos(layer_ang[i]) sin(layer_ang[i]); -sin(layer_ang[i]) cos(layer_ang[i])]
            @show ε_tens = [εxxd εxy; εxy εyyd]
            τ0_tens      = [τxx_ve0 τxy_ve0; τxy_ve0 τyy_ve0]
            @show ε_rot  = Q*ε_tens*Q'
            τ0_rot       = Q*τ0_tens*Q'
            @show cos(layer_ang[i])*sin(layer_ang[i])
            @show nx     = cos(aniso_ang[i])
            @show ny     = sin(aniso_ang[i])

            # ε        = [εxx; εyy; εxy]
            # τ0       = [τxx_ve0; τyy_ve0; τxy_ve0]
            # ε_eff    = ε .+ inv(Dani_e) * τ0 ./(2*ηe)

            ############### VISCO-ELASTIC TRIAL
            # Let's suppose all deformation is only visco-elastic (i.e no plastic deformation)
            τ_ve, ηve, sqrt_a_ve, ηpwl = LocalViscoElasticTrialStress_Newton_ηve( ε_rot, τ0_rot, ηe, sqrt(anifac_e), sqrt(anifac_v), Bpwl, Cpwl, npwl, tol, nitmax, noisy )
            anifac_ve = sqrt_a_ve^2
            τxx_ve    = τ_ve[1,1]; τxy_ve   = τ_ve[1,2]; τyy_ve   =  τ_ve[2,2]; τzz_ve   = -τxx_ve-τyy_ve
            J2 = 0.5*(τxx_ve^2 + τyy_ve^2 + τzz_ve^2) + τxy_ve^2*anifac_ve
            εxx_ve    = ε_rot[1,1]   + τ0_rot[1,1]/(2*ηe)
            εyy_ve    = ε_rot[2,2]   + τ0_rot[2,2]/(2*ηe)
            εxy_ve    = ε_rot[1,2]   + τ0_rot[1,2]/(2*ηe)*anifac_e

            εii_ve = sqrt( 1/2*(εxx_ve^2+εyy_ve^2) + εxy_ve^2/anifac_ve )
            @printf("VE:  Check J2 definition: %2.2e\n", (J2 - (2*ηve)^2 * ( εii_ve^2 ))/J2*100 )
            @printf("VE:  Check T2 definition: %2.2e\n", (sqrt(J2) - 2*ηve * εii_ve  ) / sqrt(J2)*100 )
            @printf("VE:  sqrt(J2) = %2.2e\n", sqrt(J2))
            @printf("VE:  ani_ve   = %2.2e\n", anifac_ve)
            # initialize vep variables (also if non-plastic case)
            ηvep        = ηve
            anifac_vep  = anifac_ve
            Pc          = P
            ############### PLASTIC CORRECTION
            if plasticity == true
                # Check yield condition
                J2      = 0.5*(τ_ve[1,1]^2 + τ_ve[2,2]^2) + τ_ve[1,2]^2*anifac_p
                J2t     = J2
                J2_corr = J2
                am      = (a1 + a2 + a3)/3
                F       = sqrt(J2) - C*cos(fric) - P*sin(fric)/3*(a1+a2+a3) +  sin(fric)*am
                @printf("Plasticity:  F   = %2.2e\n", F)
                # Initial guess 
                γ̇   = 0.
                Fc  = 0.
                if F>0 
                    eta_ve = ηve
                    τzz_ve = -τyy_ve-τxx_ve
                    τyx_ve = τxy_ve
                    am     = (a1 + a2 + a3)/3
                    Ji     = sin(fric)/3*( a1*τxx_ve + a2*τyy_ve + a3*τzz_ve) # WH: why devide by 3?
                    Tn2    = τxx_ve^2+τyy_ve^2+τzz_ve^2
                    Ts2    = (τxy_ve^2+τyx_ve^2)
                    Tii    = sqrt( 1/2*(Tn2 + anifac_p*Ts2)) # WH: compare with fcts_VE_trial.jl line 79 (Eii) -> not consistent?
                    for iter=1:nitmax
                        # It would be nice to have this section derived with Symbolics.jl
                        gdot     = γ̇
                        ap_n     = 1 -             eta_ve .* gdot ./  sqrt(J2)
                        ap_s_old = 1 - anifac_p .* eta_ve .* gdot ./ (sqrt(J2) .* anifac_ve ) # WH: not used, just for comparison
                        ap_s     = ap_n#1 - anifac_p .* eta_ve .* gdot ./ (sqrt(J2) .* anifac_ve )
                        Jii      = J2
                        dt       = Δt
                        eta_vp   = ηvp
                        
                        # With out of plane
                        Pc   = P + K*Δt*γ̇*sin(dil)*am

                        J1_corr   = ap_n*Ji
                        J2_corr   = 1/2*( Tn2*ap_n^2 + anifac_p*Ts2*ap_s^2)
                        Fc        =  sqrt(J2_corr) - C*cos(fric) - Pc.*am.*sin(fric) - eta_vp.*gdot;
                
                        dapndgdot = -eta_ve/Tii
                        dapsdgdot = dapndgdot#*ani_rat
                        dFdgdot   = -eta_vp + 1/2/sqrt(J2_corr) * ( anifac_p*Ts2*ap_s*dapsdgdot + Tn2*ap_n*dapndgdot ) - am^2*K*dt*sin(dil)*sin(fric) + dapndgdot*Ji
        
                        
                        dFdγ̇ = dFdgdot
                        # It would be nice to have this section derived with Symbolics.jl
                        if noisy>2 @printf("Newton VP It. %02d, r abs. = %2.2e r rel. = %2.2e\n", iter, Fc, Fc/F) end
                        if abs(Fc)/F < tol || abs(Fc)<tol
                            # error("sto")
                            @printf("Plasticity:  ap_n   = %2.2e\nPlasticity:  ap_s   = %2.2e\nPlasticity:  ap_s_old   = %2.2e\n", ap_n,ap_s,ap_s_old)
                            break
                        end
                        γ̇    -= Fc/dFdγ̇
                    end
                    # Post-process plasticity
                    if abs(Fc)/F>1e-6 Fc, error("Yield failed 3")  end
                    εxx_p      = γ̇*( (a1/3-a3/3)*sin(fric) + τxx_ve/sqrt(J2)/2)
                    εyy_p      = γ̇*( (a2/3-a3/3)*sin(fric) + τyy_ve/sqrt(J2)/2)
                    εxy_p      = γ̇*τxy_ve/sqrt(J2)/2*anifac_p
                    εyx_p      = γ̇*τxy_ve/sqrt(J2)/2*anifac_p
                    D          = 2*ηve * [1 0 0 0; 0 1 0 0; 0 0 1/anifac_ve 0; 0 0 0 1/anifac_ve]
                    τ_vec      = D*([εxx_ve; εyy_ve; εxy_ve; εxy_ve] .- [εxx_p; εyy_p; εxy_p; εyx_p])
                    τ_ve       = [τ_vec[1] τ_vec[3]; τ_vec[4] τ_vec[2]]
                    J2_corr_p  = 0.5*(τ_ve[1,1]^2 + τ_ve[2,2]^2) + τ_ve[1,2]^2*anifac_p
                    @printf("VEP: sqrt(J2_corr_p ) = %2.2e\n", sqrt(J2_corr_p) )
                    @printf("VEP: sqrt(J2_corr   ) = %2.2e\n", sqrt(J2_corr)   )
                    # Compute anifac_vep 
                    axx       = (sqrt(J2t)        -       ηve*γ̇)/(sqrt(J2t))
                    # axy       = (sqrt(J2t)*anifac_ve - anifac_p*ηve*γ̇)/(sqrt(J2t)*anifac_ve)
                    axy       = axx
                    η_mat     = ηve * [axx 0 0 0; 0 axx 0 0; 0 0 axy 0; 0 0 0 axy]
                    D         = 2*η_mat 
                    τ_vec     = D*([εxx_ve; εyy_ve; εxy_ve; εxy_ve]) 
                    # --- 
                    anifac_vep = anifac_ve*axx/axy
                    ηvep      = ηve*axx
                    η_mat     = ηvep * [1 0 0 0; 0 1 0 0; 0 0 1/anifac_vep 0; 0 0 0 1/anifac_vep] 
                    D         = 2*η_mat
                    τ_vec     = D*([εxx_ve; εyy_ve; εxy_ve; εxy_ve])
                    # ---
                    τ_ve      = [τ_vec[1] τ_vec[3]; τ_vec[4] τ_vec[2]]
                    τxx_vec      = τ_ve[1,1]; τyy_vec = τ_ve[2,2]; τzz_vec = -τxx_vec-τyy_vec
                    J2_corr_p = 0.5*(τ_ve[1,1]^2 + τ_ve[2,2]^2) + τ_ve[1,2]^2*anifac_p
                    J1_corr_p = a1*τxx_vec + a2*τyy_vec + a3*τzz_vec
                    Pc        = P + K*Δt*γ̇*sin(dil)/3*(a1+a2+a3)
                    Fchk      = sqrt(J2_corr_p) - C*cos(fric) - Pc*sin(fric)*am +  sin(fric)/3*J1_corr_p - ηvp*γ̇
                end
            end
    
            ##############################
            # Rotate back
            # τ   = τ_ve
            τ   = Q'*τ_ve*Q
            τii = sqrt( 0.5*(τ[1,1]^2 + τ[2,2]^2) + τ[1,2]^2)
            τxx_ve = τ[1,1]; τyy_ve = τ[2,2]; τxy_ve = τ[1,2]
            @printf("VEP: τii      = %2.6e\n", τii)
    
            ###############
            # Store
            τii_Ray[it,i]    = τii
            τxx_ve_Ray[it,i] = τxx_ve
            τyy_ve_Ray[it,i] = τyy_ve
            τxy_ve_Ray[it,i] = τxy_ve
            p_Ray[it,i]      = Pc
    
            # ###################################################################
            # THIS IS WHAT WILL HAPPEN IN MDOODZ: It still uses the C_ANI business we derived earlier...
            two      = 2.
            nx       = cos(aniso_ang[i]); ny = sin(aniso_ang[i])
            N        = [nx; ny]
            d0       = 2*N[1]^2*N[2]^2
            d1       = N[1]*N[2] * (-N[1]^2 + N[2]^2)
            C_ANI    = [-d0 d0 -two*d1; d0 -d0 two*d1; -d1 d1 -two*(1/2-d0)]      ###### !!! - sign in D33 -a0
            C_ISO    = [1 0 0; 0 1 0; 0 0 two*1//2]
            ani      = 1. - 1.  ./ anifac_vep
            Dani_vep = C_ISO .+ ani*C_ANI # Achtung: factor 2 removed from Dani
            ani      = 1. - 1.  ./ anifac_e
            Dani_e   = C_ISO .+ ani*C_ANI # Achtung: factor 2 removed from Dani
            ani      = 1. - 1.  ./ anifac_v
            Dani_v   = C_ISO .+ ani*C_ANI # Achtung: factor 2 removed from Dani
  
            # 1 ---------------------------------------------------------------
            ε        = [εxx; εyy; εxy]
            τ0       = [τxx_ve0; τyy_ve0; τxy_ve0]
            ε_eff    = ε .+ inv(Dani_e) * τ0 ./(2*ηe)
            τ_MD7_v1 = 2*ηvep*Dani_vep*ε_eff
         
            # 3 ---------------------------------------------------------------
            # same as 2) but without mistake. Doesn't look very sexy... Can maybe be simplified...
            E        = ηe*[1 0 0; 0 1 0; 0 0 1] .+ ηpwl*Dani_v*inv(Dani_e)
            τ_MD7_v2 = 2*ηvep*(ηpwl+ηe)*inv(E)*(Dani_v*ε + Dani_v*inv(Dani_e) * (τ0 ./(2*ηe)))
            
            # @printf("max. rel. VE tangent : %2.2e\n ", maximum( (τ_MD7_v1.-τ_MD7_v2)./τ_MD7_v1)*100 ) 
                
            τii_MD7[it,i]    = sqrt(0.5*(τ_MD7_v1[1]^2 + τ_MD7_v1[2]^2) + τ_MD7_v1[3]^2)
            @printf("VEP: τii      = %2.6e\n", τii_MD7[it,i])
            τxx_ve_MD7[it,i] = τ_MD7_v1[1]
            τyy_ve_MD7[it,i] = τ_MD7_v1[2]
            τxy_ve_MD7[it,i] = τ_MD7_v1[3]
            p_MD7[it,i]      = Pc
    
        end
    end
    
    if visu
        if length(aniso_ang)==1
        p=plot([1:nt].*(Δt/My), τii_Ray[:,1], linewidth=1, marker =:circle, color=:green,  label=:none)
        p=plot!([1:nt].*(Δt/My), τii_MD7[:,1], linewidth=0, marker =:xcross, color=:green,  label=:none)
        else
        st1 = Int64(floor(  nt/3))
        st2 = Int64(floor(2*nt/3))
        st3 = nt
        p1=plot( aniso_ang.*(180/π), τii_Ray[st1,:], color=:green, label="1/3 loading")
        p1=plot!(aniso_ang.*(180/π), τii_Ray[st2,:], color=:blue,  label="2/3 loading")
        p1=plot!(aniso_ang.*(180/π), τii_Ray[st3,:], color=:red,   label="3/3 loading")
        if showMD7
        p1=plot!(aniso_ang.*(180/π), τii_MD7[st1,:], marker =:circle, linewidth=0, color=:green, label=:none)
        p1=plot!(aniso_ang.*(180/π), τii_MD7[st2,:], marker =:circle, linewidth=0, color=:blue,  label=:none)
        p1=plot!(aniso_ang.*(180/π), τii_MD7[st3,:], marker =:circle, linewidth=0, color=:red,   label=:none)
        end
        p1=plot!(xlabel="n angle [°]", ylabel="τII [MPa]")
        
        i1 = 1  # angle 1
        i2 = Int64(floor(  length(aniso_ang)*1/2))  # angle 2
        i3 = Int64(floor(  length(aniso_ang)*2/3))  # angle 2
        p2=plot( [1:nt].*(Δt/My), τii_Ray[:,i1], color=:red,   label=@sprintf("n angle %2.2f°", aniso_ang[i1]*180/π) )
        p2=plot!([1:nt].*(Δt/My), τii_Ray[:,i2], color=:blue,  label=@sprintf("n angle %2.2f°", aniso_ang[i2]*180/π) )
        p2=plot!([1:nt].*(Δt/My), τii_Ray[:,i3], color=:green, label=@sprintf("n angle %2.2f°", aniso_ang[i3]*180/π) )
        if showMD7
        p2=plot!([1:nt].*(Δt/My), τii_MD7[:,i1], linewidth=0, marker =:circle, color=:red,    label=:none)
        p2=plot!([1:nt].*(Δt/My), τii_MD7[:,i2], linewidth=0, marker =:circle, color=:blue,   label=:none)
        p2=plot!([1:nt].*(Δt/My), τii_MD7[:,i3], linewidth=0, marker =:circle, color=:green,  label=:none)
        end
        p2=plot!(xlabel="t [My]", ylabel="τII [MPa]")

        i1 = 1  # angle 1
        i2 = Int64(floor(  length(aniso_ang)*1/2))  # angle 2
        i3 = Int64(floor(  length(aniso_ang)*2/3))  # angle 2
        p3=plot( [1:nt].*(Δt/My), p_Ray[:,i1], color=:red,   label=@sprintf("n angle %2.2f°", aniso_ang[i1]*180/π) )
        p3=plot!([1:nt].*(Δt/My), p_Ray[:,i2], color=:blue,  label=@sprintf("n angle %2.2f°", aniso_ang[i2]*180/π) )
        p3=plot!([1:nt].*(Δt/My), p_Ray[:,i3], color=:green, label=@sprintf("n angle %2.2f°", aniso_ang[i3]*180/π) )
        if showMD7
        p3=plot!([1:nt].*(Δt/My), p_MD7[:,i1], linewidth=0, marker =:circle, color=:red,    label=:none)
        p3=plot!([1:nt].*(Δt/My), p_MD7[:,i2], linewidth=0, marker =:circle, color=:blue,   label=:none)
        p3=plot!([1:nt].*(Δt/My), p_MD7[:,i3], linewidth=0, marker =:circle, color=:green,  label=:none)
        end
        p3=plot!(xlabel="t [My]", ylabel="p [MPa]")

        p4=plot(  aniso_ang.*(180/π), τxx_ve_Ray[st3,:], color=:blue,  label="τxx_ve Ray")
        p4=plot!( aniso_ang.*(180/π), τyy_ve_Ray[st3,:], color=:green, label="τyy_ve Ray")
        p4=plot!( aniso_ang.*(180/π), τxy_ve_Ray[st3,:], color=:red,   label="τxy_ve Ray")
        if showMD7
        p4=plot!( aniso_ang.*(180/π), τxx_ve_MD7[st3,:], linewidth=0, marker =:circle, color=:blue,  label="τxx_ve MD7")
        p4=plot!( aniso_ang.*(180/π), τyy_ve_MD7[st3,:], linewidth=0, marker =:circle, color=:green, label="τyy_ve MD7")
        p4=plot!( aniso_ang.*(180/π), τxy_ve_MD7[st3,:], linewidth=0, marker =:circle, color=:red,   label="τxy_ve MD7")
        end
        display(plot(p1,p2,p4,p3, size=(1200,900), legend=:none))
        end
    end
    
end

main()  
