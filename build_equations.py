def gen_strain_growth_eq(strain):
    general_mu_strain = "( mu_max_#STRAIN_NAME# * S_#SUB_NAME# / ( K_#STRAIN_NAME# + S_#SUB_NAME# ) ) "
    general_sensitivity_term = " - omega_#MICROCIN_NAME# * B_#MICROCIN_NAME# "

    strain_sensitivty = ""

    # Build strain specific sensitivity term
    for microcin in strain.microcin_sensitivity:
        strain_sensitivty = strain_sensitivty + general_sensitivity_term.replace("#MICROCIN_NAME#", microcin)

    dN =  "( " + general_mu_strain + strain_sensitivty +  " - D) * N_#STRAIN_NAME# "
    dN = dN.replace("#STRAIN_NAME#", strain.name)
    dN = dN.replace("#SUB_NAME#", strain.chassis.substrate)

    return dN

def gen_substrate_eq(substrate_name, bioreactor):
    # Init string with dilution term. Sub in the name of the substrate
    dS = (" D * (S0_#SUB_NAME# - S_#SUB_NAME# ) ")

    general_mu_strain = "( mu_max_#STRAIN_NAME# * S_#SUB_NAME# / ( K_#STRAIN_NAME# + S_#SUB_NAME# ) ) "
    general_S_yield = " * N_#STRAIN_NAME# / g_#STRAIN_NAME# )"

    # Iterate through bioreactor strains, adding the yield term for each consumer
    for strain in bioreactor.strains:
        if strain.chassis.substrate == substrate_name:
            # Consumption term
            strain_name = strain.name
            strain_mu = general_mu_strain.replace("#STRAIN_NAME#", strain_name)
            strain_yield = general_S_yield.replace("#STRAIN_NAME#", strain_name)
            dS = dS + "- ( " + strain_mu + strain_yield
            dS = dS.replace("#SUB_NAME#", substrate_name)

    return dS


def gen_AHL_eq(AHL_name, bioreactor):
    # AHL production general term
    general_AHL_production = (" + kA_#AHL_NAME# * N_#STRAIN_NAME# ").replace("#AHL_NAME#", AHL_name)

    # Create empty string to build on
    dAHL = ""

    for strain in bioreactor.strains:
        for AHL in strain.AHL_expression:
            if AHL.name == AHL_name:
                dAHL = dAHL + general_AHL_production.replace("#STRAIN_NAME#", strain.name)

    dAHL = dAHL + " - D * A_" + AHL_name + " "

    return(dAHL)

def gen_microcin_eq(microcin_name, bioreactor):
    gen_induced_microcin = " * ( powf( A_#AHL_NAME# , nB_#MICROCIN_NAME# ) / " \
                           "( powf( KB_#MICROCIN_NAME#, nB_#MICROCIN_NAME# ) + powf( A_#AHL_NAME# , nB_#MICROCIN_NAME# ) ) )"

    gen_repressed_microcin = " * ( powf( KB_#MICROCIN_NAME# , nB_#MICROCIN_NAME# ) / " \
                             "( powf( KB_#MICROCIN_NAME# , nB_#MICROCIN_NAME# ) + powf( A_#AHL_NAME# , nB_#MICROCIN_NAME# ) ) )"

    gen_induced_microcin = gen_induced_microcin.replace("#MICROCIN_NAME#", microcin_name)
    gen_repressed_microcin = gen_repressed_microcin.replace("#MICROCIN_NAME#", microcin_name)

    dB = ""

    for strain in bioreactor.strains:
        strain_production = " + kBmax_#MICROCIN_NAME# ".replace("#MICROCIN_NAME#", microcin_name)
        for microcin in strain.microcin_expression:
            if microcin.name == microcin_name:
                # Check for induced or repressed. It is possible for microcin to be both induced and repressed simultaneously.
                inducer = microcin.inducer
                repressor = microcin.repressor

                if inducer is not None:
                    strain_production = strain_production + gen_induced_microcin.replace("#AHL_NAME#", inducer)

                if repressor is not None:
                    strain_production = strain_production + gen_repressed_microcin.replace("#AHL_NAME#", repressor)


            dB = dB + strain_production + " * N_#STRAIN_NAME# "
            dB = dB.replace("#STRAIN_NAME#", strain.name)

    dB = dB + " - D * B_" + microcin_name + " "
    return dB


def build_bioreactor_equations(bioreactor):
    model_eqs = {}

    # build substrate equations for each substrate present in the bioreactor.
    for n in bioreactor.strains:
        n_name = "N_" + n.name
        model_eqs[n_name] = gen_strain_growth_eq(n)

    for s in bioreactor.substrate_names:
        s_name = "S_" + s
        model_eqs[s_name] = gen_substrate_eq(s, bioreactor)

    for a in bioreactor.AHL_names:
        a_name = "A_" + a
        model_eqs[a_name] = gen_AHL_eq(a, bioreactor)

    for m in bioreactor.microcin_names:
        m_name = "B_" + m
        model_eqs[m_name] = gen_microcin_eq(m, bioreactor)

    return model_eqs