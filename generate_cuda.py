from collections import OrderedDict

def generate_model_params_file(bioreactor, model_eqs, model_name, output_dir):
    model_eqs_dict = OrderedDict(sorted(model_eqs.items()))
    species = sorted(bioreactor.get_species_list())

    f = open(output_dir + "input_" + str(model_name) + ".xml", 'w')
    output_string = ""

    #  Generate input file
    output_string = output_string + "<name> " + "model_cuda_" + str(model_name) + " </name>\n"
    output_string = output_string + "<source> </source>\n"
    output_string = output_string + "<type> ODE </type>\n"
    output_string = output_string + "<fit>"

    for idx, s in enumerate(species):
        if "N_" in s:
            output_string = output_string + "species" + str(idx+1) + " "
    output_string = output_string + "</fit>\n"

    output_string = output_string + "<logp> False </logp>\n"
    output_string = output_string + "<initial>\n"

    prior_string = "<#PARAM#> #PRIOR# </#PARAM#>\n"

    for idx, eq in enumerate(model_eqs_dict):
        prior_val = ""
        if "A_" in species[idx]:
            prior_val = "constant 1e-10"

        elif "B_" in species[idx]:
            prior_val = "constant 0"

        elif "N_" in species[idx]:
            prior_val = "constant 1000000000000.0"

        elif "S_" in species[idx]:
            prior_val = "constant 4"

        else:
            print("error")
            exit()

        output_string = output_string + prior_string.replace("#PARAM#", eq).replace("#PRIOR#", prior_val)

    output_string = output_string + "</initial>\n\n"

    all_params = bioreactor.get_all_parameters()
    print(all_params)
    output_string = output_string + "<parameters>\n"

    for param_key in OrderedDict(sorted(all_params.items())):
        output_string = output_string + prior_string.replace("#PARAM#", param_key).replace("#PRIOR#", all_params[param_key])

    output_string = output_string + "</parameters>\n"

    f.write(output_string)
    return(output_string)

def generate_model_cuda(bioreactor, model_eqs, model_name, output_dir):
    model_eqs_dict = OrderedDict(sorted(model_eqs.items()))
    all_params = bioreactor.get_all_parameters()
    param_list = sorted(all_params)

    f = open(output_dir + "model_cuda_" + str(model_name) + ".cu", 'w')
    cuda_output = ''

    # Write includes
    cuda_output = cuda_output + "#include \"stdio.h\"\n"
    cuda_output = cuda_output + "#include <sstream>\n"
    cuda_output = cuda_output + "\n"

    # define number of species, params and react
    NSPECIES = len(model_eqs_dict)
    NPARAM = len(param_list)
    NREACT = 2

    cuda_output = cuda_output + "#define NSPECIES " + str(NSPECIES) + "\n"
    cuda_output = cuda_output + "#define NPARAM " + str(NPARAM) + "\n"
    cuda_output = cuda_output + "#define NREACT " + str(NREACT) + "\n"
    cuda_output = cuda_output + "\n"

    # Write parameter definitions
    parameter_string = "#define $PARAM$ tex2D(param_tex, #INDEX#, tid)\n"
    for idx, param in enumerate(param_list):
        write_param = parameter_string.replace('$PARAM$', param)
        write_param = write_param.replace('#INDEX#', str(idx))
        cuda_output = cuda_output + write_param

    cuda_output = cuda_output + "\n"

    # Begin myFex structure
    cuda_output = cuda_output + "struct myFex{\n__device__ void operator()(int *neq, double *t, double *y, double *ydot){\n\n"

    # Init TID
    cuda_output = cuda_output + "int tid = blockDim.x * blockIdx.x + threadIdx.x;\n\n"


    # Comment species order
    cuda_output = cuda_output + "// Order is: "
    for eq_key in model_eqs_dict:
        cuda_output = cuda_output + eq_key + ", "
    cuda_output = cuda_output + "\n\n"

    # Write equations
    equation_string = "ydot[#INDEX#] = #EQ#;\n"
    for idx, eq_key in enumerate(model_eqs_dict):
        print(eq_key)
        write_eq = equation_string.replace("#INDEX#", str(idx))
        write_eq = write_eq.replace("#EQ#", model_eqs_dict[eq_key])
        for idx, eq_key in enumerate(model_eqs_dict):
            write_eq = write_eq.replace(" " + eq_key + " ", " y[" + str(idx) + "] ")

        cuda_output = cuda_output + write_eq


    # Close myFex structure
    cuda_output = cuda_output + "}\n};"

    # Begin myJex
    cuda_output = cuda_output + "struct myJex{\n    __device__ void operator()(int *neq, double *t, double *y, int ml, " \
                                "int mu, double *pd, int nrowpd){\n"
    cuda_output = cuda_output + "return;\n"

    # Close myJex
    cuda_output = cuda_output + "} \n };"
    f.write(cuda_output)

def generate_input_file(model_list, input_file_output):
    output_string = ""
    # file_header
    output_string = output_string + "<input>\n"
    output_string = output_string + "<modelnumber> " + str(len(model_list)) + " </modelnumber>\n"
    output_string = output_string + "<restart> False </restart>\n"
    output_string = output_string + "<autoepsilon>\n"
    output_string = output_string + "<finalepsilon> 0.1 0.1 0.00001 </finalepsilon>\n"
    output_string = output_string + "<alpha> 0.7 </alpha>\n"
    output_string = output_string + "</autoepsilon>\n"
    output_string = output_string + "<particles> 100 </particles>\n"
    output_string = output_string + "<beta> 1 </beta>\n"
    output_string = output_string + "<dt> -1 </dt>\n"
    output_string = output_string + "<kernel> uniform </kernel>\n\n"

    times = range(0, 2500, 5)
    t_str = str(list(times)).replace(",", "").replace("[", "").replace("]", "")
    output_string = output_string + "<data>\n"
    output_string = output_string + "<times>\n"
    output_string = output_string + t_str + "\n"
    output_string = output_string + "</times>\n"
    output_string = output_string + "<variables> <var1> </var1> </variables>\n"
    output_string = output_string + "</data>\n\n"

    output_string = output_string + "<models>\n"
    for idx, model in enumerate(model_list):
        model_num = idx + 1
        output_string = output_string + "<model" + str(model_num) + ">\n"
        output_string = output_string + model
        output_string = output_string + "</model" + str(model_num) + ">\n\n"

    output_string = output_string + "</models>\n"
    output_string = output_string + "</input>\n"

    f = open(input_file_output, 'w')
    f.write(output_string)