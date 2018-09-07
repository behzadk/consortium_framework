import numpy as np
import xml.etree.ElementTree as ET
import build_equations
import generate_cuda

def gen_param_dict(species_name, params_spec):
    """
    Generates dictionary given a parameter specification xml file. Adds the species specific string to the end
    of the general param spec. Also uses the default values inside the param spec compatible with ABCsysbio.

    The '$' tags are added around the base parameter name so we can identify where it needs to be subbed into the
    general equation

    :param species_name: Name of species to append to the general param namses
    :param params_spec: Parameter specification file
    :return: A dictionary specific to the input, with the information needed to produced an input file for ABCsysbio
    """
    tree = ET.parse(params_spec)
    root = tree.getroot()
    param_dict = {}
    for child in root:
        if species_name == "":
            param_name = child.tag + species_name
            param_dict[param_name] = child.text

        else:
            param_name = child.tag + "_" +  species_name
            param_dict[param_name] = child.text

    return param_dict

class AHL_Expression:
    def __init__(self, AHL_name, params_spec):
        self.name = AHL_name
        self.params_dict = gen_param_dict(self.name, params_spec)

class Microcin_Expression:
    def __init__(self, microcin_name, params_spec, inducer, repressor):
        """
        Generates a microcin object containing parameter specs and defined inducers or repressors that refer to
        another species in the system by name.


        :param microcin_name:
        :param params_spec:
        :param inducer:
        :param repressor:
        """
        self.name = microcin_name
        self.params_dict = gen_param_dict(self.name, params_spec)
        self.inducer = inducer
        self.repressor = repressor

class Chassis:
    def __init__(self, strain_name, substrate_name, chassis_params_spec):
        self.name = strain_name
        self.params_dict = gen_param_dict(strain_name, chassis_params_spec)
        self.substrate = substrate_name


class Bioreactor:
    def __init__(self, bioreactor_params_spec):
        self.strains = []
        self.substrates = []

        self.strain_names = []
        self.substrate_names = []
        self.AHL_names = []
        self.microcin_names = []

        self.bioreactor_params = gen_param_dict("", bioreactor_params_spec)

    def add_new_strain(self, name, chassis,
                     AHL_expression,  microcin_expression,
                     microcin_sensitivity):

        new_strain = self.Strain(name, chassis,
                     AHL_expression,  microcin_expression,
                     microcin_sensitivity)

        self.strains.append(new_strain)
        self.strain_names.append(name)

        # Update the species lists for species that are inside the bioreactor
        for AHL in AHL_expression:
            AHL_name = AHL.name
            self.AHL_names.append(AHL_name)

        for microcin in microcin_expression:
            microcin_name = microcin.name
            self.microcin_names.append(microcin_name)

        self.substrate_names.append(chassis.substrate)

    def add_new_substrate(self, name, substrate_params_spec):
        new_substrate = self.Substrate(name, substrate_params_spec)
        self.substrates.append(new_substrate)

    def update_species_list(self):
        self.strain_names = sorted(list(set(self.strain_names)))
        self.substrate_names = sorted(list(set(self.substrate_names)))
        self.AHL_names = sorted(list(set(self.AHL_names)))
        self.microcin_names = sorted(list(set(self.microcin_names)))

    def get_species_list(self):
        strain_species =  list(map(lambda i: "N_" + i, self.strain_names))
        substrate_species = list(map(lambda i: "S_" + i, self.substrate_names))
        AHL_species = list(map(lambda i:  "A_" + i, self.AHL_names))
        microcin_species = list(map(lambda i:  "B_" + i, self.microcin_names))

        all_species = strain_species + substrate_species + AHL_species + microcin_species
        return all_species

    def list_species(self):
        print("Substrates: ")
        for s in self.substrate_names:
            print(s)
        print("")
        print("AHLs:")
        for a in self.AHL_names:
            print(a)
        print("\n")
        print("Microcins:")
        for m in self.microcin_names:
            print(m)


    def get_all_parameters(self):
        all_params = {}
        bioreactor_params = self.bioreactor_params
        all_params.update(bioreactor_params)

        for strain in self.strains:
            for AHL in strain.AHL_expression:
                all_params.update(AHL.params_dict)

            for microcin in strain.microcin_expression:
                all_params.update(microcin.params_dict)

            all_params.update(strain.chassis.params_dict)

        for substrate in self.substrates:
            all_params.update(substrate.params_dict)

        print(len(all_params.keys()))
        return all_params

    class Strain:
        def __init__(self, strain_name, chassis,
                     AHL_expression,  microcin_expression,
                     microcin_sensitivity):
            """
            Strains are initiated with a name, and arrays which contain objects that describe the expression
            of several species and describe strain properties.

            :param strain_name: name of string
            :param substrate_consumption: array containing substrate consumption object
            :param AHL_expression: array containing substrate consumption objects
            :param microcin_expression: array containing microcin expression objects
            :param microcin_sensitivity: array containing names of microcins this strain is sensitive to.
            """

            self.name = strain_name
            self.chassis = chassis
            self.AHL_expression = AHL_expression
            self.microcin_expression = microcin_expression
            self.microcin_sensitivity = microcin_sensitivity

        def list_parameters(self):
            print(self.chassis.params_dict)

            for AHL in self.AHL_expression:
                print(AHL.params_dict)

            for microcin in self.microcin_expression:
                print(microcin.params_dict)

    class Substrate:
        def __init__(self, substrate_name, substrate_params_spec):
            """
            Substrates are members of the bioreactor class because they are accessed by all strains in the system.

            :param substrate_name:
            :param substrate_params_spec:
            """
            self.name = substrate_name
            self.params_dict = gen_param_dict(substrate_name, substrate_params_spec)


def strain_generator(AHL_params_spec, microcin_params_spec):
    # 3 Strains
    # Each strain expresses maximum 1 AHL, maximum 1 microcin
    # Sensitive always to one microcin or none

    AHL_names = ["1", "2"]
    microcin_names = ["mccX", "mccY"]
    strain_names = ["i", "o", "p"]


    # Generate AHL objects
    AHL_objects = []
    for name in AHL_names:
        AHL_objects.append([AHL_Expression(name, AHL_params_spec)])

    AHL_objects.append([])


    microcin_objects = []
    # Generate microcin objects driven by different AHLs
    for A_idx, a_name in enumerate(AHL_names):
        for M_idx, m_name in enumerate(microcin_names):
            new_induced_microcin = Microcin_Expression(m_name, microcin_params_spec, a_name, None)
            new_repressed_microcin = Microcin_Expression(m_name, microcin_params_spec, None, a_name)
            microcin_objects.extend( ([new_induced_microcin], [new_repressed_microcin]) )

    microcin_objects.append([])

    microcin_sensitivity = []
    for name in microcin_names:
        microcin_sensitivity.append([name])

    microcin_sensitivity.append([])

    part_combinations = []
    for a in AHL_objects:
        for m in microcin_objects:
            for s in microcin_sensitivity:
                part_combinations.append([a, m, s])

    print(len(part_combinations))
    strain_combinations = []

    for idx_1, n_1 in enumerate(part_combinations):
        for idx_2, n_2 in enumerate(part_combinations):
            if idx_2 == idx_1:
                continue

            chassis_1 = Chassis("i", "glu", chassis_params_spec)
            chassis_2 = Chassis("o", "glu", chassis_params_spec)

            n_1 = [chassis_1] + n_1
            n_2 = [chassis_2] + n_2
            new_strain_set = [n_1, n_2]
            strain_combinations.append(1)
            print(len(strain_combinations))

    np.shape(strain_combinations)


def know_two_species_gen(bioreactor_params_spec, substrate_params_spec, AHL_params_spec, microcin_params_spec, chassis_params_spec):
    chassis_1 = Chassis("x", "glu", chassis_params_spec)
    chassis_2 = Chassis("c", "glu", chassis_params_spec)

    bioreactor_a = Bioreactor(bioreactor_params_spec)
    bioreactor_b = Bioreactor(bioreactor_params_spec)
    bioreactor_c = Bioreactor(bioreactor_params_spec)

    bioreactor_a.add_new_substrate("glu", substrate_params_spec)
    bioreactor_b.add_new_substrate("glu", substrate_params_spec)
    bioreactor_c.add_new_substrate("glu", substrate_params_spec)


    AHL_names = ["a", "b"]
    microcin_names = ["y", "z"]
    strain_names = ["x", "c"]

    AHL_a = AHL_Expression(AHL_names[0], AHL_params_spec)
    AHL_b = AHL_Expression(AHL_names[1], AHL_params_spec)

    # Model A
    x_microcin_exp = Microcin_Expression("y", microcin_params_spec, "a", None)
    c_microcin_exp = Microcin_Expression("z", microcin_params_spec, "b", None)

    x_sensitivity = ["y"]
    c_sensitivity = ["z"]
    bioreactor_a.add_new_strain("x", chassis_1, [AHL_a], [x_microcin_exp], x_sensitivity)
    bioreactor_a.add_new_strain("c", chassis_2, [AHL_b], [c_microcin_exp], c_sensitivity)

    # Model B
    x_microcin_exp = Microcin_Expression("y", microcin_params_spec, None, "b")
    c_microcin_exp = Microcin_Expression("z", microcin_params_spec, "a", None)

    x_sensitivity = ["y"]
    c_sensitivity = ["z"]

    bioreactor_b.add_new_strain("x", chassis_1, [AHL_a], [x_microcin_exp], x_sensitivity)
    bioreactor_b.add_new_strain("c", chassis_2, [AHL_b], [c_microcin_exp], c_sensitivity)


    # Model_C
    x_microcin_exp = Microcin_Expression("y", microcin_params_spec, None, "a")

    x_sensitivity = []
    c_sensitivity = ["y"]

    bioreactor_c.add_new_strain("x", chassis_1, [AHL_a], [x_microcin_exp], x_sensitivity)
    bioreactor_c.add_new_strain("c", chassis_2, [], [], c_sensitivity)

    # Model d
    bioreactor_d = Bioreactor(bioreactor_params_spec)
    bioreactor_d.add_new_substrate("glu", substrate_params_spec)

    x_microcin_exp = Microcin_Expression("y", microcin_params_spec, "a", None)

    x_sensitivity = []
    c_sensitivity = ["y"]

    bioreactor_d.add_new_strain("x", chassis_1, [], [x_microcin_exp], x_sensitivity)
    bioreactor_d.add_new_strain("c", chassis_2, [AHL_a], [], c_sensitivity)



    bioreactor_a.update_species_list()
    bioreactor_b.update_species_list()
    bioreactor_c.update_species_list()
    bioreactor_d.update_species_list()

    model_1_eqs = build_equations.build_bioreactor_equations(bioreactor_a)
    model_1_params = generate_cuda.generate_model_params_file(bioreactor_a, model_1_eqs, "model_1", "./new_models/")
    generate_cuda.generate_model_cuda(bioreactor_a, model_1_eqs, "model_1", "./new_models/")

    model_2_eqs = build_equations.build_bioreactor_equations(bioreactor_b)
    model_2_params = generate_cuda.generate_model_params_file(bioreactor_b, model_2_eqs, "model_2", "./new_models/")
    generate_cuda.generate_model_cuda(bioreactor_b, model_2_eqs, "model_2", "./new_models/")

    model_3_eqs = build_equations.build_bioreactor_equations(bioreactor_c)
    model_3_params = generate_cuda.generate_model_params_file(bioreactor_c, model_3_eqs, "model_3", "./new_models/")
    generate_cuda.generate_model_cuda(bioreactor_c, model_3_eqs, "model_3", "./new_models/")

    model_4_eqs = build_equations.build_bioreactor_equations(bioreactor_d)
    model_4_params = generate_cuda.generate_model_params_file(bioreactor_d, model_4_eqs, "model_4", "./new_models/")
    generate_cuda.generate_model_cuda(bioreactor_d, model_4_eqs, "model_4", "./new_models/")

    generate_cuda.generate_input_file([model_1_params, model_2_params, model_3_params, model_4_params], "./new_models/input_multi_model.xml")
    generate_cuda.generate_input_file([model_4_params], "./new_models/input_file_model_4.xml")

def repressilator(bioreactor_params_spec, substrate_params_spec, AHL_params_spec, microcin_params_spec, chassis_params_spec):
    # Initiate three chassis required
    chassis_1 = Chassis("x", "glu", chassis_params_spec)
    chassis_2 = Chassis("c", "glu", chassis_params_spec)
    chassis_3 = Chassis("f", "glu", chassis_params_spec)


    # Iniitiate bioreactor and envirnoment
    bioreactor_a = Bioreactor(bioreactor_params_spec)
    bioreactor_a.add_new_substrate('glu', substrate_params_spec)


    # Three constitutively expressed microcins
    microcin_x_exp = Microcin_Expression("x", microcin_params_spec, None, None)
    microcin_c_exp = Microcin_Expression("c", microcin_params_spec, None, None)
    microcin_f_exp = Microcin_Expression("f", microcin_params_spec, None, None)

    x_sensitivity = ["x"]
    c_sensitivity = ["c"]
    f_sensitivity = ["f"]

    bioreactor_a.add_new_strain("x", chassis_1, [], [microcin_x_exp], x_sensitivity)
    bioreactor_a.add_new_strain("c", chassis_2, [], [microcin_c_exp], c_sensitivity)
    bioreactor_a.add_new_strain("f", chassis_3, [], [microcin_f_exp], f_sensitivity)

    population_rpr_eqs = build_equations.build_bioreactor_equations(bioreactor_a)
    population_rpr_params = generate_cuda.generate_model_params_file(bioreactor_a, population_rpr_eqs, "ring_rpr_model_1", "./new_models/")
    generate_cuda.generate_model_cuda(bioreactor_a, population_rpr_eqs, "ring_rpr_model_1", "./new_models/")

if __name__ == "__main__":
    # Paths to parameter spec files
    AHL_params_spec = "./species_param_specs/AHL_params.xml"
    microcin_params_spec = "./species_param_specs/microcin_params.xml"
    chassis_params_spec = "./species_param_specs/chassis_params.xml"
    bioreactor_params_spec = "./species_param_specs/bioreactor_params.xml"
    substrate_params_spec = "./species_param_specs/substrate_params.xml"

    # know_two_species_gen(bioreactor_params_spec, substrate_params_spec, AHL_params_spec, microcin_params_spec, chassis_params_spec)
    repressilator(bioreactor_params_spec, substrate_params_spec, AHL_params_spec, microcin_params_spec, chassis_params_spec)

    exit()
    # Create objects which essentially represent parts.
    AHL_1 = AHL_Expression("1", AHL_params_spec)

    chassis_1 = Chassis("CV514", "leu", chassis_params_spec)
    chassis_2 = Chassis("MG", "glu", chassis_params_spec)
    chassis_3 = Chassis("XX", "glu", chassis_params_spec)

    microcin_1 = Microcin_Expression("1", microcin_params_spec, "1", None)

    # Create arrays for initiating a strain N1
    N1_AHL = [AHL_1]
    N1_chassis = chassis_1
    N1_microcin_expr = [microcin_1]
    N1_microcin_sens = ["1"]
    N1_substrate_consumption = ["leu"]

    N2_AHL = [AHL_1]
    N2_chassis = chassis_2
    N2_microcin_expr = [microcin_1]
    N2_microcin_sens = []
    N2_substrate_consumption = ["glu"]

    N3_AHL = []
    N3_chassis = chassis_3
    N3_microcin_expr = []
    N3_microcin_sens = []
    N3_substrate_consumption = ["glu"]


    # Initiate new bioreactor and load up a strain and a substrate
    b = Bioreactor(bioreactor_params_spec)
    b.add_new_strain("CV514", N1_chassis, N1_AHL, N1_microcin_expr, N1_microcin_sens)
    b.add_new_substrate("leu", substrate_params_spec)

    b.add_new_strain("MG", N2_chassis, N2_AHL, N2_microcin_expr, N2_microcin_sens)
    b.add_new_strain("XX", N3_chassis, N3_AHL, N3_microcin_expr, N3_microcin_sens)

    b.add_new_substrate("glu", substrate_params_spec)

    b.update_species_list()
    # b.list_species()
    # b.strains[0].list_parameters()
    model_1_eqs = build_equations.build_bioreactor_equations(b)
    model_2_eqs = build_equations.build_bioreactor_equations(b)

    model_1_params = generate_cuda.generate_model_params_file(b, model_1_eqs, 1, "./new_models/")
    model_2_params = generate_cuda.generate_model_params_file(b, model_2_eqs, 2, "./new_models/")

    generate_cuda.generate_model_cuda(b, model_1_eqs, 1, "./new_models/")
    generate_cuda.generate_model_cuda(b, model_2_eqs, 2, "./new_models/")

    generate_cuda.generate_input_file([model_1_params, model_2_params], "./new_models/input_multi_model.xml")