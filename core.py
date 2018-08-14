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


def part_generator():
    AHL_1 = AHL_Expression("1", AHL_params_spec)
    AHL_2 = AHL_Expression("2", AHL_params_spec)
    AHL_3 = AHL_Expression("3", AHL_params_spec)

    microcin_1 = Microcin_Expression("1", microcin_params_spec, "AHL1", None)
    microcin_2 = Microcin_Expression("1", microcin_params_spec, "AHL1", None)
    microcin_3 = Microcin_Expression("1", microcin_params_spec, "AHL1", None)

if __name__ == "__main__":
    # Paths to parameter spec files
    AHL_params_spec = "./species_param_specs/AHL_params.xml"
    microcin_params_spec = "./species_param_specs/microcin_params.xml"
    chassis_params_spec = "./species_param_specs/chassis_params.xml"
    bioreactor_params_spec = "./species_param_specs/bioreactor_params.xml"
    substrate_params_spec = "./species_param_specs/substrate_params.xml"

    strain_names = ["A", "B", "C"]
    AHL_names = ["1", "2", "3"]
    microcin_names = ["1", "2", "3"]



    # Create objects which essentially represent parts.
    chassis_1 = Chassis("CV514", "leu", chassis_params_spec)
    chassis_2 = Chassis("MG", "glu", chassis_params_spec)

    microcin_1 = Microcin_Expression("1", microcin_params_spec, "AHL1", None)

    # Create arrays for initiating a strain N1
    N1_AHL = [AHL_1]
    N1_chassis = chassis_1
    N1_microcin_expr = [microcin_1]
    N1_microcin_sens = ["mccI"]
    N1_substrate_consumption = ["leu"]

    N2_AHL = [AHL_1]
    N2_chassis = chassis_2
    N2_microcin_expr = [microcin_1]
    N2_microcin_sens = []
    N2_substrate_consumption = ["glu"]

    N3_AHL = []
    N3_chassis = chassis_2
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
    model_eqs = build_equations.build_bioreactor_equations(b)
    generate_cuda.generate_model_params_file(b, model_eqs, 1, "./new_models/")
    generate_cuda.generate_model_cuda(b, model_eqs, 1, "./new_models/")