import numpy as np
import xml.etree.ElementTree as ET



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
        param_name = '$' + child.tag + '$' + "_" +  species_name
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
        :param constitutive:
        """
        self.name = microcin_name
        self.params_dict = gen_param_dict(self.name, params_spec)
        self.inducer = inducer
        self.repressor = repressor
        

class Bioreactor:
    def __init__(self):
        self.strains = []

    def add_new_strain(self, name, substrate_consumption,
                     AHL_expression,  microcin_expression,
                     microcin_sensitivity):
        s = self.Strain(name, substrate_consumption,
                     AHL_expression,  microcin_expression,
                     microcin_sensitivity)
        self.strains.append(s)

    def iterate_strains(self):
        for s in self.strains:
            print(s.name)

    class Strain:
        def __init__(self, strain_name, substrate_consumption,
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





            self.substrate_consumption = substrate_consumption
            self.microcin_sensitivities = microcin_sensitivity

            self.AHL_expression = AHL_expression
            self.microcin_expression = microcin_expression




if __name__ == "__main__":
    AHL_params_spec = "/species_param_specs/AHL_params.xml"
    microcin_params_spec = "/species_param_specs/microcin_params.xml"

    # AHL_names = ["1", "2", "3"]
    # microcin_names = ["mccI", "mccV"]
    # substrate_names = ["1", "2", "3"]

    AHL_1 = AHL_Expression("AHL_1", AHL_params_spec)
    microcin_1 = Microcin_Expression("mccI", microcin_params_spec)

    N1_AHL = [AHL_1]
    N1_microcin_expr = [microcin_1]
    N1_microcin_sens = ["mccI"]
    N1_substrate_consumption = ["1"]

    b = Bioreactor()
    b.add_new_strain("1", N1_substrate_consumption, N1_AHL, N1_microcin_expr, N1_microcin_sens)

    b.iterate_strains()