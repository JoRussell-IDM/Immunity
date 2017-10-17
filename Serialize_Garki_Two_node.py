import os

from dtk.generic.climate import set_climate_constant
from dtk.interventions.input_EIR import add_InputEIR
from dtk.tools.serialization.serialization_tools import zero_infections
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment

init_flag = True
pickup_flag = False
analyze_flag = False

duration = 25550

exp_id = ''

# This block will be used unless overridden on the command-line
SetupParser.default_block = 'HPC'

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')
cb.update_params({'Demographics_Filenames': ['two_node_birth_cohort_demographics.json'],
                  'Vector_Species_Names': [],
                  'Simulation_Duration': 25550,
                  'Antigen_Switch_Rate': pow(10, -9.116590124),
                  'Base_Gametocyte_Mosquito_Survival_Rate': 0.002011099,
                  'Base_Gametocyte_Production_Rate': 0.06150582,
                  "Falciparum_MSP_Variants": 32,
                  "Falciparum_Nonspecific_Types": 76,
                  "Falciparum_PfEMP1_Variants": 1070,
                  "Gametocyte_Stage_Survival_Rate": 0.588569307,
                  "MSP1_Merozoite_Kill_Fraction": 0.511735322,
                  "Max_Individual_Infections": 3,
                  "Nonspecific_Antigenicity_Factor": 0.415111634,
                  'Serialization_Time_Steps': [25550],
                "Birth_Rate_Dependence": "FIXED_BIRTH_RATE"
                  #"Serialized_Population_Filenames": [
                  #    "state-00365.dtk"
                  #],
                  #"Serialized_Population_Path": r"\\\\internal.idm.ctr\\IDM\\Home\\jgerardin\\output\\Karen_ABV_newLS_Burnin_20170908_172301\\902\\bc1\\8fb\\902bc18f-ba94-e711-9401-f0921c16849d/output",


                  })
set_climate_constant(cb)

cb.set_input_files_root('Garki')

add_InputEIR(cb,[1.9375, 3.375, 7.5, 1.9375, 0.5, 0.5,0.5, 0.5, 0.25, 0.5, 0.5, 1.0],
             nodes={
            "Node_List": [
               1
            ],
            "class": "NodeSetNodeList"
         })
add_InputEIR(cb,[3.875, 7.75, 15.0, 3.875, 1, 1,1, 1, 0.5, 1, 1, 2 ],
             nodes={
            "Node_List": [
               2
            ],
            "class": "NodeSetNodeList"
         })

if __name__ == "__main__":
    SetupParser.init()
    init_manager = None
    if init_flag:
        init_manager = ExperimentManagerFactory.init()
        init_manager.run_simulations(exp_name='70_yr_burn_RafinMarke_fixed_births', config_builder=cb)
    # Wait for the simulations to be done
        init_manager.wait_for_finished(verbose=True)


