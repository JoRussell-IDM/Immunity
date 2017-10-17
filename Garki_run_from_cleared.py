import os
import numpy as np

from dtk.generic.climate import set_climate_constant
from dtk.interventions.input_EIR import add_InputEIR
from dtk.tools.serialization.serialization_tools import zero_infections
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
from dtk.utils.builders.sweep import GenericSweepBuilder
from malaria.interventions.malaria_challenge import add_challenge_trial
from simtools.ModBuilder import  ModFn, ModBuilder
from malaria.reports.MalariaReport import add_survey_report
from malaria.reports.MalariaReport import add_summary_report

burnin_duration = 25550

monthlyEIRs_1 = [1.9375, 3.375, 7.5, 1.9375, 0.5, 0.5, 0.5, 0.5, 0.25, 0.5, 0.5,1.0]
monthlyEIRs_2 = [3.875, 7.75, 15.0, 3.875, 1, 1, 1, 1, 0.5, 1, 1,2]

exp_name  = 'Challenge_bite_20_350_Dan_fix'

builder = ModBuilder.from_list(
    [[ModFn(add_challenge_trial, x, disable_vitals = False),
      ModFn(add_survey_report,[x],reporting_interval = 730),
      ModFn(add_summary_report,x,interval = 365/12,description = 'Monthly Report',
             age_bins = [1.0, 4.0, 8.0, 18.0, 28.0, 43.0, 400000.0],
             parasitemia_bins = [0.0, 16.0, 409.0, 4000000.0])
     ]
      for x in range(730,740,10)])


# builder = ModBuilder.from_list(
#     [[ModFn(add_InputEIR, monthlyEIRs_1, start_day=x, nodes={
#             "Node_List": [
#                1
#             ],
#             "class": "NodeSetNodeList"
#          }),
#       ModFn(add_InputEIR, monthlyEIRs_2, start_day=x, nodes={
#             "Node_List": [
#                2
#             ],
#             "class": "NodeSetNodeList"
#          }),
#       ModFn(add_summary_report,x,interval = 365/12,description = 'Monthly Report',
#             age_bins = [1.0, 4.0, 8.0, 18.0, 28.0, 43.0, 400000.0],
#             parasitemia_bins = [0.0, 16.0, 409.0, 4000000.0])
#
#       ]
#
#             for x in range(730,740,10)]
#
#     )
# This block will be used unless overridden on the command-line
SetupParser.default_block = 'HPC'

init_path = r'\\internal.idm.ctr\IDM\home\jorussell\output\70_yr_burn_RafinMarke_fixed_births_20171005_214046\6de\475\d81\6de475d8-15aa-e711-9414-f0921c16b9e5\output'
state_name = 'state-{:05d}.dtk'.format(burnin_duration)
state_cleared = 'cleared-state-{:05d}.dtk'.format(burnin_duration)

# zero_infections(os.path.join(init_path, state_name), os.path.join(init_path, state_cleared))

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')
cb.update_params({'Demographics_Filenames': ['two_node_birth_cohort_demographics.json'],
                  'Vector_Species_Names': [],
                  'Simulation_Duration': 2000,
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

                  "Serialized_Population_Filenames": [
        state_cleared
    ],
                  "Serialized_Population_Path": init_path

                  })


set_climate_constant(cb)

cb.set_input_files_root('Garki')

run_sim_args = {
    'exp_name': exp_name,
    'exp_builder': builder,
    'config_builder':cb
}


if __name__ == "__main__":
    SetupParser.init()
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())




