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
from malaria.reports.MalariaReport import add_survey_report, add_filtered_report
from malaria.reports.MalariaReport import add_summary_report
from dtk.interventions.health_seeking import add_health_seeking
from dtk.interventions.outbreakindividual import recurring_outbreak
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.utils.builders.sweep import GenericSweepBuilder
from dtk.vector.species import set_species_param
from dtk.generic.serialization import add_SerializationTimesteps
from dtk.generic.climate import set_climate_constant
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser

import os

if __name__ == "__main__":


        burnin_duration = 365 * 70
        sim_duration = 365*51
        num_seeds = 100

        exp_name  = 'Outbreak_healthseek_51yr_filtered'

        def add_my_outbreak(cb,start_day):
            recurring_outbreak(cb,start_day = start_day,repetitions=1)
            return {'outbreak_start_day': start_day}

        builder = ModBuilder.from_list(
            [[ModFn(add_my_outbreak, start_day=x*365),
              ModFn(DTKConfigBuilder.set_param, 'Run_Number', y),
              ModFn(DTKConfigBuilder.set_param, 'Simulation_Duration',(x*365)+(3*365)),
              ModFn(add_filtered_report,start = x*365,end = (x*365)+(3*365)),
              ModFn(add_health_seeking, start_day=0, targets=[
                {'trigger': 'NewClinicalCase', 'coverage': 1, 'agemin': 0, 'agemax': 200, 'seek': 1,'rate': 0.3},
                  {'trigger': 'NewSevereCase', 'coverage': 1, 'seek': 1, 'rate': 0.5}])

              ]

             for y in range(num_seeds)for x in range(51)]

            )
            #ModFn(add_my_outbreak, start_day=x),
            #ModFn(add_health_seeking, start_day=0, targets=[
                  #{'trigger': 'NewClinicalCase', 'coverage': 1, 'agemin': 0, 'agemax': 200, 'seek': 1,'rate': 0.3},
                  #{'trigger': 'NewSevereCase', 'coverage': 1, 'seek': 1, 'rate': 0.5}])

#              ]
        SetupParser.default_block = 'HPC'

        init_path = r'\\internal.idm.ctr\IDM\home\jorussell\output\Matsari_closed_loop__burnin_20171011_221421\439\9de\98d\4399de98-d1ae-e711-9414-f0921c16b9e5\output'
        state_name = 'state-{:05d}.dtk'.format(burnin_duration)
        state_cleared = 'cleared-state-{:05d}.dtk'.format(burnin_duration)

        zero_infections(os.path.join(init_path, state_name), os.path.join(init_path, state_cleared))


        cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')

        set_climate_constant(cb)

        set_species_param(cb, 'gambiae', 'Larval_Habitat_Types',
                          {"LINEAR_SPLINE": {
                              "Capacity_Distribution_Per_Year": {
                                  "Times": [0.0, 30.417, 60.833, 91.25, 121.667, 152.083,
                                            182.5, 212.917, 243.333, 273.75, 304.167, 334.583],
                                  # "Values": [1, 0.25, 0.1, 1, 1, 0.5, 12, 5, 3, 2, 1.5, 1.5]
                                  "Values": [3, 0.8, 1.25, 0.1, 2.7, 8, 4, 35, 6.8, 6.5, 2.6, 2.1]
                              },
                              "Max_Larval_Capacity": 1e9
                          }})
        set_species_param(cb, "gambiae", "Indoor_Feeding_Fraction", 1.0)

        cb.update_params({

        'Age_Initialization_Distribution_Type': 'DISTRIBUTION_SIMPLE',
        'Antigen_Switch_Rate': pow(10, -9.116590124),
        'Base_Gametocyte_Production_Rate': 0.06150582,
        'Base_Gametocyte_Mosquito_Survival_Rate': 0.002011099,
        'Base_Population_Scale_Factor': 1,  # Change x_Temporary_Larval_Habitat by same factor if changing Base_Population_Scale_Factor
        'x_Temporary_Larval_Habitat': 1,

        'Birth_Rate_Dependence': 'FIXED_BIRTH_RATE',
        "Death_Rate_Dependence": "NONDISEASE_MORTALITY_BY_AGE_AND_GENDER",

        "Demographics_Filenames": ['BitingRisk/single_node_demographics.json'],

        "Disable_IP_Whitelist": 1,

        'Enable_Birth': 1,
        'Enable_Default_Reporting': 1,
        # 'logLevel_default': 'ERROR',
        'Enable_Vital_Dynamics': 1,
        # 'Enable_Property_Output': 1,

        'Falciparum_MSP_Variants': 32,
        'Falciparum_Nonspecific_Types': 76,
        'Falciparum_PfEMP1_Variants': 1070,
        'Gametocyte_Stage_Survival_Rate': 0.588569307,

        'MSP1_Merozoite_Kill_Fraction': 0.511735322,
        'Max_Individual_Infections': 3,
        'Nonspecific_Antigenicity_Factor': 0.415111634,

        "Simulation_Duration": sim_duration + 1,
        "Run_Number": 0,
        "Vector_Species_Names": ['gambiae'],
        "Serialized_Population_Filenames": [
                                               state_cleared
                                           ],
        "Serialized_Population_Path": init_path

        })
        SetupParser.init()
        run_sim_args = {'config_builder': cb,
                        'exp_name': exp_name,
                        'exp_builder': builder}

        exp_manager = ExperimentManagerFactory.from_setup()
        exp_manager.run_simulations(**run_sim_args)
        exp_manager.wait_for_finished(verbose=True)

