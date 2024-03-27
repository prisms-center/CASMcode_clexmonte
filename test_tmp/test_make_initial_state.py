import json
import numpy as np
import libcasm.clexmonte as clexmonte
import libcasm.clexmonte.semigrand_canonical as sgc
import libcasm.xtal as xtal

with open("data/Clex_ZrO_Occ/system.json", "r") as f:
    system_data = json.load(f)

system = clexmonte.System.from_dict(
    data=system_data,
    search_path=["data/Clex_ZrO_Occ"],
)

###
import libcasm.configuration as casmconfig
import libcasm.enumerate as casmenum

supercells = casmconfig.SupercellSet(prim=system.prim)
configurations = casmconfig.ConfigurationSet()
enum = casmenum.ConfigEnumAllOccupations(prim=system.prim, supercell_set=supercells)
for config in enum.by_supercell(
    supercells={"max": 6},
    skip_non_canonical=True,
    skip_non_primitive=True,
):
    configurations.add(config)
###
print(len(configurations))
# assert len(configurations) == 336  # max=4

f = clexmonte.MonteCarloConfiguration.from_config_with_sym_info


def mol_composition(config):
    return system.composition_calculator.mean_num_each_component(config.occupation)


def param_composition(config):
    return system.composition_converter.param_composition(mol_composition(config))


#
# for record in configurations:
#     configuration = record.configuration
#     name = record.configuration_name
#
#     mc_state = clexmonte.MonteCarloState(
#         configuration=f(configuration),
#     )
#     formation_energy = system.clex(mc_state, "formation_energy")
#     value = formation_energy.per_unitcell()
#
#     comp_n = mol_composition(configuration)
#     comp = param_composition(configuration)
#     print(f"- name={name}, ef={value}, comp_n={comp_n}, comp={comp}")
#
# exit()


motif_list = []
id_list = []

y = np.arange(-4, 0.01, step=0.1)
for x in y:
    initial_state, motif, id = sgc.make_initial_state(
        system=system,
        conditions={
            "temperature": 300.0,
            "param_chem_pot": [x],
        },
        dirs="abc",
        min_volume=1000,
        transformation_matrix_to_super=np.array(
            [
                [1, 0, 0],
                [0, 3, 0],
                [0, 0, 2],
            ],
            dtype="int",
        ),
        motif=None,
        configurations=configurations,
    )
    print(f"- x={x}, id={id}")
    if id not in id_list:
        print(f"motif: (id={id})")
        print(xtal.pretty_json(motif.to_dict()))
        print(xtal.pretty_json(motif.to_structure().to_dict()))
        motif_list.append(motif)
        id_list.append(id)

assert isinstance(initial_state, clexmonte.MonteCarloState)
print(initial_state.configuration.transformation_matrix_to_super)
assert False
