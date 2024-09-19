from pymol import cmd
from pymol.cmd import util
import sys
import os

sys.path.append(os.getcwd())
from get_raw_distances import get_raw_distances

def get_residue_info(selection):
    model = cmd.get_model(selection)
    if model.atom:
        atom = model.atom[0]
        return f"{atom.resn}{atom.resi}"
    return None

def parse_dist(distances):
    donor_acceptor_pairs = []
    for dist in distances:
        donor_index = dist[0]  # Index of the donor atom
        acceptor_index = dist[1]  # Index of the acceptor atom
        #
        donor_selection = f"{donor_index[0]} and index {donor_index[1]}"
        acceptor_selection = f"{acceptor_index[0]} and index {acceptor_index[1]}"
        #
        donor_residue = get_residue_info(donor_selection)
        acceptor_residue = get_residue_info(acceptor_selection)
        #
        if donor_residue and acceptor_residue:
            donor_acceptor_pairs.append((donor_residue, acceptor_residue))
    return donor_acceptor_pairs

def print_tuple(name, pairs):
    print(f"\n{name}, {len(pairs)}")
    print(f"donor\tacceptor")
    for p in pairs:
        print(f"{p[0]}\t{p[1]}")

def HBond(ChA,ChB):
    name = f'Hbond'
    dist = cmd.distance(name, ChA,ChB, '3.5', '2')
    #
    distances = get_raw_distances(name)
    do_ac = parse_dist(distances)
    print_tuple(name, do_ac)
    return name

def Pi_Pi_interaction(ChA,ChB):
    name = f'PI_PI'
    cmd.distance(name, ChA, ChB,'5.5','6')
    return name

def Cation_Pi_interaction(ChA,ChB):
    name = f'PI_Cation'
    cmd.distance(name, ChA, ChB, '5.5', '7')
    return name

def Salt_bridge_positive(ChA,ChB):
    name = f'SBP'
    selection_1 = f'{ChB} and ((resn LYS and name NZ) or (resn ARG and name NE+NH*))'
    selection_2 = f'{ChA} and resn ASP+GLU and name OD*+OE*'
    cmd.distance(name,selection_1,selection_2, '5.0','0')
    return name

def Salt_bridge_negative(ChA,ChB):
    name = f'SBN'
    selection_1 = f'{ChA} and ((resn LYS and name NZ) or (resn ARG and name NE+NH*))'
    selection_2 = f'{ChB} and resn ASP+GLU and name OD*+OE*'
    cmd.distance(name,selection_1,selection_2,'5.0','0')
    return name

def run(pdb,ChA,ChB,dist):
    cmd.delete('all')
    cmd.load(pdb)
    #pdb_name = pdb[:-4]    # when af3
    pdb_name = 'ranked_0'    # when TCRmodel2, G0503_final_models/ranked_0.pdb
    print(pdb_name)
    cmd.remove('not polymer')
    cmd.select('inter_1',f'{pdb_name} and chain {ChA} and byres(chain {ChB}) around {dist}')
    cmd.select('inter_2', f'{pdb_name} and chain {ChB} and byres(chain {ChA}) around {dist}')
    cmd.show('line','inter_1')     # 周围残基设置为stick/line
    cmd.show('stick', 'inter_2')    # 周围残基设置为stick/line
    cmd.label(f'name CA and (inter_1 or inter_2)', 'oneletter+resi')
    cmd.set('label_bg_color', -7, '', 0)
    cmd.set('label_size', 18, '', 0)
    util.cba(33, f'chain {ChA}')

    if ChB in ['H+L','L+H']:
        util.cba(156, 'chain H')
        util.cba(20, 'chain L')
    else:
        util.cba(20, f'chain {ChB}')

    cmd.set('stick_radius', '0.1', 'all')
    cmd.set('cartoon_transparency', 0.7)
    # interaction
    Hbond = HBond('inter_1', 'inter_2')                 # cutoff=3.5 & mode=2 (polar contact)
    PI = Pi_Pi_interaction('inter_1', 'inter_2')        # cutoff=5.5 & mode=6 (pi-pi)
    C_PI = Cation_Pi_interaction('inter_1', 'inter_2')  # cutoff=5.5 & mode=7 (pi-cation)
    SBP = Salt_bridge_positive('inter_1', 'inter_2')    # B(LYS+ARG)-A(ASP+GLU) & cutoff=5.0 & mode=0 (all)
    SBN = Salt_bridge_negative('inter_1', 'inter_2')    # B(ASP+GLU)-A(LYS+ARG) & cutoff=5.0 & mode=0 (all)
    # color interaction
    cmd.set('dash_color', 'yellow', Hbond)
    cmd.set('dash_color', 'blue', PI)
    cmd.set('dash_color', 'blue', C_PI)
    cmd.set('dash_color', 'orange', SBP)
    cmd.set('dash_color', 'orange', SBN)
    cmd.set('dash_radius', '0.07')
    cmd.set('dash_gap', '0.3')
    cmd.remove('e. h and  neighbor(name C*)')
    cmd.zoom('Hbond')
    cmd.save(f'{pdb_name}_interaction.pse')

def main():
    pdb_file = sys.argv[1]
    CHA = sys.argv[2]
    CHB = sys.argv[3]
    dist = sys.argv[4]
    cmd.set('bg_rgb','white')
    run(pdb_file,CHA,CHB,dist)

main()
