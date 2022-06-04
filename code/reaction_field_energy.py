def reaction_field_energy_DNA(reaction_field_energy_file):
    reaction_field_energy=open(reaction_field_energy_file).readlines()[0].strip('\n').split()
    DNA=float(reaction_field_energy[2])
    return DNA

def change_reaction_field_energy(reaction_field_energy_file):
    reaction_field_energy=open(reaction_field_energy_file).readlines()[0].strip('\n').split()
    complex=float(reaction_field_energy[0])
    protein=float(reaction_field_energy[1])
    DNA=float(reaction_field_energy[2])
    change=protein+DNA-complex
    return change