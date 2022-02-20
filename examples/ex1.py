from materia.structure import Structure

s = Structure.generate("smiles", "O")
print(s.center_of_mass)
print(s.multiplicity)
print(s.name)
print(s.fragments)
s.fragments = [[0], [1, 2]]
s.fragment_charges = [-2, 2]
print(s.fragments)
print(s.fragment_charges)
s.obmol.spin = 20
print(s.multiplicity)
