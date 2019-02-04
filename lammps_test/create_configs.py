import random

class wpacket:

    def __init__(self, mass, charge, width, coord):
        self.mass = mass
        self.charge = charge
        self.width = width
        self.x, self.y, self.z = coord

    def eopt_str(self, i, type):
        text = "[{0}{1}]\n".format(type, i) + \
               "mass={0}\n".format(self.mass) + \
               "charge={0}\n".format(self.charge) + \
               "gamma={0}\n".format(self.width) + \
               "x={0}\ny={1}\nz={2}\n\n".format(self.x, self.y, self.z)
        return text

    def lammps_str(self, i, type, etag):
        text = " ".join(str(v) for v in [
           i, type, self.charge, type-1, self.width, etag, 1.0, 0.0, self.x, self.y, self.z
        ])

        return text+"\n"


class eopt_io:

    def __init__(self, path):
        self.path = path

    def store_system(self, ions, electrons, bound):

        with open(self.path, "w") as out:
            header = "[METHOD]\nApprox = HARTREE\nDFTApprox = Void\nMeshSize = 100\nSpaceSize = {0}\n[VM]\nUpper=1\nLower=1\n".format(bound[0]*2)
#            header = "[METHOD]\nApprox = UHF\nDFTApprox = Void\nMeshSize = 100\nSpaceSize = {0}\n[VM]\nUpper=1\nLower=1\n".format(bound[0]*2)
            out.write(header)
            out.write("[OPTIMIZATION]\nMinimize=Total\n\n")

            pindex = 0
            while pindex < len(ions):
                out.write(ions[pindex].eopt_str(pindex, 'I'))
                out.write(electrons[pindex].eopt_str(pindex, 'E'))

                pindex += 1


class lammps_io:

    def __init__(self, path):
        self.path = path

    def store_system(self, ions, electrons, bound):
        with open(self.path, "w") as out:
            out.writelines("Create by autotest system\n")
            out.writelines("{0} atoms\n".format(len(ions) + len(electrons)))
            out.writelines("2 atom types\n")

            out.writelines("{0} {1} xlo xhi\n".format(-bound[0], bound[0]))
            out.writelines("{0} {1} ylo yhi\n".format(-bound[1], bound[1]))
            out.writelines("{0} {1} zlo zhi\n".format(-bound[2], bound[2]))

            out.writelines("Masses\n\n")
            out.writelines("1 1.000794\n2 0.000544616997098749\n\n")
            out.writelines("Atoms\n\n")

            atom_index = 1
            for ion in ions:
                out.write(ion.lammps_str(atom_index, 1, 0))
                atom_index += 1

            etag = 1
            for electron in electrons:
                out.write(electron.lammps_str(atom_index, 2, etag))
                atom_index += 1
                etag += 1

def create_single_ion(bound):
    mass = 1836
    charge = 1
    w = 1
    c = [random.uniform(-bound[i], bound[i]) for i in [0, 1, 2]]

    return wpacket(mass, charge, w, c)


def create_single_electron(bound):
    mass = 1
    charge = -1
    w = random.uniform(0, bound[3])
    c = [random.uniform(-bound[i], bound[i]) for i in [0, 1, 2]]

    return wpacket(mass, charge, w, c)


def create_system(size, bound, id):
    ions = []
    electrons = []

    while size > 0:
        ions.append(create_single_ion(bound))
        electrons.append(create_single_electron(bound))
        size -= 1

    e_io = eopt_io("inconfig/{0}.ini".format(id))
    e_io.store_system(ions, electrons, bound)

    l_io = lammps_io("inconfig/{0}.data".format(id))
    l_io.store_system(ions, electrons, bound)


import sys

create_system(int(sys.argv[1]), [5, 5, 5, 5], sys.argv[1])
