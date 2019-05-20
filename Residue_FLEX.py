#!/Users/stephanedecesco/anaconda/bin/python



# ------------------------ IMPORT LIBRARY --------------------------------#
import sys, rlcompleter, readline, glob, os, shutil, math, copy
import numpy as np
from math import cos, sin

# set the print option to print 3 decimals even if last one is a zero
np.set_printoptions(formatter={'float': '{: 0.3f}'.format})


# ----------------------- TAB: AUTOCOMPLETE ------------------------------#
def complete(text, state):
    return (glob.glob(text + '*') + [None])[state]


readline.set_completer_delims(' \t\n;')
readline.set_completer(complete)
readline.parse_and_bind('tab: complete')


# ---------------------------- ATOM CLASS --------------------------------#
class Atom(object):
    def __init__(self, atom_ID, atom_type, x, y, z, rotate=False, residue=False):
        self.ID = atom_ID
        self.name = atom_type
        self.residue = residue
        self.coordinate = np.array([x, y, z])
        self.coordinate = np.around(self.coordinate, 3)
        self.rotate = rotate

    # ---------------- PRINT ATOM INFO (dev_purpose)--------------------------#
    def show(self):
        return "Atom ID: " + str(self.ID) + " Atom name ID: " + str(self.name) + " Residue: " + str(
            self.residue) + " Coordinate: " + str(self.coordinate) + " Will rotate: " + str(self.rotate) + ""


# ------ Functions to define objects and variables and modify them -------#
def dihedral_atoms_list(in_chi):
    # --------------------- DIHEDRAL ATOMS DICTIONARY -----------------------#
    # List of atoms involved in dihedral angle calculation for every amino acid and for all
    # the different CHI angle possible
    chi1_dictionary = {'ARG': ['N', 'CA', 'CB', 'CG'],
                       'ASN': ['N', 'CA', 'CB', 'CG'],
                       'ASP': ['N', 'CA', 'CB', 'CG'],
                       'GLN': ['N', 'CA', 'CB', 'CG'],
                       'GLU': ['N', 'CA', 'CB', 'CG'],
                       'HIS': ['N', 'CA', 'CB', 'CG'],
                       'LEU': ['N', 'CA', 'CB', 'CG'],
                       'LYS': ['N', 'CA', 'CB', 'CG'],
                       'MET': ['N', 'CA', 'CB', 'CG'],
                       'PHE': ['N', 'CA', 'CB', 'CG'],
                       'PRO': ['N', 'CA', 'CB', 'CG'],
                       'TRP': ['N', 'CA', 'CB', 'CG'],
                       'TYR': ['N', 'CA', 'CB', 'CG'],
                       'CYS': ['N', 'CA', 'CB', 'SG'],
                       'ILE': ['N', 'CA', 'CB', 'CG1'],
                       'SER': ['N', 'CA', 'CB', 'OG'],
                       'THR': ['N', 'CA', 'CB', 'OG1'],
                       'VAL': ['N', 'CA', 'CB', 'CG1']}
    chi2_dictionary = {'ARG': ['CA', 'CB', 'CG', 'CD'],
                       'GLN': ['CA', 'CB', 'CG', 'CD'],
                       'GLU': ['CA', 'CB', 'CG', 'CD'],
                       'LYS': ['CA', 'CB', 'CG', 'CD'],
                       'PRO': ['CA', 'CB', 'CG', 'CD'],
                       'ASN': ['CA', 'CB', 'CG', 'OD1'],
                       'ASP': ['CA', 'CB', 'CG', 'OD1'],
                       'HIS': ['CA', 'CB', 'CG', 'ND1'],
                       'ILE': ['CA', 'CB', 'CG1', 'CD'],
                       'LEU': ['CA', 'CB', 'CG', 'CD1'],
                       'MET': ['CA', 'CB', 'CG', 'SD'],
                       'PHE': ['CA', 'CB', 'CG', 'CD1'],
                       'TYR': ['CA', 'CB', 'CG', 'CD1'],
                       'TRP': ['CA', 'CB', 'CG', 'CD1']}
    chi3_dictionary = {'ARG': ['CB', 'CG', 'CD', 'NE'],
                       'GLN': ['CB', 'CG', 'CD', 'OE1'],
                       'GLU': ['CB', 'CG', 'CD', 'OE1'],
                       'LYS': ['CB', 'CG', 'CD', 'CE'],
                       'MET': ['CB', 'SG', 'CD', 'CE']}
    chi4_dictionary = {'ARG': ['CG', 'CD', 'NE', 'CZ'],
                       'LYS': ['CG', 'CD', 'CE', 'NZ']}
    chi5_dictionary = {'ARG': ['CD', 'NE', 'CZ', 'NH1']}
    # -------------- SETTING ATOMS PART OF DIHEDRAL ANGLE -------------------#
    chi_check = 0
    while chi_check == 0:
        if in_chi == "CHI1":
            dihedral_list = chi1_dictionary[residue_name]
            chi_check = 1
        elif in_chi == "CHI2":
            dihedral_list = chi2_dictionary[residue_name]
            chi_check = 1
        elif in_chi == "CHI3":
            dihedral_list = chi3_dictionary[residue_name]
            chi_check = 1
        elif in_chi == "CHI4":
            dihedral_list = chi4_dictionary[residue_name]
            chi_check = 1
        elif in_chi == "CHI5":
            dihedral_list = chi5_dictionary[residue_name]
            chi_check = 1
        else:
            in_chi = raw_input(
                'Chi angle not well inputed. Correct format -->(CHI1 or CHI2 or CHI3 or CHI4 or CHI5):').upper()
    dihedral_list.append(in_chi)
    return dihedral_list


def get_chi(residue_name):
    # --------------------- DIHEDRAL ATOMS DICTIONARY -----------------------#
    # List of atoms involved in dihedral angle calculation for every amino acid and for all
    # the different CHI angle possible
    chi1_dictionary = {'ARG': ['N', 'CA', 'CB', 'CG'],
                       'ASN': ['N', 'CA', 'CB', 'CG'],
                       'ASP': ['N', 'CA', 'CB', 'CG'],
                       'GLN': ['N', 'CA', 'CB', 'CG'],
                       'GLU': ['N', 'CA', 'CB', 'CG'],
                       'HIS': ['N', 'CA', 'CB', 'CG'],
                       'LEU': ['N', 'CA', 'CB', 'CG'],
                       'LYS': ['N', 'CA', 'CB', 'CG'],
                       'MET': ['N', 'CA', 'CB', 'CG'],
                       'PHE': ['N', 'CA', 'CB', 'CG'],
                       'PRO': ['N', 'CA', 'CB', 'CG'],
                       'TRP': ['N', 'CA', 'CB', 'CG'],
                       'TYR': ['N', 'CA', 'CB', 'CG'],
                       'CYS': ['N', 'CA', 'CB', 'SG'],
                       'ILE': ['N', 'CA', 'CB', 'CG1'],
                       'SER': ['N', 'CA', 'CB', 'OG'],
                       'THR': ['N', 'CA', 'CB', 'OG1'],
                       'VAL': ['N', 'CA', 'CB', 'CG1']}
    chi2_dictionary = {'ARG': ['CA', 'CB', 'CG', 'CD'],
                       'GLN': ['CA', 'CB', 'CG', 'CD'],
                       'GLU': ['CA', 'CB', 'CG', 'CD'],
                       'LYS': ['CA', 'CB', 'CG', 'CD'],
                       'PRO': ['CA', 'CB', 'CG', 'CD'],
                       'ASN': ['CA', 'CB', 'CG', 'OD1'],
                       'ASP': ['CA', 'CB', 'CG', 'OD1'],
                       'HIS': ['CA', 'CB', 'CG', 'ND1'],
                       'ILE': ['CA', 'CB', 'CG1', 'CD'],
                       'LEU': ['CA', 'CB', 'CG', 'CD1'],
                       'MET': ['CA', 'CB', 'CG', 'SD'],
                       'PHE': ['CA', 'CB', 'CG', 'CD1'],
                       'TYR': ['CA', 'CB', 'CG', 'CD1'],
                       'TRP': ['CA', 'CB', 'CG', 'CD1']}
    chi3_dictionary = {'ARG': ['CB', 'CG', 'CD', 'NE'],
                       'GLN': ['CB', 'CG', 'CD', 'OE1'],
                       'GLU': ['CB', 'CG', 'CD', 'OE1'],
                       'LYS': ['CB', 'CG', 'CD', 'CE'],
                       'MET': ['CB', 'SG', 'CD', 'CE']}
    chi4_dictionary = {'ARG': ['CG', 'CD', 'NE', 'CZ'],
                       'LYS': ['CG', 'CD', 'CE', 'NZ']}
    chi5_dictionary = {'ARG': ['CD', 'NE', 'CZ', 'NH1']}
    # -------------- SETTING ATOMS PART OF DIHEDRAL ANGLE -------------------#
    chi_list = []
    if residue_name in chi1_dictionary:
        chi_list.append(chi1_dictionary[residue_name])
        chi_list[0].append("CHI1")
    if residue_name in chi2_dictionary:
        chi_list.append(chi2_dictionary[residue_name])
        chi_list[1].append("CHI2")
    if residue_name in chi3_dictionary:
        chi_list.append(chi3_dictionary[residue_name])
        chi_list[2].append("CHI3")
    if residue_name in chi4_dictionary:
        chi_list.append(chi4_dictionary[residue_name])
        chi_list[3].append("CHI4")
    if residue_name in chi5_dictionary:
        chi_list.append(chi5_dictionary[residue_name])
        chi_list[4].append("CHI5")
    return chi_list


def dihedral_atoms_objects(dihedral_list, residue_atoms):
    dihedral_atoms = []
    for i in range(len(residue_atoms)):
        if dihedral_list[0] == str(residue_atoms[i].name):
            dihedral_atoms.append(residue_atoms[i])
            for j in range(len(residue_atoms)):
                if dihedral_list[1] == str(residue_atoms[j].name):
                    dihedral_atoms.append(residue_atoms[j])
                    for k in range(len(residue_atoms)):
                        if dihedral_list[2] == str(residue_atoms[k].name):
                            dihedral_atoms.append(residue_atoms[k])
                            for l in range(len(residue_atoms)):
                                if dihedral_list[3] == str(residue_atoms[l].name):
                                    dihedral_atoms.append(residue_atoms[l])
    return dihedral_atoms


def set_rotation(residue_atoms, dihedral_list, reset=0):
    backbone_list = ['C', 'N', 'CA', 'O']
    check = 0
    after = False
    while check == 0:
        if reset == 0:
            for i in range(len(residue_atoms)):
                if residue_atoms[i].name == dihedral_list[0] or after:
                    if str(residue_atoms[i].name) not in str(dihedral_list[0:3]) and str(
                            residue_atoms[i].name) not in str(backbone_list):
                        residue_atoms[i].rotate = True
                    after = True
            check = 1
        elif reset == 1:
            for i in range(len(residue_atoms)):
                residue_atoms[i].rotate = False
            check = 1
        else:
            reset = int(raw_input('Please enter a valid value of reset (Yes=1 // No = 0): '))


def create_atoms_list(file_name):
    with open(file_name, 'r') as pdb:
        # store all the lines
        lines = pdb.readlines()
        atoms_list = []
        for line in lines:
            atoms_list.append(line.split(" "))
    return atoms_list


def create_residue(atoms_list, residue_name, residue_number, in_chain):
    residue_To_Move = []
    residue_atoms = []
    for i in range(len(atoms_list)):
        if residue_name in atoms_list[i] and residue_number in atoms_list[i] and "ATOM" in atoms_list[i]:
            residue_To_Move.append(filter(None, atoms_list[i]))
    # ------------------ CHECK IF RESIDUE EXISTS ----------------------------#
    while not residue_To_Move:
        in_residue = raw_input(
            'Residue not found, make sure to use proper format -->(3 Letter code + residue number attached. e.g. TYR568): ')
        residue_name = in_residue[0:3]
        residue_number = in_residue[3:]
        for i in range(len(atoms_list)):
            if residue_name in atoms_list[i] and residue_number in atoms_list[i] and "ATOM" in atoms_list[i]:
                residue_To_Move.append(filter(None, atoms_list[i]))
    # ------------------------ SETTING "residue_atoms"-----------------------#
    for i in range(len(residue_To_Move)):
        if in_chain == residue_To_Move[i][4]:
            residue_atoms.append(
                Atom(int(residue_To_Move[i][1]), str(residue_To_Move[i][2]), float(residue_To_Move[i][6]),
                     float(residue_To_Move[i][7]), float(residue_To_Move[i][8])))
            residue_atoms[i].residue = residue_To_Move[i][3] + "" + residue_To_Move[i][5]
    return (residue_atoms, residue_name, residue_number)


def create_protein(atoms_list):
    protein_atoms = []
    all_atoms = []
    for i in range(len(atoms_list)):
        if "ATOM" in atoms_list[i][0]:
            all_atoms.append(filter(None, atoms_list[i]))
    for i in range(len(all_atoms)):
        protein_atoms.append(
            Atom(int(all_atoms[i][1]), str(all_atoms[i][2]), float(all_atoms[i][6]), float(all_atoms[i][7]),
                 float(all_atoms[i][8])))
        protein_atoms[i].residue = all_atoms[i][3] + "" + all_atoms[i][5]
    return protein_atoms


# -------- METHODS USED TO ROTATE AND MEASURE DIHEDRAL ANGLE -------------#
def get_dihedral(atom):
    """Calculate dihedral angle.

		Calculate dihedral angle between the vectors list[0]->list[1]
		and list[2]->list[3], where list contains the atomic indexes
		in question.
		"""

    # vector 0->1, 1->2, 2->3 and their normalized cross products:
    a = atom[1].coordinate - atom[0].coordinate
    b = atom[2].coordinate - atom[1].coordinate
    c = atom[3].coordinate - atom[2].coordinate
    bxa = np.cross(b, a)
    bxa /= np.linalg.norm(bxa)
    cxb = np.cross(c, b)
    cxb /= np.linalg.norm(cxb)
    angle = np.vdot(bxa, cxb)
    # check for numerical trouble due to finite precision:
    if angle < -1:
        angle = -1
    if angle > 1:
        angle = 1
    angle = np.arccos(angle)
    if np.vdot(bxa, c) > 0:
        angle = 2 * np.pi - angle
    """
		Return an angle in radians
		Use math.degrees(angle_in_radian) to convert
		"""
    return math.degrees(angle)


def set_dihedral(residue_at, dihedral_at, angle):
    """Set the dihedral angle between vectors list[0]->list[1] and
	list[2]->list[3] by changing the atom indexed by list[3]
	if mask is not None, all the atoms described in mask
	(read: the entire subgroup) are moved. Alternatively to the mask,
	the indices of the atoms to be rotated can be supplied.

	example: the following defines a very crude
	ethane-like molecule and twists one half of it by 30 degrees.

	>>> atoms = Atoms('HHCCHH', [[-1, 1, 0], [-1, -1, 0], [0, 0, 0],
								 [1, 0, 0], [2, 1, 0], [2, -1, 0]])
	>>> atoms.set_dihedral([1,2,3,4],7*pi/6,mask=[0,0,0,1,1,1])
	"""
    # compute necessary in dihedral change, from current value
    current = math.radians(get_dihedral(dihedral_at))
    diff = math.radians(angle) - current
    axis = dihedral_at[2].coordinate - dihedral_at[1].coordinate
    center = dihedral_at[2].coordinate

    return _masked_rotate(center, axis, diff, residue_at)


def scan_dihedral(residue_at, dihedral_at, angle):
    """Set the dihedral angle between vectors list[0]->list[1] and
	list[2]->list[3] by changing the atom indexed by list[3]
	if mask is not None, all the atoms described in mask
	(read: the entire subgroup) are moved. Alternatively to the mask,
	the indices of the atoms to be rotated can be supplied.

	example: the following defines a very crude
	ethane-like molecule and twists one half of it by 30 degrees.

	>>> atoms = Atoms('HHCCHH', [[-1, 1, 0], [-1, -1, 0], [0, 0, 0],
								 [1, 0, 0], [2, 1, 0], [2, -1, 0]])
	>>> atoms.set_dihedral([1,2,3,4],7*pi/6,mask=[0,0,0,1,1,1])
	"""

    # compute necessary in dihedral change, from current value
    current = math.radians(get_dihedral(dihedral_at))
    diff = math.radians(angle)
    axis = dihedral_at[2].coordinate - dihedral_at[1].coordinate
    center = dihedral_at[2].coordinate
    # 	print current
    # 	print math.degrees(diff)
    # 	print axis
    # 	print center
    return _masked_rotate(center, axis, diff, residue_at)


def scan_loop(residue_atoms, dihedral_list, dihedral_atoms, scan_angle, previous_run, number=0):
    index = len(dihedral_list) - number
    if number > 1:
        total = 0
        while total < 360:
            # ----------------------- ROTATION OF THE RESIDUE -----------------------#
            set_rotation(residue_atoms, dihedral_list[index], reset=1)
            set_rotation(residue_atoms, dihedral_list[index])
            dihedral_atoms[index] = dihedral_atoms_objects(dihedral_list[index], residue_atoms)
            scan_dihedral(residue_atoms, dihedral_atoms[index], scan_angle[index])
            total += int(scan_angle[index])
            scan_loop(residue_atoms, dihedral_list, dihedral_atoms, scan_angle, previous_run, number=number - 1)
    else:
        total = 0
        while total < 360:
            # ----------------------- ROTATION OF THE RESIDUE -----------------------#
            set_rotation(residue_atoms, dihedral_list[index], reset=1)
            set_rotation(residue_atoms, dihedral_list[index])
            dihedral_atoms[index] = dihedral_atoms_objects(dihedral_list[index], residue_atoms)
            scan_dihedral(residue_atoms, dihedral_atoms[index], scan_angle[index])
            total += int(scan_angle[index])
            # ----------------------- CHECK FOR CLASHES------------------------------#
            if in_contact(residue_atoms, protein_atoms):
                for i in range(len(dihedral_list)):
                    print "The dihedral angle for " + str(dihedral_list[i]) + " is = " + str(
                        get_dihedral(dihedral_atoms[i]))
                print "/!\ ------------ STERIC CLASH -----------------/!\\"
                print "--------------- STRUCTURE SKIPPED ----------------\n\n\n"
            else:
                # ----------------------- IF NO CLASH CREATE FILE------------------------#
                for i in range(len(dihedral_list)):
                    print "The dihedral angle for " + str(dihedral_list[i]) + " is = " + str(
                        get_dihedral(dihedral_atoms[i]))
                print " ---------------- NO CLASH DETECTED ------------------------"

                name = []
                for i in range(len(dihedral_list)):
                    name.append(
                        "_" + str(dihedral_list[i][4]) + "_" + str(int(get_dihedral(dihedral_atoms[i]))) + "DEG")
                desc = ''.join(name)
                file_name = output_name + "_" + str(residue_name) + "" + str(residue_number) + "" + str(desc) + ".pdb"
                similar = False
                for i in range(len(previous_run)):
                    if RMSD(previous_run[i], residue_atoms) < 0.1:
                        similar = True
                if similar:
                    print "FILE NOT CREATED : CONFORMATION ALREADY EXISTING"
                    print "Cause: Symmetry of " + str(residue_name)
                else:
                    if write_file(file_name, atoms_list, residue_atoms):
                        print "FILE CREATED: (" + file_name + ")\n\n\n"
                        previous_run.append(copy.deepcopy(residue_atoms))
                    else:
                        print "An ERROR occured during file creation process"


def _masked_rotate(center, axis, diff, residue_at):
    # do rotation of subgroup by copying it to temporary atoms object
    # and then rotating that
    #
    # recursive object definition might not be the most elegant thing,
    # more generally useful might be a rotation function with a mask?
    group = []
    for i in range(len(residue_at)):
        if residue_at[i].rotate == True:
            group.append(residue_at[i])
    group = translate(group, -center)
    group = rotate(group, axis, diff)
    group = translate(group, center)
    # set positions in original atoms object
    j = 0
    for i in range(len(residue_atoms)):
        if residue_at[i].rotate == True:
            residue_at[i].coordinate = np.around(group[j].coordinate, 3)
            j += 1
    return residue_at


def rotate(group, v, a=None, center=(0, 0, 0)):
    """Rotate atoms based on a vector and an angle, or two vectors.

	Parameters:

	v:
		Vector to rotate the atoms around. Vectors can be given as
		strings: 'x', '-x', 'y', ... .

	a = None:
		Angle that the atoms is rotated around the vecor 'v'. If an angle
		is not specified, the length of 'v' is used as the angle
		(default). The angle can also be a vector and then 'v' is rotated
		into 'a'.

	center = (0, 0, 0):
		The center is kept fixed under the rotation. Use 'COM' to fix
		the center of mass, 'COP' to fix the center of positions or
		'COU' to fix the center of cell.

	Examples:

	Rotate 90 degrees around the z-axis, so that the x-axis is
	rotated into the y-axis:

	>>> a = pi / 2
	>>> atoms.rotate('z', a)
	>>> atoms.rotate((0, 0, 1), a)
	>>> atoms.rotate('-z', -a)
	>>> atoms.rotate((0, 0, a))
	>>> atoms.rotate('x', 'y')
	"""

    norm = np.linalg.norm
    if a is None:
        a = norm(v)
    if isinstance(a, (float, int)):
        v /= norm(v)
        c = cos(a)
        s = sin(a)
    for i in range(len(group)):
        p = group[i].coordinate - center
        group[i].coordinate[:] = (c * p -
                                  np.cross(p, s * v) +
                                  np.outer(np.dot(p, v), (1.0 - c) * v) +
                                  center)
    return group


def translate(group, displacement):
    """Translate atomic positions.

	The displacement argument can be a float an xyz vector or an
	nx3 array (where n is the number of atoms)."""
    for i in range(len(group)):
        group[i].coordinate += np.array(displacement)
    return group


def RMSD(residue_new, residue_old):
    number_atoms = len(residue_new)
    sumd_dist_sq = 0
    residue_name = residue_new[1].residue[0:3]
    residue_symmetry = {'ARG': 1,
                        'ASN': 0,
                        'ASP': 1,
                        'GLN': 0,
                        'GLU': 1,
                        'HIS': 0,
                        'LEU': 1,
                        'LYS': 0,
                        'MET': 0,
                        'PHE': 1,
                        'PRO': 0,
                        'TRP': 0,
                        'TYR': 1,
                        'CYS': 0,
                        'ILE': 0,
                        'SER': 0,
                        'THR': 0,
                        'VAL': 1}

    for i in range(len(residue_new)):
        if residue_symmetry[residue_name] == 0:
            sumd_dist_sq += distance(residue_new[i].coordinate, residue_old[i].coordinate) * distance(
                residue_new[i].coordinate, residue_old[i].coordinate)
        else:
            if "1" in residue_new[i].name:
                sumdist = distance(residue_new[i].coordinate, residue_old[i].coordinate) * distance(
                    residue_new[i].coordinate, residue_old[i].coordinate)
                sumdist_2 = distance(residue_new[i].coordinate, residue_old[i + 1].coordinate) * distance(
                    residue_new[i].coordinate, residue_old[i + 1].coordinate)
                if sumdist_2 < sumdist:
                    sumd_dist_sq += sumdist_2
                else:
                    sumd_dist_sq += sumdist
            elif "2" in residue_new[i].name:
                sumdist = distance(residue_new[i].coordinate, residue_old[i].coordinate) * distance(
                    residue_new[i].coordinate, residue_old[i].coordinate)
                sumdist_2 = distance(residue_new[i].coordinate, residue_old[i - 1].coordinate) * distance(
                    residue_new[i].coordinate, residue_old[i - 1].coordinate)
                if sumdist_2 < sumdist:
                    sumd_dist_sq += sumdist_2
                else:
                    sumd_dist_sq += sumdist
            else:
                sumd_dist_sq += distance(residue_new[i].coordinate, residue_old[i].coordinate) * distance(
                    residue_new[i].coordinate, residue_old[i].coordinate)

    rmsd = sumd_dist_sq / number_atoms
    if rmsd > 0:
        rmsd = math.sqrt(rmsd)
    else:
        rmsd = 0
    return rmsd


# --------------------- METHODS USED FOR CLASH ---------------------------#
def mag(vector):
    """
  Returns the magnitude of a vector
  """
    return np.sqrt(np.dot(vector, vector))


def distance(p1, p2):
    """
  Returns distance between the two points
  """
    return mag(p1 - p2)


def in_contact(residue, protein):
    backbone_list = ['C', 'N', 'CA', 'O']
    for i in range(len(residue)):
        for j in range(len(protein)):
            if residue[i].residue != protein[j].residue:
                if distance(residue[i].coordinate, protein[j].coordinate) < 2 and str(residue[i].name) not in str(
                        backbone_list):
                    return True
    return False


# ------------------- METHOD USED FOR WRITING FILE -----------------------#
def write_file(file_name, atoms_list, residue_atoms):
    with open(file_name, 'w') as output:
        non_empty_param = []
        for i in range(len(atoms_list)):
            non_empty_param.append(list(np.nonzero(atoms_list[i])[0]))
        for i in range(len(residue_atoms)):
            for j in range(len(atoms_list)):
                if "ATOM" == str(atoms_list[j][0]) and int(residue_atoms[i].ID) == int(
                        atoms_list[j][non_empty_param[j][1]]):
                    atoms_list[j][non_empty_param[j][6]] = "{0:.3f}".format(residue_atoms[i].coordinate[0])
                    atoms_list[j][non_empty_param[j][7]] = "{0:.3f}".format(residue_atoms[i].coordinate[1])
                    atoms_list[j][non_empty_param[j][8]] = "{0:.3f}".format(residue_atoms[i].coordinate[2])
        for j in range(len(atoms_list)):
            write = ' '.join(atoms_list[j])
            output.write(write)
    return True


# ---------------------------- INPUTS ------------------------------------#
in_mode = 0
# Asking for input files name and output name
while True:
    in_PDB = raw_input('Enter the PDB file name you want to mutate the residue: ')
    if os.path.isfile(in_PDB):
        break

in_chain = raw_input('Enter the chain of the residue (e.g. A,B or C): ').upper()
in_residue = raw_input('Enter the residue you want to screen conformations (e.g. TYR568): ')

# ------------------ SET RESIDUE NAME AND NUMBER  -----------------------#
residue_name = in_residue[0:3]
residue_number = in_residue[3:]
# ------------------------ OPENING PDB FILE -----------------------------#
atoms_list = create_atoms_list(in_PDB)
# extracting from the pdb file only the atoms associated to the wanted residue to rotate into a list of object atoms
# ------------------------ SETTING "residue_atoms"-----------------------#
residue_atoms, residue_name, residue_number = create_residue(atoms_list, residue_name, residue_number, in_chain)
# -------------- PUTTING ALL PROTEIN ATOMS IN OBJECTS--------------------#
protein_atoms = create_protein(atoms_list)
# ------------- RESET RESIDUE NAME AND NUMBER  --------------------------#
chi_list = get_chi(residue_name)

# -------------------- ASKING THE MODE OF ROTATION ----------------------#
# ----------------------- MODE 1 : SET DIHEDRAL  ------------------------#
# ----------------------- MODE 2 : SCAN DIHEDRAL ------------------------#
# ----------------------- MODE 3 : SCAN MULTIPLE CHI --------------------#
# loop to make sure the mode is well set
while in_mode != 1 and in_mode != 2 and in_mode != 3:
    print "------------------------------------- ROTATION MODES --------------------------------------------"
    print "Mode 1: You pick one CHI angle and SET its value to a choosen number"
    print "Mode 2: You pick one CHI angle and SCAN its value by a choosen increment"
    print "Mode 3: You pick multiple CHI angles and SCAN all the possible combination by a choosen increment"
    print "--------------------------------------------------------------------------------------------------"
    print "List of available chi angle"
    print "---------------------------"
    for i in range(len(chi_list)):
        print chi_list[i]
    in_mode = int(raw_input('Choose your mode :(SET = 1 // SCAN = 2 // MULTI_SCAN = 3)'))
    if in_mode == 1:
        # ------------------ ASK FOR THE DIHEDRAL TO CHANGE ---------------------#
        in_chi = raw_input('Input the Chi angle (CHI1 or CHI2 or CHI3 or CHI4 or CHI5):').upper()
        # -------------- SETTING ATOMS PART OF DIHEDRAL ANGLE -------------------#
        dihedral_list = dihedral_atoms_list(in_chi)
        in_chi = dihedral_list[4]
        # ----------------------- SETTING "dihedral_atoms" ----------------------#
        dihedral_atoms = dihedral_atoms_objects(dihedral_list, residue_atoms)
        # ----------------------- SETTING ATOMS THAT ROTATE----------------------#
        set_rotation(residue_atoms, dihedral_list)
        # ----------------------- PRINTING INITIAL DIHEDRAL ANGLE ---------------#
        print "The initial dihedral angle for " + str(dihedral_list) + " is = " + str(get_dihedral(dihedral_atoms))
        # ----------------------- ASK FOR THE ANGLE TO SET ----------------------#
        set_angle = float(raw_input('You choosed to SET the angle to the value (in degrees) = '))
    if in_mode == 2:
        # ------------------ ASK FOR THE DIHEDRAL TO CHANGE ---------------------#
        in_chi = raw_input('Input the Chi angle (CHI1 or CHI2 or CHI3 or CHI4 or CHI5):').upper()
        # -------------- SETTING ATOMS PART OF DIHEDRAL ANGLE -------------------#
        dihedral_list = dihedral_atoms_list(in_chi)
        in_chi = dihedral_list[4]
        # ----------------------- SETTING "dihedral_atoms" ----------------------#
        dihedral_atoms = dihedral_atoms_objects(dihedral_list, residue_atoms)
        # ----------------------- SETTING ATOMS THAT ROTATE----------------------#
        set_rotation(residue_atoms, dihedral_list)
        # ----------------------- PRINTING INITIAL DIHEDRAL ANGLE ---------------#
        print "The initial dihedral angle for " + str(dihedral_list) + " is = " + str(get_dihedral(dihedral_atoms))
        # ----------------------- ASK FOR THE ANGLE TO SET ----------------------#
        scan_angle = float(raw_input('You choosed to SCAN the angle with an increment of (in degrees) = '))
    if in_mode == 3:
        dihedral_list = []
        in_chi = []
        dihedral_atoms = []
        for i in range(len(chi_list)):
            in_chi.append(chi_list[i][4])
            dihedral_list.append(dihedral_atoms_list(in_chi[i]))
            dihedral_atoms.append(dihedral_atoms_objects(dihedral_list[i], residue_atoms))

        # ----------------------- PRINTING INITIALS DIHEDRAL ANGLE ---------------#
        for i in range(len(dihedral_list)):
            print "The initial dihedral angle for " + str(dihedral_list[i]) + " is = " + str(
                int(get_dihedral(dihedral_atoms[i])))
        # ----------------------- ASK FOR THE ANGLE TO SET ----------------------#
        same_angle = raw_input('Do you want the same increment for the two angles (y/n): ')
        scan_angle = []
        while same_angle != "y" and same_angle != "n":
            same_angle = raw_input(
                'Do you want the same increment for the two angles ? (Wrong input:please enter "y" or "n") (y/n): ')
        if same_angle == "y":
            scan_angle_in = float(raw_input('Choose an increment (in degrees) = '))
            for i in range(len(dihedral_list)):
                scan_angle.append(scan_angle_in)
        elif same_angle == "n":
            for i in range(len(dihedral_list)):
                scan_angle.append(
                    float(raw_input('Choose an increment for ' + str(dihedral_list[i]) + ' (in degrees) = ')))

# ----------------------- SETTING OUTPUT NAME ---------------------------#
output_name = raw_input('Enter export file name :')

# atoms_list_2 = create_atoms_list("C.pdb")
# residue_atoms_2 = create_residue(atoms_list_2,residue_name,residue_number)
# print RMSD(residue_atoms,residue_atoms_2)

if in_mode == 1:
    # -----------------------------------------------------------------------#
    # ------------------- MODE 1 : SET DIHEDRAL -----------------------------#
    # -----------------------------------------------------------------------#
    print "-----------------------------------------------------------------------"
    print "------------------- ROTATION WILL NOW START ---------------------------"
    print "-----------------------------------------------------------------------"
    print "\n"
    # ----------------------- ROTATION OF THE RESIDUE -----------------------#
    residue_atoms_old = copy.deepcopy(residue_atoms)
    set_dihedral(residue_atoms, dihedral_atoms, set_angle)

    # ----------------------- CHECK FOR CLASHES------------------------------#
    if in_contact(residue_atoms, protein_atoms):
        print str(in_chi) + " is = " + str(get_dihedral(dihedral_atoms))
        print "/!\ ------------ STERIC CLASH -----------------/!\\"
        print "--------------- STRUCTURE SKIPPED ----------------\n\n\n"
    else:
        # ----------------------- IF NO CLASH CREATE FILE------------------------#
        file_name = output_name + ".pdb"
        if write_file(file_name, atoms_list, residue_atoms):
            print "The new dihedral angle for " + str(in_chi) + " is = " + str(
                get_dihedral(dihedral_atoms)) + "\nFile created (" + str(file_name) + ")"
            print "RMSD with original structure = " + str(RMSD(residue_atoms, residue_atoms_old))
        else:
            print "An ERROR occured during file creation process"
# -----------------------------------------------------------------------#
# -------------------------- END MODE 1 ---------------------------------#
# -----------------------------------------------------------------------#
if in_mode == 2:
    # -----------------------------------------------------------------------#
    # --------------------- MODE 2 : SCAN DIHEDRAL --------------------------#
    # -----------------------------------------------------------------------#
    print "-----------------------------------------------------------------------"
    print "------------------- ROTATION WILL NOW START ---------------------------"
    print "-----------------------------------------------------------------------"
    print "\n"
    total = 0
    previous_run = []
    previous_run.append(copy.deepcopy(residue_atoms))
    while total < 360:
        # ----------------------- ROTATION OF THE RESIDUE -----------------------#
        scan_dihedral(residue_atoms, dihedral_atoms, scan_angle)
        total += int(scan_angle)
        # ----------------------- CHECK FOR CLASHES------------------------------#
        if in_contact(residue_atoms, protein_atoms):
            print str(in_chi) + " is = " + str(get_dihedral(dihedral_atoms))
            print "/!\ ------------ STERIC CLASH -----------------/!\\"
            print "--------------- STRUCTURE SKIPPED ----------------\n\n\n"
        else:
            # ----------------------- IF NO CLASH CREATE FILE------------------------#
            print " ---------------- NO CLASH DETECTED ------------------------"
            print str(in_chi) + " is = " + str(get_dihedral(dihedral_atoms))
            file_name = output_name + "_" + str(residue_name) + "" + str(residue_number) + "_" + str(total) + "DEG.pdb"
            similar = False
            for i in range(len(previous_run)):
                if RMSD(previous_run[i], residue_atoms) < 0.1:
                    similar = True
            if similar:
                print "FILE NOT CREATED : CONFORMATION ALREADY EXISTING"
                print "Cause: Symmetry of " + str(residue_name)
            else:
                if write_file(file_name, atoms_list, residue_atoms):
                    print "FILE CREATED (" + file_name + ")\n\n\n"
                    previous_run.append(copy.deepcopy(residue_atoms))
                else:
                    print "An ERROR occured during file creation process"
# -----------------------------------------------------------------------#
# -------------------------- END MODE 2 ---------------------------------#
# -----------------------------------------------------------------------#
if in_mode == 3:
    # -----------------------------------------------------------------------#
    # ------------------ MODE 3 : SCAN MULTIPLE CHI -------------------------#
    # -----------------------------------------------------------------------#
    print "-----------------------------------------------------------------------"
    print "------------------- ROTATION WILL NOW START ---------------------------"
    print "-----------------------------------------------------------------------"
    print "\n"
    previous_run = []
    previous_run.append(copy.deepcopy(residue_atoms))
    scan_loop(residue_atoms, dihedral_list, dihedral_atoms, scan_angle, previous_run, number=len(dihedral_list))

    # -----------------------------------------------------------------------#
    # -------------------------- END MODE 3 ---------------------------------#
    # -----------------------------------------------------------------------#
