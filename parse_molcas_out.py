import re

def readlog(f):
    in_basis_set = None
    in_mulliken = None
    in_dipole = None
    in_fc_matrix = None
    in_gradients = None
    re_coord = re.compile("^\s*\S+(?:\s+-?\d*\.\d*){3}\s*$")
    re_gradient = re.compile("^\s*\S+(?:\s+-?\d*\.\d*[DdEe][+-]?\d+){3}\s*$")
    re_a_tens_term = re.compile("^(?:\s+-?\d*\.\d*){3}\s*$")
    labels = []
    energy = None
    chg_nuc = []
    chg_ele = []
    dipole = []
    a_tens_term = []
    fermi_contact = []
    gradients = []

    for L in f:
        if L.startswith("::    RASSCF root number  1 Total energy:"):
            energy = float(L.rstrip().split()[-1])

        if in_basis_set:
            if "--" in L:
                in_basis_set = None
            if re_coord.match(L):
                token = L.rstrip().split()
                chg_nuc.append(chg_nuc_tmp)
        if "Associated Effective Charge" in L:
            chg_nuc_tmp = float(L.split()[3])
            in_basis_set = True
            
        if in_mulliken:
            if "Total" in L:
                in_mulliken = None
                chg_ele = [float(q) for q in L.split()[1:]]
        if L.startswith("      Mulliken population analysis for root number:  1"):
            in_mulliken = True

        if in_dipole:
            if "Total" in L:
                in_dipole = None
                token = L.split()
                dipole = [float(q) for q in (token[1], token[3], token[5], token[7])]
        if L.startswith("      Expectation values of various properties for root number:  1"):
            in_dipole = True

        if in_fc_matrix:
            if "A-TENSOR:" in L:
                in_fc_matrix = None
                # 同じ値のはずだが数値誤差でずれるので平均化
                fermi_contact.append((
                    a_tens_term[0][0] +
                    a_tens_term[1][1] +
                    a_tens_term[2][2]) / 3.)
                a_tens_term = []
            if re_a_tens_term.match(L):
                a_tens_term.append([float(term) for term in L.split()])
        if L.startswith("   A (FC)-Matrix for center:"):
            in_fc_matrix = True

        if in_gradients:
            if "Stop" in L:
                in_gradients = None
            if re_gradient.match(L):
                gradients.append([float(g) for g in L.split()[1:]])
        if L.startswith(" *              Molecular gradients               *"):
            in_gradients = True

    with open("punch", r"w") as f:
         f.write(f"{'Energy':<16}{1:>8d}{energy:>24.12f}\n")
         for i, q in enumerate(qn-qe for (qn, qe) in zip(chg_ele, chg_nuc)):
             f.write(f"{'Charge':<16}{i+1:>8d}{q:>24.12f}\n")
         f.write(f"{'DipoleX':<16}{1:>8d}{dipole[0]:>24.12f}\n")
         f.write(f"{'DipoleY':<16}{1:>8d}{dipole[1]:>24.12f}\n")
         f.write(f"{'DipoleZ':<16}{1:>8d}{dipole[2]:>24.12f}\n")
         f.write(f"{'DipoleT':<16}{1:>8d}{dipole[3]:>24.12f}\n")
         for i, fc in enumerate(fermi_contact):
             f.write(f"{'FermiContact':<16}{i+1:>8d}{fc:>24.12f}\n")
         for i, gr in enumerate(gradients):
             f.write(f"{'ForceX':<16}{i+1:>8d}{-gr[0]:>24.12f}\n")
         for i, gr in enumerate(gradients):
             f.write(f"{'ForceY':<16}{i+1:>8d}{-gr[1]:>24.12f}\n")
         for i, gr in enumerate(gradients):
             f.write(f"{'ForceZ':<16}{i+1:>8d}{-gr[2]:>24.12f}\n")

