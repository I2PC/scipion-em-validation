import os
import numpy
import qscores
import mmcif
import chimera
from chimera.resCode import protein3to1, nucleic3to1

def ReadMol(fpath, log=False):

    from random import random

    cif, loops = mmcif.ReadCif(fpath, log)

    # descriptions by chain id:
    descrByEntityId = mmcif.GetEntityDescr(cif, loops)

    try:
        atoms = loops['_atom_site']['data']
        print(" - {} atom records".format(len(atoms)))
    except:
        print(" - no atoms in cif?")
        return None

    labels = loops['_atom_site']['labels']
    if 0:
        print("Labels:")
        for l in labels:
            print(" : {}".format(l))

    import time
    start = time.time()

    rmap = {}

    nmol = chimera.Molecule()
    from os import path
    # nmol.name = path.splitext(path.split(fpath)[1])[0]
    nmol.name = path.split(fpath)[1]
    nmol.openedAs = [fpath, []]
    nmol.cif = cif
    nmol.cifLoops = loops

    nmol.chainColors = {}
    nmol.chainDescr = {}

    numQ = 0
    first = True
    for at in atoms:
        mp = at['asMap']

        if log and first:
            for li, label in enumerate(labels):
                print("   {} : {} : {}".format(li+1, label, mp[label]))

        first = False

        atType = mp['type_symbol']
        atName = mp['label_atom_id']
        rtype = mp['label_comp_id']
        chainId = mp['auth_asym_id']
        chainEId = mp['label_entity_id']
        px = mp['Cartn_x']
        py = mp['Cartn_y']
        pz = mp['Cartn_z']
        occ = mp['occupancy']
        bfactor = mp['B_iso_or_equiv']
        altLoc = mp['label_alt_id']
        if altLoc == ".":
            altLoc = ''

        if chainEId in descrByEntityId:
            nmol.chainDescr[chainId] = descrByEntityId[chainEId]

        resId = mmcif.ResId(mp)
        if resId is None:
            continue

        ris = "{}{}".format(chainId, resId)
        res = None
        if ris not in rmap:
            res = nmol.newResidue(rtype, chimera.MolResId(chainId, resId))
            rmap[ris] = res
        else:
            res = rmap[ris]

        clr = None
        if chainId not in nmol.chainColors:
            clr = chimera.MaterialColor(random(), random(), random(), 1.0)
            nmol.chainColors[chainId] = clr
            if 0 and log:
                print(" - chain {}".format(chainId))
        else:
            clr = nmol.chainColors[chainId]

        nat = nmol.newAtom(atName, chimera.Element(atType))

        drawRib = rtype in protein3to1 or rtype in nucleic3to1

        # aMap[at] = nat
        res.addAtom(nat)
        nat.setCoord(chimera.Point(float(px), float(py), float(pz)))
        nat.altLoc = altLoc
        nat.occupancy = float(occ)
        nat.bfactor = float(bfactor)

        if 'Q-score' in mp:
            try:
                Q = float(mp['Q-score'])
                nat.Q = Q
                numQ += 1
            except:
                # print(f" - q score is {mp['Q-score']}")
                pass

    end = time.time()
    print(" - created {} atoms, {:.1f}s, {} q-scores".format(len(nmol.atoms), end-start, numQ))

    return nmol


def SaveQStats(mol, chainId, sigma, RES=3.0):
    '''
    Adapted from SaveQStats() from the MapQ plugin version 1.8.2. 
    https://github.com/gregdp/mapq/
    '''

    if chainId is None:
        chainId = "All"

    cres = {}
    for r in mol.residues:
        if r.id.chainId == chainId or chainId == "All":
            cres.setdefault(r.id.chainId, []).append([r.id.position, r])

    molPath = os.path.splitext(mol.openedAs[0])[0]
    nname = molPath + "__Q__" + chainId + ".txt"

    print("\nSaving per-chain & per-residue Q-scores:")
    print(" -> res = %s" % str(RES))
    print(" -> chain: %s" % chainId)
    print("\nMapQ Statistics saved at %s" % nname)

    with open(nname, "w") as fp:
        fp.write("\n")
        fp.write("Resolution entered (RES): %g\n" % RES)
        fp.write("Model: %s\n" % mol.name)
        fp.write("Sigma: %g\n" % sigma)
        fp.write("\n")

        avgQrna, eq_nucleic = qscores.eQ_nucleic(RES, sigma)
        avgQprot, eq_protein = qscores.eQ_protein(RES, sigma)
        avgQIon, eq_ion =  qscores.eQ_ion(RES, sigma)
        avgQWater, eq_water =  qscores.eQ_water(RES, sigma)

        fp.write("Protein: expectedQ = %s\n" % eq_protein)
        fp.write("Nucleic: expectedQ = %s\n" % eq_nucleic)
        fp.write("Ion: expectedQ = %s\n" % eq_ion)
        fp.write("Water: expectedQ = %s\n" % eq_water)
        fp.write("\n")
        fp.write("Chain\tType\t# residues\tAvg. Q\tExpectedQ@%.2f\tEst.Res.\n" % RES)

        chains = list(cres.keys())
        chains.sort()

        for cid in chains:
            ress = cres[cid]
            type_ats = {}
            type_ress = {}
            resAtoms = []

            for ri, r in ress:
                tp = ""
                if r.isProt:
                    tp = "Protein"
                elif r.isNA:
                    tp = "Nucleic"
                elif r.type.upper() in qscores.chargedIons:
                    tp = "Ion"
                elif r.type.upper() == "HOH":
                    tp = "Water"
                else:
                    tp = r.type

                if tp in type_ats:
                    type_ats[tp].extend(r.atoms)
                else:
                    type_ats[tp] = r.atoms[:]

                if tp in type_ress:
                    type_ress[tp].append(r)
                else:
                    type_ress[tp] = [r]

            for rtype, atoms in type_ats.items():
                qs = [at.Q for at in atoms if (at.element.name != "H" and hasattr(at,'Q'))]
                if len(qs) == 0:
                    continue
                avgQ = numpy.average(qs)
                numR = len(type_ress[rtype])

                formula, estRes = None, None
                if "Protein" in rtype:
                    formula = "=" + qscores.eQ_protein(RES,sigma)[1].replace("RES",'%.2f') % RES
                    estRes = (avgQ - 1.1192) / -0.1775
                elif "Nucleic" in rtype:
                    formula ="=" + qscores.eQ_nucleic(RES,sigma)[1].replace("RES",'%.2f') % RES
                    estRes = (avgQ - 0.9973) / -0.1377
                elif "Ion" in rtype:
                    formula = "=" + qscores.eQ_ion(RES,sigma)[1].replace("RES",'%.2f') % RES
                    estRes = (avgQ - 1.0795) / -0.1103
                elif "Water" in rtype:
                    formula ="=" + qscores.eQ_water(RES,sigma)[1].replace("RES",'%.2f') % RES
                    estRes = (avgQ - 1.0001) / -0.0895
                else:
                    formula = "?"
                    estRes = 0.0

                fp.write("%s\t%s\t%d\t%.2f\t%s\t%.2f\n" % (cid, rtype, numR, avgQ, formula, estRes))

        print("\n")






