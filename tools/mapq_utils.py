import os
import numpy
import qscores


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
                elif r.type.upper() in chargedIons:
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






