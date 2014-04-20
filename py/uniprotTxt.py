"""Getting uniprot txt files for proteins
   parsing them
   loading parsed data."""
import os, time, sys
import files
from collections import defaultdict

def parseProteinName(protein):
    with open(mkTxtFile(protein)) as f:
        for line in f:
            # DE   RecName: Full=Fatty acid synthase;
            if line[0:2] == 'DE' and 'RecName' in line:
                name = line.split(':')[1].split(';')[0].split('=')[1]
                return name
            elif line[0:2] == 'DE' and 'SubName' in line:
                name = line.split(':')[1].split(';')[0].split('=')[1]
                return name
            elif line[0:2] == 'GE' and 'Name' in line:
                name = line.split('Name=')[1].split(';')[0]
                return name

def dumpProteinNames(proteins, nameFile):
    """ '../working/neal/proteins.names' """
    with open(nameFile, 'w') as fout:
        for protein in proteins:
            print >> fout, '\t'.join((protein, parseProteinName(protein)))

def loadProteinNames(nameFile):
    protein2name = {}
    with open(nameFile) as f:
        for line in f:
            protein, name = line.strip('\n').split('\t')
            protein2name[protein] = name
    return protein2name

def mkTxtFile(proteinID):
    txtFile = files.dataDir + '/uniprot/txt/' + proteinID + '.txt'
    return txtFile

def mkProteinToProsite(proteins):
    protein2prosite = defaultdict(dict)
    for protein in proteins:
        prosites = parseProsite(protein)
        for prosite in prosites:
            protein2prosite[protein][prosite] = prosites[prosite]
    return protein2prosite

def parseProsite(protein):
    prosites = {}
    with open(mkTxtFile(protein)) as f:
        for line in f:
            # DR   PROSITE; PS51257; PROKAR_LIPOPROTEIN; 1.
            # DR   GO; GO:0046872; F:metal ion binding; IEA:UniProtKB-KW.
            if ('GO;' in line or 'PROSITE' in line or 'Pfam' in line or 'SMART' in line) and line[0:2] == 'DR':
                sp = [x.strip() for x in line.strip('\n').split(';')]
                prositeID, prositeName = sp[1:3]
                prosites[prositeID] = prositeName
            elif line[0:2] == 'KW':
                # KW   NAD; NADP; Oxidoreductase; Phosphopantetheine; Phosphoprotein;
                sp = [x.strip().strip('.') for x in line.strip('\n').split(';')]
                for prosite in sp[1:]:
                    if prosite:
                        prosites['KW:'+prosite] = 'KW:' + prosite
    return prosites

def download(proteinID):
    txtFile = mkTxtFile(proteinID)
    if not os.path.exists(txtFile):
        os.system('wget "http://www.uniprot.org/uniprot/%s.txt" -O %s' % (proteinID, txtFile))
        time.sleep(2)

def getProteins():
    proteins = {}
    with open('../working/neal/proteins.xls') as f:
        for line in f:
            mod, proteinLs = line.strip('\n').split('\t')
            for p in proteinLs.split(';'):
                proteins[p] = True
    return proteins

def updateTxtFiles():
    proteins = getProteins()
    for p in proteins:
        download(p)

def dumpProsite(proteins, prositeFile):
    """ '../working/neal/proteins.prosite' """
    protein2prosite = mkProteinToProsite(proteins)
    with open(prositeFile, 'w') as fout:
        for protein in protein2prosite:
            for prosite in protein2prosite[protein]:
                print >> fout, '\t'.join( (protein, prosite, protein2prosite[protein][prosite]) )

def loadProteinToProsite(prositeFile):
    """ '../working/neal/proteins.prosite' """
    protein2prosite = defaultdict(dict)
    with open(prositeFile) as f:
        for line in f:
            protein, prositeID, prositeName = line.strip('\n').split('\t')
            protein2prosite[protein][prositeID] = prositeName
    return protein2prosite

def loadProteinToGOCC(prositeFile):
    """ '../working/neal/proteins.prosite' """
    protein2prosite = defaultdict(dict)
    with open(prositeFile) as f:
        for line in f:
            protein, prositeID, prositeName = line.strip('\n').split('\t')
            if 'C:' in prositeName:
                protein2prosite[protein][prositeID] = prositeName
    return protein2prosite

def loadGoCounts(proteinFile, goFile):
    """C,P,F counts for experimental proteins per mod
    '../working/neal/proteins.xls'
    '../working/neal/proteins.prosite'
"""
    proteins = {'SNO':defaultdict(dict),
                'RSG':defaultdict(dict),
                'SPAL':defaultdict(dict),
                'SOH':defaultdict(dict)}
    mod2proteins = defaultdict(dict)
    with open(proteinFile) as f:
        for line in f:
            mod, proteinLs = line.strip('\n').split('\t')
            for p in proteinLs.split(';'):
                mod2proteins[mod][p] = proteinLs
    with open(goFile) as f:
        for line in f:
            protein, prositeID, prositeName = line.strip('\n').split('\t')
            goType = prositeName.split(':')[0]
            if goType in ('C', 'F', 'P'):
                for mod in mod2proteins:
                    if protein in mod2proteins[mod]:
                        proteins[mod][goType][mod2proteins[mod][protein]] = True
    return proteins

if __name__ == "__main__":
    dumpProteinNames()
#    dumpProsite()

