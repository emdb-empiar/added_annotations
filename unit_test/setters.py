from models import Protein, GO, Interpro, Cath, Pdbekb, Supra, EMDB_complex, Ligand

def set_protein(emdb_id, sample_id, sample_name, sample_organism, pdb, sample_complexes, uniprot_id, provenance,
            sequence, sample_copies, go, interpro, pfam, cath, scop, scop2, scop2B, pdbekb, alphafold):
    prot = Protein(emdb_id, sample_id)
    prot.emdb_id = emdb_id
    prot.sample_id = sample_id
    prot.sample_name = sample_name
    prot.sample_organism = sample_organism
    prot.pdb = pdb
    prot.sample_complexes = sample_complexes
    prot.uniprot_id = uniprot_id
    prot.provenance = provenance
    prot.sequence = sequence
    prot.sample_copies = sample_copies
    prot.go = go
    prot.interpro = interpro
    prot.pfam = pfam
    prot.cath = cath
    prot.scop = scop
    prot.scop2 = scop2
    prot.scop2B = scop2B
    prot.pdbekb = pdbekb
    prot.alphafold = alphafold
    return prot


def set_go(id, namespace, type, unip_id, provenance):
    go = GO()
    go.id = id
    go.namespace = namespace
    go.type = type
    go.unip_id = unip_id
    go.provenance = provenance
    return go


def set_ipro(id, namespace, unip_id, provenance, start, end, unp_start, unp_end):
    intpfam = Interpro()
    intpfam.id = id
    intpfam.namespace = namespace
    intpfam.unip_id = unip_id
    intpfam.provenance = provenance
    intpfam.start = start
    intpfam.end = end
    intpfam.unp_start = unp_start
    intpfam.unp_end = unp_end
    return intpfam


def set_csss(id, unip_id, provenance, start, end, unp_start, unp_end):
    csss = Cath()
    csss.id = id
    csss.unip_id = unip_id
    csss.provenance = provenance
    csss.start = start
    csss.end = end
    csss.unp_start = unp_start
    csss.unp_end = unp_end
    return csss


def set_kbalpha(unip_id, provenance):
    kbalpha = Pdbekb(unip_id, provenance)
    kbalpha.unip_id = unip_id
    kbalpha.provenance = provenance
    return kbalpha


def set_supra(emdb_id, supra_id, supra_name, kind):
    supra = Supra(emdb_id, supra_id)
    supra.emdb_id = emdb_id
    supra.supra_id = supra_id
    supra.supra_name = supra_name
    supra.kind = kind
    return supra


def set_EMDB_comp(emdb_id, sample_id, supra_name, sample_copies, complex_sample_id, cpx_list, proteins, provenance, score):
    em_cpx = EMDB_complex(emdb_id, sample_id, supra_name, sample_copies, complex_sample_id)
    em_cpx.emdb_id = emdb_id
    em_cpx.sample_id = sample_id
    em_cpx.supra_name = supra_name
    em_cpx.sample_copies = sample_copies
    em_cpx.complex_sample_id = complex_sample_id
    em_cpx.cpx_list = cpx_list
    em_cpx.proteins = proteins
    em_cpx.provenance = provenance
    em_cpx.score = score
    return em_cpx

def set_ligand(emdb_id, sample_id, HET, lig_name, lig_copies, chembl_id, chebi_id, drugbank_id, provenance_chembl, provenance_chebi, provenance_drugbank):
    lig = Ligand(emdb_id, sample_id)
    lig.HET = HET
    lig.lig_name = lig_name
    lig.lig_copies = lig_copies
    lig.chembl_id = chembl_id
    lig.chebi_id = chebi_id
    lig.drugbank_id = drugbank_id
    lig.provenance_chembl = provenance_chembl
    lig.provenance_chebi = provenance_chebi
    lig.provenance_drugbank = provenance_drugbank
    return lig
