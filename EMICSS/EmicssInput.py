import itertools

class EmicssInput:
    """
    Writing all the annotations as dictionary for generateDS input for every entry
    """

    def __init__(self, mapping_list):
        self.mapping_list = mapping_list

    def execute(self):
        emicss_annotation = self.dict_emicss(self.mapping_list)
        return emicss_annotation

    def dict_emicss(self, mapping_list):
        """
        Converts dictionary individual annotation to deeply nested dictionary of all the added annotations
        """

        emicss_dict = {}

        for db in range(0, len(mapping_list), 2):
            if mapping_list[db] == "UNIPROT":
                self.unip_map = mapping_list[db+1]
            if mapping_list[db] == "COMPLEX":
                self.cpx_map = mapping_list[db + 1]
            if mapping_list[db] == "LIGANDS":
                self.lig_map = mapping_list[db+1]
            if mapping_list[db] == "MODEL":
                self.mw_map = mapping_list[db+1]
            if mapping_list[db] == "WEIGHT":
                self.sw_map = mapping_list[db + 1]
            if mapping_list[db] == "EMPIAR":
                self.empiar_map = mapping_list[db+1]
            if mapping_list[db] == "CITATION":
                self.pmc_map = mapping_list[db+1]
            if mapping_list[db] == "PROTEIN-TERMS":
                self.proteins_map = mapping_list[db+1]
            if mapping_list[db] == "PDBeKB":
                self.KB_map = mapping_list[db+1]
            if mapping_list[db] == "ALPHAFOLD":
                self.AF_map = mapping_list[db+1]

        try:
            if self.mw_map:
                for mw in self.mw_map:
                    if mw.emdb_id not in emicss_dict:
                        emicss_dict[mw.emdb_id] = {}
                    if mw.emdb_id not in emicss_dict[mw.emdb_id]:
                        emicss_dict[mw.emdb_id][mw.pdb_id] = mw.__dict__
                    else:
                        emicss_dict[mw.emdb_id][mw.pdb_id] += mw.__dict__
        except AttributeError as e:
            print("MODEL mapping doesn't exist", e)

        try:
            if self.sw_map:
                for sw in self.sw_map:
                    if sw.emdb_id not in emicss_dict:
                        emicss_dict[sw.emdb_id] = {}
                    if sw.emdb_id not in emicss_dict[sw.emdb_id]:
                        emicss_dict[sw.emdb_id][sw.method] = sw.__dict__
                    else:
                        emicss_dict[sw.emdb_id][sw.method] += sw.__dict__
        except AttributeError as e:
            print("WEIGHT mapping doesn't exist", e)

        try:
            if self.empiar_map:
                for empiar in self.empiar_map:
                    if empiar.emdb_id not in emicss_dict:
                        emicss_dict[empiar.emdb_id] = {}
                    if empiar.emdb_id not in emicss_dict[empiar.emdb_id]:
                        emicss_dict[empiar.emdb_id][empiar.empiar_id] = empiar.__dict__
                    else:
                        emicss_dict[empiar.emdb_id][empiar.empiar_id] += empiar.__dict__
        except AttributeError as e:
            print("EMPIAR mapping doesn't exist", e)

        try:
            if self.pmc_map:
                for pmc in self.pmc_map:
                    if pmc.emdb_id not in emicss_dict:
                        emicss_dict[pmc.emdb_id] = {}
                    if pmc.emdb_id not in emicss_dict[pmc.emdb_id]:
                        emicss_dict[pmc.emdb_id]["PMC"] = pmc.__dict__
                    else:
                        emicss_dict[pmc.emdb_id]["PMC"] += pmc.__dict__
        except AttributeError as e:
            print("CITATION mapping doesn't exist", e)

        try:
            if self.unip_map:
                for unip in self.unip_map:
                    if unip.emdb_id not in emicss_dict:
                        emicss_dict[unip.emdb_id] = {}
                    if unip.emdb_id not in emicss_dict[unip.emdb_id]:
                        emicss_dict[unip.emdb_id][unip.uniprot_id] = unip.__dict__
                    else:
                        emicss_dict[unip.emdb_id][unip.uniprot_id] += unip.__dict__
                    try:
                        if self.unip_map and self.proteins_map:
                            for ProTerm in self.proteins_map:
                                if ProTerm.uniprot_id == unip.uniprot_id:
                                    ind = 0
                                    if ProTerm.emdb_id not in emicss_dict.keys():
                                        emicss_dict[ProTerm.emdb_id] = {}
                                    if ProTerm.emdb_id not in emicss_dict[ProTerm.emdb_id].keys():
                                        emicss_dict[ProTerm.emdb_id][ProTerm.uniprot_id] = {}
                                    if ProTerm.sample_id not in emicss_dict[ProTerm.emdb_id]:
                                        emicss_dict[ProTerm.emdb_id][ProTerm.uniprot_id] = unip.__dict__
                                    for go2 in ProTerm.go:
                                        go = go2.__dict__
                                        if ProTerm.uniprot_id == go2.unip_id:
                                            for k in go.keys():
                                                new_key = '{}_{}'.format(k, ind)
                                                emicss_dict[ProTerm.emdb_id][ProTerm.uniprot_id][new_key] = go[k]
                                            ind = ind + 1

                                    for ipr2 in ProTerm.interpro:
                                        ipr = ipr2.__dict__
                                        if ProTerm.uniprot_id == ipr2.unip_id:
                                            for k in ipr.keys():
                                                new_key = '{}_{}_{}'.format("ipr", k, ind)
                                                emicss_dict[ProTerm.emdb_id][ProTerm.uniprot_id][new_key] = ipr[k]
                                            ind = ind + 1

                                    for pfam2 in ProTerm.pfam:
                                        pfam = pfam2.__dict__
                                        if ProTerm.uniprot_id == pfam2.unip_id:
                                            for k in pfam.keys():
                                                new_key = '{}_{}_{}'.format("pfam", k, ind)
                                                emicss_dict[ProTerm.emdb_id][ProTerm.uniprot_id][new_key] = pfam[k]
                                            ind = ind + 1
                                            
                                    for cath2 in ProTerm.cath:
                                        cath = cath2.__dict__
                                        if ProTerm.uniprot_id == cath2.unip_id:
                                            for k in cath.keys():
                                                new_key = '{}_{}_{}'.format("cath", k, ind)
                                                emicss_dict[ProTerm.emdb_id][ProTerm.uniprot_id][new_key] = cath[k]
                                            ind = ind + 1

                                    for scopi in ProTerm.scop:
                                        scop = scopi.__dict__
                                        if ProTerm.uniprot_id == scopi.unip_id:
                                            for k in scop.keys():
                                                new_key = '{}_{}_{}'.format("scop", k, ind)
                                                emicss_dict[ProTerm.emdb_id][ProTerm.uniprot_id][new_key] = scop[k]
                                            ind = ind + 1

                                    for scopj in ProTerm.scop2:
                                        scop2 = scopj.__dict__
                                        if ProTerm.uniprot_id == scopj.unip_id:
                                            for k in scop2.keys():
                                                new_key = '{}_{}_{}'.format("scop2", k, ind)
                                                emicss_dict[ProTerm.emdb_id][ProTerm.uniprot_id][new_key] = scop2[k]
                                            ind = ind + 1

                                    emicss_dict[ProTerm.emdb_id][ProTerm.uniprot_id]["ind"] = ind
                                    self.ind = ind
                    except AttributeError as e:
                        print("PROTEIN-TERMS mapping doesn't exist", e)

                    try:
                        if self.unip_map and self.KB_map:
                            for KB in self.KB_map:
                                if KB.uniprot_id == unip.uniprot_id:
                                    if self.ind:
                                        ind = self.ind
                                    if not self.ind:
                                        ind = 0
                                    for pdbekb2 in KB.pdbekb:
                                        pdbekb = pdbekb2.__dict__
                                        if KB.uniprot_id == pdbekb2.unip_id:
                                            for k in pdbekb.keys():
                                                new_key = '{}_{}_{}'.format("kb", k, ind)
                                                emicss_dict[unip.emdb_id][KB.uniprot_id][new_key] = pdbekb[k]
                                            ind = ind + 1
                                    emicss_dict[unip.emdb_id][KB.uniprot_id]["ind"] = ind
                                    self.ind = ind
                    except AttributeError as e:
                        print("PDBeKB mapping doesn't exist", e)

                    try:
                        if self.unip_map and self.AF_map:
                            for AF in self.AF_map:
                                if AF.uniprot_id == unip.uniprot_id:
                                    if self.ind:
                                        ind = self.ind
                                    if not self.ind:
                                        ind = 0
                                    for alphafold2 in AF.alphafold:
                                        alphafold = alphafold2.__dict__
                                        if AF.uniprot_id == alphafold2.unip_id:
                                            for k in alphafold.keys():
                                                new_key = '{}_{}_{}'.format("af", k, ind)
                                                emicss_dict[AF.emdb_id][AF.uniprot_id][new_key] = alphafold[k]
                                            ind = ind + 1
                                    emicss_dict[AF.emdb_id][AF.uniprot_id]["ind"] = ind
                    except AttributeError as e:
                        print("ALPHAFOLD DB mapping doesn't exist", e)
        except AttributeError as e:
            print("UNIPROT mapping doesn't exist", e)

        try:
            if self.cpx_map:
                for emcpx in self.cpx_map:
                    if emcpx:
                        for cpx in emcpx.cpx_list:
                            if emcpx.emdb_id not in emicss_dict.keys():
                                emicss_dict[emcpx.emdb_id] = {}
                            if emcpx.sample_id not in emicss_dict[emcpx.emdb_id].keys():
                                emicss_dict[emcpx.emdb_id][emcpx.sample_id] = {}
                                ind = 0
                            lcpx = ["supra_name", emcpx.supra_name, "sample_copies", emcpx.sample_copies, "cpx_id" + "_" + str(ind),
                                    cpx.cpx_id, "cpx_name" + "_" + str(ind), cpx.name, "provenance" + "_" + str(ind),
                                    emcpx.provenance, "score" + "_" + str(ind), emcpx.score]
                            dcpx = dict(itertools.zip_longest(*[iter(lcpx)] * 2, fillvalue=""))
                            for k in dcpx.keys():
                                emicss_dict[emcpx.emdb_id][emcpx.sample_id][k] = dcpx[k]
                            ind = ind + 1
                        emicss_dict[emcpx.emdb_id][emcpx.sample_id]["ind"] = ind
        except AttributeError as e:
            print("COMPLEX mapping doesn't exist", e)

        try:
            if self.lig_map:
                for ligand in self.lig_map:
                    if ligand.emdb_id not in emicss_dict:
                        emicss_dict[ligand.emdb_id] = {}
                    if ligand.emdb_id not in emicss_dict[ligand.emdb_id]:
                        emicss_dict[ligand.emdb_id]["ligand_"+ligand.sample_id] = ligand.__dict__
                    else:
                        emicss_dict[ligand.emdb_id]["ligand_"+ligand.sample_id] += ligand.__dict__
        except AttributeError as e:
            print("LIGAND mapping doesn't exist", e)

        return emicss_dict

