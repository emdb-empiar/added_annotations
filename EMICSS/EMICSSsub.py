#!/usr/bin/env python

#
# Generated Wed Jul 14 21:04:48 2021 by generateDS.py version 2.38.6.
# Python 3.7.6 (default, Dec 30 2019, 19:38:28)  [Clang 11.0.0 (clang-1100.0.33.16)]
#
# Command line options:
#   ('-o', 'EMICSS.py')
#   ('-s', 'EMICSSsub.py')
#
# Command line arguments:
#   EMDB_EMICSS.xsd
#
# Command line:
#   /usr/local/bin/generateDS.py -o "EMICSS.py" -s "EMICSSsub.py" EMDB_EMICSS.xsd
#
# Current working directory (os.getcwd()):
#   EMICSS
#

import os
import sys
from lxml import etree as etree_

import ??? as supermod

def parsexml_(infile, parser=None, **kwargs):
    if parser is None:
        # Use the lxml ElementTree compatible parser so that, e.g.,
        #   we ignore comments.
        parser = etree_.ETCompatXMLParser()
    try:
        if isinstance(infile, os.PathLike):
            infile = os.path.join(infile)
    except AttributeError:
        pass
    doc = etree_.parse(infile, parser=parser, **kwargs)
    return doc

def parsexmlstring_(instring, parser=None, **kwargs):
    if parser is None:
        # Use the lxml ElementTree compatible parser so that, e.g.,
        #   we ignore comments.
        try:
            parser = etree_.ETCompatXMLParser()
        except AttributeError:
            # fallback to xml.etree
            parser = etree_.XMLParser()
    element = etree_.fromstring(instring, parser=parser, **kwargs)
    return element

#
# Globals
#

ExternalEncoding = ''
SaveElementTreeNode = True

#
# Data representation classes
#


class emicssSub(supermod.emicss):
    def __init__(self, emdb_id=None, dbs=None, empiars=None, europe_pmc=None, molecular_weight=None, sample=None, **kwargs_):
        super(emicssSub, self).__init__(emdb_id, dbs, empiars, europe_pmc, molecular_weight, sample,  **kwargs_)
supermod.emicss.subclass = emicssSub
# end class emicssSub


class dbsTypeSub(supermod.dbsType):
    def __init__(self, db=None, **kwargs_):
        super(dbsTypeSub, self).__init__(db,  **kwargs_)
supermod.dbsType.subclass = dbsTypeSub
# end class dbsTypeSub


class dbTypeSub(supermod.dbType):
    def __init__(self, db_source=None, db_version=None, **kwargs_):
        super(dbTypeSub, self).__init__(db_source, db_version,  **kwargs_)
supermod.dbType.subclass = dbTypeSub
# end class dbTypeSub


class empiarsTypeSub(supermod.empiarsType):
    def __init__(self, empiar=None, **kwargs_):
        super(empiarsTypeSub, self).__init__(empiar,  **kwargs_)
supermod.empiarsType.subclass = empiarsTypeSub
# end class empiarsTypeSub


class empiarTypeSub(supermod.empiarType):
    def __init__(self, empiar_id=None, provenance=None, **kwargs_):
        super(empiarTypeSub, self).__init__(empiar_id, provenance,  **kwargs_)
supermod.empiarType.subclass = empiarTypeSub
# end class empiarTypeSub


class europe_pmcTypeSub(supermod.europe_pmcType):
    def __init__(self, pmc=None, **kwargs_):
        super(europe_pmcTypeSub, self).__init__(pmc,  **kwargs_)
supermod.europe_pmcType.subclass = europe_pmcTypeSub
# end class europe_pmcTypeSub


class pmcTypeSub(supermod.pmcType):
    def __init__(self, pubmed_id=None, doi=None, issn=None, provenance=None, **kwargs_):
        super(pmcTypeSub, self).__init__(pubmed_id, doi, issn, provenance,  **kwargs_)
supermod.pmcType.subclass = pmcTypeSub
# end class pmcTypeSub


class molecular_weightTypeSub(supermod.molecular_weightType):
    def __init__(self, models=None, weights=None, **kwargs_):
        super(molecular_weightTypeSub, self).__init__(models, weights,  **kwargs_)
supermod.molecular_weightType.subclass = molecular_weightTypeSub
# end class molecular_weightTypeSub


class modelsTypeSub(supermod.modelsType):
    def __init__(self, model=None, **kwargs_):
        super(modelsTypeSub, self).__init__(model,  **kwargs_)
supermod.modelsType.subclass = modelsTypeSub
# end class modelsTypeSub


class modelTypeSub(supermod.modelType):
    def __init__(self, pdb_id=None, assemblies=None, weight=None, units=None, provenance=None, **kwargs_):
        super(modelTypeSub, self).__init__(pdb_id, assemblies, weight, units, provenance,  **kwargs_)
supermod.modelType.subclass = modelTypeSub
# end class modelTypeSub


class weightsTypeSub(supermod.weightsType):
    def __init__(self, weight=None, **kwargs_):
        super(weightsTypeSub, self).__init__(weight,  **kwargs_)
supermod.weightsType.subclass = weightsTypeSub
# end class weightsTypeSub


class weightTypeSub(supermod.weightType):
    def __init__(self, kind=None, method=None, weight=None, units=None, provenance=None, **kwargs_):
        super(weightTypeSub, self).__init__(kind, method, weight, units, provenance,  **kwargs_)
supermod.weightType.subclass = weightTypeSub
# end class weightTypeSub


class sampleTypeSub(supermod.sampleType):
    def __init__(self, name=None, cross_ref_dbs=None, supramolecules=None, macromolecules=None, **kwargs_):
        super(sampleTypeSub, self).__init__(name, cross_ref_dbs, supramolecules, macromolecules,  **kwargs_)
supermod.sampleType.subclass = sampleTypeSub
# end class sampleTypeSub


class cross_ref_dbsTypeSub(supermod.cross_ref_dbsType):
    def __init__(self, cross_ref_db=None, **kwargs_):
        super(cross_ref_dbsTypeSub, self).__init__(cross_ref_db,  **kwargs_)
supermod.cross_ref_dbsType.subclass = cross_ref_dbsTypeSub
# end class cross_ref_dbsTypeSub


class cross_ref_dbTypeSub(supermod.cross_ref_dbType):
    def __init__(self, name=None, db_source=None, provenance=None, db_accession_id=None, score=None, **kwargs_):
        super(cross_ref_dbTypeSub, self).__init__(name, db_source, provenance, db_accession_id, score,  **kwargs_)
supermod.cross_ref_dbType.subclass = cross_ref_dbTypeSub
# end class cross_ref_dbTypeSub


class supramoleculesTypeSub(supermod.supramoleculesType):
    def __init__(self, supramolecule=None, **kwargs_):
        super(supramoleculesTypeSub, self).__init__(supramolecule,  **kwargs_)
supermod.supramoleculesType.subclass = supramoleculesTypeSub
# end class supramoleculesTypeSub


class supramoleculeTypeSub(supermod.supramoleculeType):
    def __init__(self, kind=None, id=None, copies=None, name=None, cross_ref_dbs=None, **kwargs_):
        super(supramoleculeTypeSub, self).__init__(kind, id, copies, name, cross_ref_dbs,  **kwargs_)
supermod.supramoleculeType.subclass = supramoleculeTypeSub
# end class supramoleculeTypeSub


class cross_ref_dbsType1Sub(supermod.cross_ref_dbsType1):
    def __init__(self, cross_ref_db=None, **kwargs_):
        super(cross_ref_dbsType1Sub, self).__init__(cross_ref_db,  **kwargs_)
supermod.cross_ref_dbsType1.subclass = cross_ref_dbsType1Sub
# end class cross_ref_dbsType1Sub


class cross_ref_dbType2Sub(supermod.cross_ref_dbType2):
    def __init__(self, name=None, db_source=None, provenance=None, db_accession_id=None, score=None, **kwargs_):
        super(cross_ref_dbType2Sub, self).__init__(name, db_source, provenance, db_accession_id, score,  **kwargs_)
supermod.cross_ref_dbType2.subclass = cross_ref_dbType2Sub
# end class cross_ref_dbType2Sub


class macromoleculesTypeSub(supermod.macromoleculesType):
    def __init__(self, macromolecule=None, **kwargs_):
        super(macromoleculesTypeSub, self).__init__(macromolecule,  **kwargs_)
supermod.macromoleculesType.subclass = macromoleculesTypeSub
# end class macromoleculesTypeSub


class macromoleculeTypeSub(supermod.macromoleculeType):
    def __init__(self, kind=None, id=None, copies=None, name=None, ccd_id=None, cross_ref_dbs=None, **kwargs_):
        super(macromoleculeTypeSub, self).__init__(kind, id, copies, name, ccd_id, cross_ref_dbs,  **kwargs_)
supermod.macromoleculeType.subclass = macromoleculeTypeSub
# end class macromoleculeTypeSub


class cross_ref_dbsType3Sub(supermod.cross_ref_dbsType3):
    def __init__(self, cross_ref_db=None, **kwargs_):
        super(cross_ref_dbsType3Sub, self).__init__(cross_ref_db,  **kwargs_)
supermod.cross_ref_dbsType3.subclass = cross_ref_dbsType3Sub
# end class cross_ref_dbsType3Sub


class cross_ref_dbType4Sub(supermod.cross_ref_dbType4):
    def __init__(self, name=None, db_source=None, provenance=None, db_accession_id=None, score=None, **kwargs_):
        super(cross_ref_dbType4Sub, self).__init__(name, db_source, provenance, db_accession_id, score,  **kwargs_)
supermod.cross_ref_dbType4.subclass = cross_ref_dbType4Sub
# end class cross_ref_dbType4Sub


def get_root_tag(node):
    tag = supermod.Tag_pattern_.match(node.tag).groups()[-1]
    rootClass = None
    rootClass = supermod.GDSClassesMapping.get(tag)
    if rootClass is None and hasattr(supermod, tag):
        rootClass = getattr(supermod, tag)
    return tag, rootClass


def parse(inFilename, silence=False):
    parser = None
    doc = parsexml_(inFilename, parser)
    rootNode = doc.getroot()
    rootTag, rootClass = get_root_tag(rootNode)
    if rootClass is None:
        rootTag = 'emicss'
        rootClass = supermod.emicss
    rootObj = rootClass.factory()
    rootObj.build(rootNode)
    # Enable Python to collect the space used by the DOM.
    if not SaveElementTreeNode:
        doc = None
        rootNode = None
    if not silence:
        sys.stdout.write('<?xml version="1.0" ?>\n')
        rootObj.export(
            sys.stdout, 0, name_=rootTag,
            namespacedef_='',
            pretty_print=True)
    return rootObj


def parseEtree(inFilename, silence=False):
    parser = None
    doc = parsexml_(inFilename, parser)
    rootNode = doc.getroot()
    rootTag, rootClass = get_root_tag(rootNode)
    if rootClass is None:
        rootTag = 'emicss'
        rootClass = supermod.emicss
    rootObj = rootClass.factory()
    rootObj.build(rootNode)
    mapping = {}
    rootElement = rootObj.to_etree(None, name_=rootTag, mapping_=mapping)
    reverse_mapping = rootObj.gds_reverse_node_mapping(mapping)
    # Enable Python to collect the space used by the DOM.
    if not SaveElementTreeNode:
        doc = None
        rootNode = None
    if not silence:
        content = etree_.tostring(
            rootElement, pretty_print=True,
            xml_declaration=True, encoding="utf-8")
        sys.stdout.write(content)
        sys.stdout.write('\n')
    return rootObj, rootElement, mapping, reverse_mapping


def parseString(inString, silence=False):
    if sys.version_info.major == 2:
        from StringIO import StringIO
    else:
        from io import BytesIO as StringIO
    parser = None
    rootNode= parsexmlstring_(inString, parser)
    rootTag, rootClass = get_root_tag(rootNode)
    if rootClass is None:
        rootTag = 'emicss'
        rootClass = supermod.emicss
    rootObj = rootClass.factory()
    rootObj.build(rootNode)
    # Enable Python to collect the space used by the DOM.
    if not SaveElementTreeNode:
        rootNode = None
    if not silence:
        sys.stdout.write('<?xml version="1.0" ?>\n')
        rootObj.export(
            sys.stdout, 0, name_=rootTag,
            namespacedef_='')
    return rootObj


def parseLiteral(inFilename, silence=False):
    parser = None
    doc = parsexml_(inFilename, parser)
    rootNode = doc.getroot()
    rootTag, rootClass = get_root_tag(rootNode)
    if rootClass is None:
        rootTag = 'emicss'
        rootClass = supermod.emicss
    rootObj = rootClass.factory()
    rootObj.build(rootNode)
    # Enable Python to collect the space used by the DOM.
    if not SaveElementTreeNode:
        doc = None
        rootNode = None
    if not silence:
        sys.stdout.write('#from ??? import *\n\n')
        sys.stdout.write('import ??? as model_\n\n')
        sys.stdout.write('rootObj = model_.rootClass(\n')
        rootObj.exportLiteral(sys.stdout, 0, name_=rootTag)
        sys.stdout.write(')\n')
    return rootObj


USAGE_TEXT = """
Usage: python ???.py <infilename>
"""


def usage():
    print(USAGE_TEXT)
    sys.exit(1)


def main():
    args = sys.argv[1:]
    if len(args) != 1:
        usage()
    infilename = args[0]
    parse(infilename)


if __name__ == '__main__':
    #import pdb; pdb.set_trace()
    main()
