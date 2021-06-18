#!/usr/bin/env python

#
# Generated Fri Jun 18 10:33:10 2021 by generateDS.py version 2.38.6.
# Python 3.7.6 (default, Dec 30 2019, 19:38:28)  [Clang 11.0.0 (clang-1100.0.33.16)]
#
# Command line options:
#   ('-o', 'EICSS.py')
#   ('-s', 'EICSSsub.py')
#
# Command line arguments:
#   EMDB_EICSS.xsd
#
# Command line:
#   /usr/local/bin/generateDS.py -o "EICSS.py" -s "EICSSsub.py" EMDB_EICSS.xsd
#
# Current working directory (os.getcwd()):
#   EICSS
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


class eicssSub(supermod.eicss):
    def __init__(self, EMDB_ID=None, DBs_list=None, models_list=None, sample_annotation=None, **kwargs_):
        super(eicssSub, self).__init__(EMDB_ID, DBs_list, models_list, sample_annotation,  **kwargs_)
supermod.eicss.subclass = eicssSub
# end class eicssSub


class DBs_listTypeSub(supermod.DBs_listType):
    def __init__(self, DB=None, **kwargs_):
        super(DBs_listTypeSub, self).__init__(DB,  **kwargs_)
supermod.DBs_listType.subclass = DBs_listTypeSub
# end class DBs_listTypeSub


class DBTypeSub(supermod.DBType):
    def __init__(self, DB_source=None, DB_version=None, **kwargs_):
        super(DBTypeSub, self).__init__(DB_source, DB_version,  **kwargs_)
supermod.DBType.subclass = DBTypeSub
# end class DBTypeSub


class models_listTypeSub(supermod.models_listType):
    def __init__(self, model_annotation=None, **kwargs_):
        super(models_listTypeSub, self).__init__(model_annotation,  **kwargs_)
supermod.models_listType.subclass = models_listTypeSub
# end class models_listTypeSub


class model_annotationTypeSub(supermod.model_annotationType):
    def __init__(self, PDBID=None, assemblies=None, weight=None, units=None, provenance=None, **kwargs_):
        super(model_annotationTypeSub, self).__init__(PDBID, assemblies, weight, units, provenance,  **kwargs_)
supermod.model_annotationType.subclass = model_annotationTypeSub
# end class model_annotationTypeSub


class sample_annotationTypeSub(supermod.sample_annotationType):
    def __init__(self, sample_name=None, list_crossRefDBs=None, list_supra_molecules=None, list_macro_molecules=None, **kwargs_):
        super(sample_annotationTypeSub, self).__init__(sample_name, list_crossRefDBs, list_supra_molecules, list_macro_molecules,  **kwargs_)
supermod.sample_annotationType.subclass = sample_annotationTypeSub
# end class sample_annotationTypeSub


class list_crossRefDBsTypeSub(supermod.list_crossRefDBsType):
    def __init__(self, crossRefDB=None, **kwargs_):
        super(list_crossRefDBsTypeSub, self).__init__(crossRefDB,  **kwargs_)
supermod.list_crossRefDBsType.subclass = list_crossRefDBsTypeSub
# end class list_crossRefDBsTypeSub


class crossRefDBTypeSub(supermod.crossRefDBType):
    def __init__(self, DB_source=None, provenance=None, DB_accession_ID=None, score=None, **kwargs_):
        super(crossRefDBTypeSub, self).__init__(DB_source, provenance, DB_accession_ID, score,  **kwargs_)
supermod.crossRefDBType.subclass = crossRefDBTypeSub
# end class crossRefDBTypeSub


class list_supra_moleculesTypeSub(supermod.list_supra_moleculesType):
    def __init__(self, supra_molecule_annotation=None, **kwargs_):
        super(list_supra_moleculesTypeSub, self).__init__(supra_molecule_annotation,  **kwargs_)
supermod.list_supra_moleculesType.subclass = list_supra_moleculesTypeSub
# end class list_supra_moleculesTypeSub


class supra_molecule_annotationTypeSub(supermod.supra_molecule_annotationType):
    def __init__(self, supra_kind=None, supra_ID=None, supra_copies=None, supra_name=None, list_crossRefDBs=None, **kwargs_):
        super(supra_molecule_annotationTypeSub, self).__init__(supra_kind, supra_ID, supra_copies, supra_name, list_crossRefDBs,  **kwargs_)
supermod.supra_molecule_annotationType.subclass = supra_molecule_annotationTypeSub
# end class supra_molecule_annotationTypeSub


class list_crossRefDBsType1Sub(supermod.list_crossRefDBsType1):
    def __init__(self, crossRefDB=None, **kwargs_):
        super(list_crossRefDBsType1Sub, self).__init__(crossRefDB,  **kwargs_)
supermod.list_crossRefDBsType1.subclass = list_crossRefDBsType1Sub
# end class list_crossRefDBsType1Sub


class crossRefDBType2Sub(supermod.crossRefDBType2):
    def __init__(self, DB_source=None, provenance=None, DB_accession_ID=None, score=None, **kwargs_):
        super(crossRefDBType2Sub, self).__init__(DB_source, provenance, DB_accession_ID, score,  **kwargs_)
supermod.crossRefDBType2.subclass = crossRefDBType2Sub
# end class crossRefDBType2Sub


class list_macro_moleculesTypeSub(supermod.list_macro_moleculesType):
    def __init__(self, macro_molecule_annotation=None, **kwargs_):
        super(list_macro_moleculesTypeSub, self).__init__(macro_molecule_annotation,  **kwargs_)
supermod.list_macro_moleculesType.subclass = list_macro_moleculesTypeSub
# end class list_macro_moleculesTypeSub


class macro_molecule_annotationTypeSub(supermod.macro_molecule_annotationType):
    def __init__(self, macro_kind=None, macro_ID=None, macro_copies=None, macro_name=None, macro_CCD_ID=None, list_crossRefDBs=None, **kwargs_):
        super(macro_molecule_annotationTypeSub, self).__init__(macro_kind, macro_ID, macro_copies, macro_name, macro_CCD_ID, list_crossRefDBs,  **kwargs_)
supermod.macro_molecule_annotationType.subclass = macro_molecule_annotationTypeSub
# end class macro_molecule_annotationTypeSub


class list_crossRefDBsType3Sub(supermod.list_crossRefDBsType3):
    def __init__(self, crossRefDB=None, **kwargs_):
        super(list_crossRefDBsType3Sub, self).__init__(crossRefDB,  **kwargs_)
supermod.list_crossRefDBsType3.subclass = list_crossRefDBsType3Sub
# end class list_crossRefDBsType3Sub


class crossRefDBType4Sub(supermod.crossRefDBType4):
    def __init__(self, DB_source=None, provenance=None, DB_accession_ID=None, score=None, **kwargs_):
        super(crossRefDBType4Sub, self).__init__(DB_source, provenance, DB_accession_ID, score,  **kwargs_)
supermod.crossRefDBType4.subclass = crossRefDBType4Sub
# end class crossRefDBType4Sub


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
        rootTag = 'eicss'
        rootClass = supermod.eicss
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
        rootTag = 'eicss'
        rootClass = supermod.eicss
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
        rootTag = 'eicss'
        rootClass = supermod.eicss
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
        rootTag = 'eicss'
        rootClass = supermod.eicss
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
