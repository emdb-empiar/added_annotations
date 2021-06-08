"""
empiar_headers_loader.py
Load EMDB database. Used in the web-based release procedure.
Loads only one entry
Copyright [2020] EMBL - European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the
"License"); you may not use this file except in
compliance with the License. You may obtain a copy of
the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on
an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
"""

__author__ = 'Andrii Iudin'
__email__ = 'andrii@ebi.ac.uk'
__date__ = '2015-01-23'


import datetime
import logging
import math
import os
import os.path
import sys
import time
import traceback
from builtins import int
from random import random

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
if sys.version_info.major == 2:
    BaseStrType_ = basestring
else:
    BaseStrType_ = str

from empiar import empiar
from django.contrib.auth.models import User
from django.db import IntegrityError
from django.db.utils import OperationalError
from six import iteritems
from empiar_da.forms import category_choices, header_format_choices, data_format_choices, voxel_type_choices
from empiar_da.models import EmpiarAnnotEntry
from empiar_da.models import (EmpiarPrInvCorrAuthor, EmpiarCorrespondingAuthor, EmpiarPrincipalInvestigator,
                              EmpiarEntry, EmpiarImageSet, EmpiarSegmentation, EmpiarAuthorOrcid, EmpiarAuthor,
                              EmpiarCrossReference, EmpiarBioStudiesReference, EmpiarIdrReference,
                              EmpiarEmpiarReference, EmpiarVersionHistory, EmpiarCitation,
                              EmpiarCitAuthor, EmpiarCitEditor, EmpiarScipion)
from empiar_da.models import (EmpiarPreRelEntryRights, EmpiarPreRelCorrespondingAuthor,
                              EmpiarPreRelPrincipalInvestigator, EmpiarPreRelEntry, EmpiarPreRelImageSet,
                              EmpiarPreRelSegmentation, EmpiarPreRelAuthor, EmpiarPreRelCrossReference,
                              EmpiarPreRelBioStudiesReference, EmpiarPreRelIdrReference, EmpiarPreRelEmpiarReference,
                              EmpiarPreRelCitation, EmpiarPreRelCitAuthor, EmpiarPreRelCitEditor, EmpiarPreRelScipion,
                              unsorted_country_choices,
                              experiment_type_choices)
from empiar_da.util import file_exists, model_get_or_create_from_id, read_file

EMPIAR_NAMESPACES = {'e': 'EMPIARSchema'}
NS = EMPIAR_NAMESPACES

# Write into the server log
logger = logging.getLogger(__name__)


def xstr(s):
    return '' if s is None else str(s)


class EMPIARLoader:
    """
    Handles loading of EMPIAR database from XML files
    """

    def __init__(self):
        self.prerel = False  # Set to true if the header is being used to create a new deposition, rather than to
        # populate the released entry table

        # Log that goes to the annotator/web interface
        self.log = StringIO()

    def _proper_country(self, country):
        result = ''
        for unsorted_country in unsorted_country_choices:
            if country == unsorted_country[1]:
                result = unsorted_country[0]
            elif country == 'USA':
                result = 'US'
            elif country == 'UK':
                result = 'GB'
        return result

    def _model_get_or_create_auto_id(self, model_class, field_name_value):
        """
        Query a database for an item with a field=value from field_name_value dictionary.
        If it does not exist, then create one.
        @param model_class: name of database model class
        @param field_name_value: dictionary of model class fields as keys and values of fields
        as corresponding values of keys
        @return: model object with the value set
        """
        # field = model_class._meta.get_field_by_name(field_name)
        try:
            ref = model_class.objects.get(**field_name_value)
        except model_class.DoesNotExist:
            num_retries = 6
            operation = 0
            retry = True
            while retry and operation < num_retries:
                try:
                    if len(model_class.objects.all()) > 0:
                        field_name_value['id'] = int(model_class.objects.latest('id').id) + 1
                    else:
                        field_name_value['id'] = 1

                    ref = model_class(**field_name_value)
                    ref.save()
                    retry = False
                except OperationalError:
                    sleep_time = math.pow(2, operation) + random()
                    time.sleep(sleep_time)

                operation += 1

        return ref

    def _load_cross_references(self, header_xml, entry):
        """
        Load list of cross references from XML file according to schema
        @param header_xml: parsed XML tree from which cross references are extracted
        @param entry: related entry - ManyToMany relationship with EmpiarEntry table
        """

        if not self.prerel:
            existingCrossReferencesList = entry.empiarcrossreference_set.all()
            if existingCrossReferencesList:
                for existingCrossReference in existingCrossReferencesList:
                    existingCrossReference.delete()

            existingCitationsList = entry.empiarcitation_set.all()
            if existingCitationsList:
                for existingCitation in existingCitationsList:
                    existingCitation.empiarcitauthor_set.all().delete()
                    existingCitation.empiarciteditor_set.all().delete()
                    existingCitation.delete()

        # Extract the information about cross references:
        crossReferences = header_xml.get_crossReferences()
        if crossReferences:
            try:
                relatedEMDBEntries = crossReferences.get_relatedEMDBEntries().get_emdbEntry()
                for relatedEMDBEntry in relatedEMDBEntries:
                    if self.prerel:
                        self._model_get_or_create_auto_id(EmpiarPreRelCrossReference,
                                                          {'name': relatedEMDBEntry, 'entry': entry})
                    else:
                        new_cross_reference = self._model_get_or_create_auto_id(EmpiarCrossReference,
                                                                                {'name': relatedEMDBEntry})
                        new_cross_reference.entry.add(entry)
                        new_cross_reference.save()
            except AttributeError:
                pass

            try:
                relatedBioStudiesEntries = crossReferences.get_relatedBioStudiesEntries().get_bioStudiesEntry()
                for relatedBioStudiesEntry in relatedBioStudiesEntries:
                    if self.prerel:
                        self._model_get_or_create_auto_id(EmpiarPreRelBioStudiesReference,
                                                          {'name': relatedBioStudiesEntry, 'entry': entry})
                    else:
                        new_biostudies_reference = self._model_get_or_create_auto_id(EmpiarBioStudiesReference,
                                                                                     {'name': relatedBioStudiesEntry})
                        new_biostudies_reference.entry.add(entry)
                        new_biostudies_reference.save()
            except AttributeError:
                pass

            try:
                relatedIdrEntries = crossReferences.get_relatedIdrEntries().get_idrEntry()
                for relatedIdrEntry in relatedIdrEntries:
                    if self.prerel:
                        self._model_get_or_create_auto_id(EmpiarPreRelIdrReference,
                                                          {'name': relatedIdrEntry, 'entry': entry})
                    else:
                        new_idr_reference = self._model_get_or_create_auto_id(EmpiarIdrReference,
                                                                                     {'name': relatedIdrEntry})
                        new_idr_reference.entry.add(entry)
                        new_idr_reference.save()
            except AttributeError:
                pass

            try:
                relatedEmpiarEntries = crossReferences.get_relatedEmpiarEntries().get_empiarEntry()
                for relatedEmpiarEntry in relatedEmpiarEntries:
                    if self.prerel:
                        self._model_get_or_create_auto_id(EmpiarPreRelEmpiarReference,
                                                          {'name': relatedEmpiarEntry, 'entry': entry})
                    else:
                        new_empiar_reference = self._model_get_or_create_auto_id(EmpiarEmpiarReference,
                                                                                     {'name': relatedEmpiarEntry})
                        new_empiar_reference.entry.add(entry)
                        new_empiar_reference.save()
            except AttributeError:
                pass

            cit_list = crossReferences.get_citationList()
            if cit_list:
                universalCitationList = cit_list.get_universalCitation()

                for universalCitation in universalCitationList:
                    citation = universalCitation.get_citationType()

                    citation_dict = {}

                    if citation.original_tagname_ == 'journalCitation':
                        citation_dict['title'] = citation.get_title()
                        citation_dict['journal'] = citation.get_journal()
                        citation_dict['journal_abbreviation'] = citation.get_journalAbbreviation()
                        citation_dict['issue'] = citation.get_issue()
                        citation_dict['j_or_nj_citation'] = True
                        citation_dict['preprint'] = citation.get_preprint()
                    elif citation.original_tagname_ == 'nonJournalCitation':
                        citation_dict['title'] = citation.get_bookTitle()
                        citation_dict['publisher'] = citation.get_publisher()
                        citation_dict['publication_location'] = citation.get_publicationLocation()
                        citation_dict['j_or_nj_citation'] = False

                    citation_dict['published'] = citation.get_published()
                    citation_dict['country'] = citation.get_country()
                    citation_dict['volume'] = citation.get_volume()
                    citation_dict['first_page'] = citation.get_firstPage()
                    citation_dict['last_page'] = citation.get_lastPage()
                    citation_dict['year'] = citation.get_year()
                    citation_dict['language'] = citation.get_language()
                    citation_dict['details'] = citation.get_details()

                    externalReferences = citation.get_externalReferences()

                    if externalReferences:
                        for externalReference in externalReferences:

                            referenceType = externalReference.get_type()
                            if referenceType == 'doi':
                                citation_dict['doi'] = externalReference.get_valueOf_()
                            elif referenceType == 'pubmed':
                                citation_dict['pubmedid'] = externalReference.get_valueOf_()

                    citation_dict['entry'] = entry

                    if self.prerel:
                        empiar_citation = self._model_get_or_create_auto_id(EmpiarPreRelCitation, citation_dict)
                    else:
                        empiar_citation = self._model_get_or_create_auto_id(EmpiarCitation, citation_dict)

                    if self.prerel:
                        citAuthorClass = EmpiarPreRelCitAuthor
                        citEditorClass = EmpiarPreRelCitEditor
                    else:
                        citAuthorClass = EmpiarCitAuthor
                        citEditorClass = EmpiarCitEditor

                    authors = citation.get_author()
                    for author in authors:
                        author_dict = {}
                        if not self.prerel:
                            author_dict['name'] = author.get_valueOf_()
                        else:
                            name = author.get_valueOf_()
                            if ' ' in name:
                                author_dict['name'] = "(u'%s', u'%s')" % tuple(author.get_valueOf_().rsplit(" ", 1))
                            else:
                                author_dict['name'] = "(u'" + name + "', u'')"
                        author_dict['author_orcid'] = author.get_authorORCID()
                        author_dict['publication'] = empiar_citation
                        author_dict['order_id'] = author.get_order()
                        self._model_get_or_create_auto_id(citAuthorClass, author_dict)

                    editors = citation.get_editor()
                    for editor in editors:
                        editor_dict = {}
                        if not self.prerel:
                            editor_dict['name'] = editor.get_valueOf_()
                        else:
                            name = editor.get_valueOf_()
                            if ' ' in name:
                                editor_dict['name'] = "(u'%s', u'%s')" % tuple(editor.get_valueOf_().rsplit(" ", 1))
                            else:
                                editor_dict['name'] = "(u'" + name + "', u'')"
                        editor_dict['author_orcid'] = editor.get_authorORCID()
                        editor_dict['publication'] = empiar_citation
                        editor_dict['order_id'] = editor.get_order()
                        self._model_get_or_create_auto_id(citEditorClass, editor_dict)

    def _load_authors(self, admin, entry):
        """
        Load list of authors from XML file according to schema
        
        @param admin: admin part of the parsed XML tree from which authors are extracted
        @param entry: related entry - foreign key to EmpiarEntry table
        """

        # Extract the information about authors:
        authorList = admin.get_authorsList().get_author()

        if not self.prerel:
            existingAuthorList = entry.empiarauthor_set.all().delete()

        i = 0
        for author in authorList:
            if self.prerel:
                name = author.get_valueOf_()
                if ' ' in name:
                    name_from_tuple = "(u'%s', u'%s')" % tuple(author.get_valueOf_().rsplit(" ", 1))
                else:
                    name_from_tuple = "(u'" + name + "', u'')"
                self._model_get_or_create_auto_id(EmpiarPreRelAuthor,
                                                  {'entry': entry, 'order_id': i, 'name': name_from_tuple,
                                                   'author_orcid': author.get_authorORCID()})
            else:
                empiar_author = self._model_get_or_create_auto_id(EmpiarAuthorOrcid,
                                                                  {'name': author.get_valueOf_(),
                                                                   'author_orcid': author.get_authorORCID()})
                empiar_author.save()
                self._model_get_or_create_auto_id(EmpiarAuthor,
                                                  {'author': empiar_author, 'entry': entry, 'order_id': i})

            i += 1

    def _load_version_history(self, admin, entry):
        """
        Load version history from the database according to schema
        @param admin: admin part of the parsed XML tree from which authors are extracted
        @param entry: related entry from the EmpiarAnnotEntry table
        """

        versionHistory = admin.get_versionHistory()
        if versionHistory:
            versions = versionHistory.get_version()
            if versions:
                eVersionHistory = entry.empiarversionhistory_set.all()
                if eVersionHistory:
                    for version in eVersionHistory:
                        version.delete()
                for xmlVersion in versions:
                    if xmlVersion.get_annotator():
                        empiar_annotator = annotator = xmlVersion.get_annotator().get_valueOf_()
                        if ' ' in annotator:
                            first_name = annotator.rsplit(' ', 1)[0]
                            last_name = annotator.rsplit(' ', 1)[1]
                            try:
                                empiar_annotator = User.objects.filter(first_name=first_name, last_name=last_name,
                                                                       groups__name='annotators')[0]
                            except IndexError:
                                pass
                        else:
                            try:
                                empiar_annotator = User.objects.filter(username=annotator, groups__name='annotators')[0]
                            except IndexError:
                                pass

                        self._model_get_or_create_auto_id(EmpiarVersionHistory,
                                                          {
                                                              'entry': entry,
                                                              'annotator': empiar_annotator,
                                                              'version_number': xmlVersion.get_versionNumber(),
                                                              'date': xmlVersion.get_date(),
                                                              'status_code': xmlVersion.get_statusCode(),
                                                              'details': xmlVersion.get_details(),
                                                          })

    def _load_corresponding_author(self, admin, entry):
        """
        Load a corresponding author from XML file according to schema
        @param admin: admin part of the parsed XML tree from which information about
        a corresponding author is extracted
        @param entry: related entry from the EmpiarAnnotEntry table
        """

        corr_auth_dict = {
            'author_orcid': admin.get_correspondingAuthor().get_authorORCID(),
            'first_name': admin.get_correspondingAuthor().get_firstName(),
            'middle_name': admin.get_correspondingAuthor().get_middleName(),
            'last_name': admin.get_correspondingAuthor().get_lastName(),
            'organization': admin.get_correspondingAuthor().get_organization().get_valueOf_(),
            'street': admin.get_correspondingAuthor().get_street(),
            'town_or_city': admin.get_correspondingAuthor().get_townOrCity(),
            'state_or_province': admin.get_correspondingAuthor().get_stateOrProvince(),
            'country': admin.get_correspondingAuthor().get_country()
        }

        if self.prerel:
            corr_auth_dict['country'] = self._proper_country(corr_auth_dict['country'])

        corr_auth_dict['post_or_zip'] = admin.get_correspondingAuthor().get_postOrZipCode()
        if self.prerel:
            corr_auth_dict['email'] = admin.get_correspondingAuthor().get_email().get_valueOf_()
        else:
            corr_auth_dict['email'] = admin.get_correspondingAuthor().get_email().get_valueOf_().replace('@', ' [at] ')

        telephone_obj = admin.get_correspondingAuthor().get_telephone()
        if telephone_obj:
            corr_auth_dict['telephone'] = xstr(telephone_obj.get_local())
        else:
            corr_auth_dict['telephone'] = telephone_obj

        fax_obj = admin.get_correspondingAuthor().get_fax()
        if fax_obj:
            corr_auth_dict['fax'] = xstr(fax_obj.get_local())
        else:
            corr_auth_dict['fax'] = fax_obj

        try:
            corresponding_author = entry.corresponding_author
            for attr, value in iteritems(corr_auth_dict):
                setattr(corresponding_author.author, attr, value)
            corresponding_author.author.save()
        except (AttributeError, EmpiarCorrespondingAuthor.DoesNotExist, EmpiarPrInvCorrAuthor.DoesNotExist):
            if self.prerel:
                corr_auth_dict['entry'] = entry
                corresponding_author = self._model_get_or_create_auto_id(EmpiarPreRelCorrespondingAuthor,
                                                                         corr_auth_dict)
            else:
                ref = self._model_get_or_create_auto_id(EmpiarPrInvCorrAuthor, corr_auth_dict)
                corresponding_author = self._model_get_or_create_auto_id(EmpiarCorrespondingAuthor, {'author': ref})

        return corresponding_author

    def _load_principal_investigators(self, admin, entry):
        """
        Load principal investigators from XML file according to schema
        @param admin: admin part of the parsed XML tree from which information about
        a principal investigator is extracted
        @param entry: related entry from the EmpiarAnnotEntry table
        """

        if not self.prerel:
            existingPrincipalInvestigatorsList = entry.empiarprincipalinvestigator_set.all()
            if existingPrincipalInvestigatorsList:
                for existingPrincipalInvestigator in existingPrincipalInvestigatorsList:
                    author = existingPrincipalInvestigator.author

                    try:
                        existingPrincipalInvestigator.delete()
                    except IntegrityError as e:
                        self.log.write("Error deleting a principal investigator:")
                        logger.error(e)
                        logger.error(traceback.format_exc())

                    try:
                        author.delete()
                    except IntegrityError as e:
                        self.log.write("Error deleting an author:")
                        logger.error(e)
                        logger.error(traceback.format_exc())

        # Extract the information about principal investigators:
        principalInvestigators = admin.get_principalInvestigator()
        if principalInvestigators:
            for principalInvestigator in principalInvestigators:
                pr_inv_dict = {
                    'author_orcid': principalInvestigator.get_authorORCID(),
                    'first_name': principalInvestigator.get_firstName(),
                    'middle_name': principalInvestigator.get_middleName(),
                    'last_name': principalInvestigator.get_lastName(),
                    'organization': principalInvestigator.get_organization().get_valueOf_(),
                    'street': principalInvestigator.get_street(),
                    'town_or_city': principalInvestigator.get_townOrCity(),
                    'state_or_province': principalInvestigator.get_stateOrProvince(),
                    'country': principalInvestigator.get_country()
                }

                if self.prerel:
                    pr_inv_dict['country'] = self._proper_country(pr_inv_dict['country'])

                pr_inv_dict['post_or_zip'] = principalInvestigator.get_postOrZipCode()
                if self.prerel:
                    pr_inv_dict['email'] = principalInvestigator.get_email().get_valueOf_()
                else:
                    pr_inv_dict['email'] = principalInvestigator.get_email().get_valueOf_().replace('@', ' [at] ')
                telephone_obj = principalInvestigator.get_telephone()
                if telephone_obj:
                    pr_inv_dict['telephone'] = xstr(telephone_obj.get_local())
                else:
                    pr_inv_dict['telephone'] = telephone_obj

                fax_obj = principalInvestigator.get_fax()
                if fax_obj:
                    pr_inv_dict['fax'] = xstr(fax_obj.get_local())
                else:
                    pr_inv_dict['fax'] = fax_obj

                if self.prerel:
                    pr_inv_dict['entry'] = entry
                    self._model_get_or_create_auto_id(EmpiarPreRelPrincipalInvestigator, pr_inv_dict)
                else:
                    ref = self._model_get_or_create_auto_id(EmpiarPrInvCorrAuthor, pr_inv_dict)
                    self._model_get_or_create_auto_id(EmpiarPrincipalInvestigator,
                                                      {'author': ref, 'entry': entry})

    def _load_image_sets(self, header_xml, entry):
        """
        Load image sets from XML file according to schema
        @param header_xml: parsed XML tree from which information about
        image sets is extracted
        @param entry: related entry - foreign key to EmpiarEntry table
        """

        if self.prerel:
            imageSetClass = EmpiarPreRelImageSet
            segmentationClass = EmpiarPreRelSegmentation
        else:
            imageSetClass = EmpiarImageSet
            segmentationClass = EmpiarSegmentation

            existingImageSetsList = entry.empiarimageset_set.all()

            for existingImageSet in existingImageSetsList:
                existingSegmentationList = existingImageSet.empiarsegmentation_set.all()
                for existingSegmentation in existingSegmentationList:
                    existingSegmentation.delete()

                existingImageSet.delete()

        # Extract the information about image sets:
        imageSets = header_xml.get_imageSet()

        for imageSet in imageSets:
            frame_range_min = frame_range_max = None
            if imageSet.get_frameRange():
                frame_range_min = imageSet.get_frameRange().frameRangeMin
                frame_range_max = imageSet.get_frameRange().frameRangeMax

            category = imageSet.get_category()
            if self.prerel:
                for categorytype in category_choices:
                    if categorytype[1] == category:
                        if categorytype[0] == 'OT':
                            category = "(u'%s', u'%s')" % (categorytype[0], category)
                        else:
                            category = "(u'%s', u'')" % categorytype[0]

            header_format = imageSet.get_headerFormat()
            if self.prerel:
                for headertype in header_format_choices:
                    if headertype[1] == header_format:
                        if headertype[0] == 'OT':
                            header_format = "(u'%s', u'%s')" % (headertype[0], header_format)
                        else:
                            header_format = "(u'%s', u'')" % headertype[0]

            data_format = imageSet.get_dataFormat()
            if self.prerel:
                for datatype in data_format_choices:
                    if datatype[1] == data_format:
                        if datatype[0] == 'OT':
                            data_format = "(u'%s', u'%s')" % (datatype[0], data_format)
                        else:
                            data_format = "(u'%s', u'')" % datatype[0]

            voxel_type = imageSet.get_voxelType()
            if self.prerel:
                for voxeltypetype in voxel_type_choices:
                    if voxeltypetype[1] == voxel_type:
                        if voxeltypetype[0] == 'OT':
                            voxel_type = "(u'%s', u'%s')" % (voxeltypetype[0], voxel_type)
                        else:
                            voxel_type = "(u'%s', u'')" % voxeltypetype[0]

            image_width = imageSet.get_dimensions().get_imageWidth()
            image_height = imageSet.get_dimensions().get_imageHeight()
            pixel_width = imageSet.get_dimensions().get_pixelWidth()
            pixel_height = imageSet.get_dimensions().get_pixelHeight()

            eImageSet = self._model_get_or_create_auto_id(
                imageSetClass, {'entry': entry, 'name': imageSet.get_name(),
                                'directory': imageSet.get_directory(),
                                'category': category,
                                'header_format': header_format,
                                'data_format': data_format,
                                'num_images_or_tilt_series': imageSet.get_numImagesOrTiltSeries(),
                                'frames_per_image': imageSet.get_framesPerImage(),
                                'frame_range_min': frame_range_min,
                                'frame_range_max': frame_range_max,
                                'voxel_type': voxel_type,
                                'image_width': None if image_width == 'variable' else image_width,
                                'pixel_width': None if pixel_width == 'variable' else pixel_width,
                                'image_height': None if image_height == 'variable' else image_height,
                                'pixel_height': None if pixel_height == 'variable' else pixel_height,
                                'micrographs_file_pattern': imageSet.get_micrographsFilePattern(),
                                'picked_particles_file_pattern': imageSet.get_pickedParticlesFilePattern(),
                                'picked_particles_directory': imageSet.get_pickedParticlesDirectory(),
                                'details': imageSet.get_details()
                                })

            seg_list = imageSet.get_segmentationList()
            if seg_list:
                segs = seg_list.get_segmentation()
                for seg in segs:
                    self._model_get_or_create_auto_id(segmentationClass,
                                                      {'imageset': eImageSet, 'id': seg.get_segmentationId(),
                                                       'file': seg.get_file(),
                                                       'description': seg.get_description(),
                                                       'original_files': seg.get_originalFiles(),
                                                       'original_format': seg.get_originalFormat()})

    def _loadEntry(self, header_xml, e):
        """
        Load entry according to the version 1.8 schema
        @param header_xml: parsed XML tree of a header file
        @param e: EmpiarEntry which is modified by this method
        """
        admin = header_xml.get_admin()
        e.title = admin.get_title()
        e.status = admin.get_currentStatus()
        experiment = admin.get_experimentType()

        if self.prerel:
            e.deposition_date = datetime.date.today()
            self._load_corresponding_author(admin, e)
            try:
                for experiment_type in experiment_type_choices:
                    if experiment_type[1] == experiment:
                        e.experiment_type = experiment_type[0]
            except (SyntaxError, ValueError):
                e.experiment_type = 1
        else:
            e.experiment_type = experiment
            keyDates = admin.get_keyDates()
            e.deposition_date = keyDates.get_depositionDate()
            e.release_date = keyDates.get_releaseDate()
            e.update_date = keyDates.get_updateDate()
            e.obsolete_date = keyDates.get_obsoleteDate()

            e.corresponding_author = self._load_corresponding_author(admin, e)
            e.dataset_size = admin.get_datasetSize().get_valueOf_() + ' ' + \
                             admin.get_datasetSize().get_units()

            # Number of files per entry is not being supported by XML. To avoid creation of redundant files, it has been
            #  decided to populate this field directly from the annotation table
            e.num_of_files = EmpiarAnnotEntry.objects.get(empiar_id=e.id).num_of_files

        # Save entry
        e.save()

        # Load tables that have the entry id as a foreign key
        self._load_authors(admin, e)
        self.log.write("\nSuccessfully loaded entry authors")
        self._load_cross_references(header_xml, e)
        self.log.write("\nSuccessfully loaded cross references")
        self._load_principal_investigators(admin, e)
        self.log.write("\nSuccessfully loaded principal investigators")
        self._load_image_sets(header_xml, e)
        self.log.write("\nSuccessfully loaded image sets")

        scipion_workflow = admin.get_scipionWorkflow()
        if scipion_workflow:
            if self.prerel:
                scipionClass = EmpiarPreRelScipion
            else:
                scipionClass = EmpiarScipion
                if hasattr(e, 'empiarannotscipion'):
                    existingScipion = e.empiarscipion
                    if existingScipion:
                        existingScipion.delete()
            self._model_get_or_create_auto_id(scipionClass, {'entry': e, 'scipion_workflow': scipion_workflow})
            self.log.write("\nSuccessfully loaded Scipion")

        if not self.prerel:
            self._load_version_history(admin, e)
            self.log.write("\nSuccessfully loaded version history")

    def loadEntry(self, headerFile):
        """
        Load a single entry into database.
        @param headerFile: Header XML file following current EMPIAR XML schema. If this is prerel loading, then
        headerFile is a string that contains the XML
        """

        if self.prerel:
            header_xml = empiar.parseString(headerFile, silence=True)

            try:
                empiar_prerel_entry = EmpiarPreRelEntry.objects.latest('id')
                empiar_annot_entry = EmpiarAnnotEntry.objects.latest('id')
                empiar_prerel_rights_entry = EmpiarPreRelEntryRights.objects.latest('entry_tmp_id')

                id = max(empiar_prerel_entry.id, empiar_annot_entry.id, empiar_prerel_rights_entry.entry_tmp_id) + 1
            except EmpiarAnnotEntry.DoesNotExist:
                id = max(empiar_prerel_entry.id, empiar_prerel_rights_entry.entry_tmp_id) + 1
            except (EmpiarPreRelEntry.DoesNotExist, EmpiarPreRelEntryRights.DoesNotExist):
                id = 1
            e, success = model_get_or_create_from_id(EmpiarPreRelEntry, {'id': id, 'locked': False})

            e.save()

        else:
            try:
                entry_info = read_file(headerFile)
                header_xml = empiar.parseString(entry_info, silence=True)
            except Exception as e:
                logger.error(e)
                logger.error(traceback.format_exc())
                self.log.write("\nUnable to parse file: %s" % headerFile)
                self.log.write("\nError: %s" % xstr(e))
                return

            id = header_xml.accessionCode
            try:
                e = EmpiarEntry.objects.get(id=id)
            except EmpiarEntry.DoesNotExist:
                e = EmpiarEntry(id=id)

        e.schema_version = header_xml.get_schemaVersion()
        try:
            self._loadEntry(header_xml, e)
            self.log.write("\nSuccessfully updated the database from the entry XML header file")
            return e.id
        except Exception as e:
            logger.error(e)
            logger.error(traceback.format_exc())
            self.log.write("\nUnable to load from file: %s" % headerFile)
            self.log.write("\nError: %s" % xstr(e))
            return

    def loadEntries(self, upload_dir, archive_dir, entry_id):
        """
        Load XML header info into database
        """

        self.log.write("***EMPIAR Loader***")

        dt1 = datetime.datetime.now()
        if self.prerel:
            headerFile = entry_id
        else:
            headerFile = os.path.join(upload_dir, '.private', entry_id + '.xml')
            if not file_exists(headerFile):
                headerFile = os.path.join(archive_dir, entry_id, entry_id + '.xml')
            self.log.write("\nLoading entry from: %s" % headerFile)
        result = self.loadEntry(headerFile)
        dt2 = datetime.datetime.now()

        self.log.write("\nTotal time: %s" % (dt2 - dt1))

        return result


def headers_loader(upload_dir, archive_dir, entry_id, prerel=False):
    """
    Load the entry XML header into the EMPIAR database
    """

    # Create EMPIARLoader object
    empiarLoader = EMPIARLoader()

    if prerel:
        empiarLoader.prerel = True

    result = empiarLoader.loadEntries(upload_dir, archive_dir, entry_id)

    log = empiarLoader.log.getvalue()

    if prerel:
        if 'Successfully updated the database' in log and result and (
                isinstance(result, int) or (isinstance(result, BaseStrType_) and result.isdigit())):
            return result
        else:
            return -1
    return log

