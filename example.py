from scidata_xwalker.scidata_xwalker import \
    flatten_json_iterative_solution,\
    datasetmodder,\
    remove_extra_metadata,\
    cleanup_flattened,\
    get_semantics,\
    crosswalker, \
    group_link_override, \
    bin_grouper, \
    binner
from scidata import SciData
import json

"""Merged_dictionary is the discrete dataset in JSON format and should have
a single root key. Central_table is the name of the root key. The
crosswalks, namespaces, and ontterms variables were inspected from MySQL
tables. In practice, these tables are accessed directly (for example using
Django models). The crosswalks, namespaces, and ontterms variables were
truncated to only included those entries relevant to this sample dataset
defined in merged_dictionary. group_overrides (used by group_link_override)
use regex to modify grouping within the document sci_links (used by
scilinker) use regex to create internal links within document. Generate
output without running scicleanup to view # key which is the target of the
group_link_override and scilinker functions """

merged_dictionary = \
    {"compounds": {
        "id": 596,
        "dsstox_compound_id": "DTXCID10596",
        "smiles": "CCC1=CC=CC=C1",
        "jchem_inchi_key": "YNQLUTRBYVCPMQ-UHFFFAOYSA-N",
        "acd_iupac_name": "Ethylbenzene",
        "mol_formula": "C8H10",
        "mol_weight": 106.168,
        "qsar_predicted_properties": [{
            "id": 5909,
            "efk_dsstox_compound_id": 596,
            "efk_qsar_model_id": 14,
            "result_value": 130.52,
            "report_filename": "DTXCID10596-TEST_Thermal_Conductivity.html",
            "models": {
                "id": 14,
                "efk_chemprop_endpoint_id": 15,
                "name": "TEST_Thermal_Conductivity"
            }}, {
            "id": 5982,
            "efk_dsstox_compound_id": 596,
            "efk_qsar_model_id": 4,
            "result_value": 0.89,
            "report_filename": "DTXCID10596-TEST_Density.html",
            "models": {
                "id": 4,
                "efk_chemprop_endpoint_id": 5,
                "name": "TEST_Density",
            }}]}}

object_data_json = json.loads(json.dumps(merged_dictionary))

central_table = 'compounds'

crosswalks = [{'id': 16,
               'table': 'compounds',
               'field': 'dsstox_compound_id',
               'ontterm_id': 323,
               'sdsection': 'system',
               'sdsubsection': 'compound',
               'sdsubsubsection': 'identifier',
               'category': None,
               'unit': None,
               'datatype': 'string'},
              {'id': 23,
               'table': 'compounds',
               'field': 'smiles',
               'ontterm_id': 80,
               'sdsection': 'system',
               'sdsubsection': 'compound',
               'sdsubsubsection': 'identifier',
               'category': None,
               'unit': None,
               'datatype': 'string'},
              {'id': 47,
               'table': 'qsar_predicted_properties',
               'field': 'result_value',
               'ontterm_id': 9999,
               'sdsection': 'dataset',
               'sdsubsection': 'exptdata',
               'sdsubsubsection': None,
               'category': None,
               'unit': None,
               'datatype': 'string'},
              {'id': 48,
               'table': 'models',
               'field': 'name',
               'ontterm_id': 9999,
               'sdsection': 'methodology',
               'sdsubsection': 'model',
               'sdsubsubsection': '',
               'category': None,
               'unit': None,
               'datatype': 'string'},
              {'id': 49,
               'table': 'qsar_predicted_properties',
               'field': 'report_filename',
               'ontterm_id': 9999,
               'sdsection': 'dataset',
               'sdsubsection': 'suppdata',
               'sdsubsubsection': None,
               'category': None,
               'unit': None,
               'datatype': 'string'}]

onttermslist = [{'id': 323,
                 'title': 'Identifier',
                 'definition': 'One or more characters used to identify, '
                               'name, or characterize the nature, '
                               'properties, or contents of a thing. ['
                               'def-source: NCI]',
                 'code': 'NCIT_C25364',
                 'url': 'obo:NCIT_C25364',
                 'nspace_id': 2,
                 'sdsection': 'metadata',
                 'sdsubsection': 'identifier'},
                {'id': 80,
                 'title': 'SMILES string',
                 'definition': 'SMILES string corresponding to drug '
                               'structure [database_cross_reference: '
                               'PMID:14755292]',
                 'code': 'MI:2039',
                 'url': 'obo:MI_2039',
                 'nspace_id': 2,
                 'sdsection': 'system',
                 'sdsubsection': 'compound'},
                {'id': 9999,
                 'title': 'TEST',
                 'definition': 'TEST',
                 'code': 'test_9999',
                 'url': 'test:test_9999',
                 'nspace_id': 19,
                 'sdsection': 'system',
                 'sdsubsection': ''}]

nspaceslist = [
    {'id': 2, 'name': 'OBO Foundry', 'ns': 'obo',
     'path': 'http://purl.obolibrary.org/obo/',
     'homepage': 'http://www.ontobee.org'},
    {'id': 19, 'name': 'SciData', 'ns': 'sdo',
     'path': 'https://stuchalk.github.io/scidata/ontology/scidata.owl#',
     'homepage': 'https://stuchalk.github.io/scidata/'}]


"""Define group overrides to modify organization Generate SciData JSON-LD
before running scicleanup to see # key and value to assist in writing
group_overrides term 1 = < regex pattern to find in value # key > term 2 = <
regex pattern to find in value of # for match > { term 1 : term 2 } If term
1 depends on enumeration use parentheses to create regex groups. Group 2
should be the enumerated value ie. (\\d{1,}) If term 2 match depends on the
enumeration, use '$!@%' in the position of enumeration
 ## indicates group_link before override"""
group_overrides = {}
group_overrides.update({
    '(compounds;qsar_predicted_properties;)(\\d{1,})(\\/exptdata)':
    'compounds;qsar_predicted_properties;$!@%/data'})
group_overrides.update({
    '(compounds;qsar_predicted_properties;)(\\d{1,})(\\/suppdata)':
    'compounds;qsar_predicted_properties;$!@%/data'})

"""Define sci_links to create internal links between sections 
Generate
SciData JSON-LD before running scicleanup to see # key and value to assist
in writing sci_links 
term 1 = regex pattern to find in value # key 
term 2 = name of key to be added to identify relationship 
term 3 = regex pattern to find in value of # for match 
{ term 1 : { term 2 : term 3 } }
If term 1 depends on enumeration use parentheses to create regex groups.
Group 2 should be the enumerated value ie. (\\d{1,}) If term 3 match depends
on the enumeration, use '$!@%' in the position of enumeration 
Multiple key/value pairs can be included in terms 2 and 3 if needed """
sci_links = {}
sci_links.update({'(compounds;qsar_predicted_properties;)(\\d{1,})(\\/data)': {
    'model': 'compounds;qsar_predicted_properties;$!@%/model',
    'compound': 'compounds/compound'
}})


"""Define sci_groups"""
sci_groups = {}
sci_groups.update({"crystal": "crystal/$!@%/"})


uid = '596'

""""""
""""""

"""Initialize SciData object"""
test = SciData(uid)

crosswalker(central_table, object_data_json, crosswalks, 'field', False)

flat_object_data_json = flatten_json_iterative_solution(object_data_json)

object_data_json_flat_filtered = cleanup_flattened(flat_object_data_json)

group_link_override(object_data_json_flat_filtered, group_overrides)

nspaces = get_semantics(
    object_data_json_flat_filtered,
    onttermslist,
    nspaceslist)

bins = binner(object_data_json_flat_filtered)

bins_grouped = bin_grouper(bins)

datasetmod = datasetmodder(bins_grouped['dataset'])

remove_extra_metadata(bins_grouped['methodology'])

remove_extra_metadata(bins_grouped['system'])

test.aspects(bins_grouped['methodology'])
test.facets(bins_grouped['system'])
test.datapoint(datasetmod)

sourcecode = 'cross'
datasetname = 'walk'
unique_id = 'er'

test.namespaces(nspaces)
test.discipline('w3i:Chemistry')
test.subdiscipline('w3i:ChemicalInformatics')
test.docid('1')
test.version('1')
test.title('SciData Crosswalker')
test.author(['Dylan Johnson'])
test.description('Example of Crosswalker Functions')
test.publisher('Chalk Lab')
test.base(
    "https://scidata.unf.edu/" +
    sourcecode +
    ":" +
    datasetname +
    ":" +
    unique_id +
    "/")
test.permalink(
    "https://scidata.unf.edu/" +
    sourcecode +
    ":" +
    datasetname +
    ":" +
    unique_id +
    "/")
test.graph_uid(sourcecode + ":" + datasetname + ":" + unique_id)

test.scilinker(sci_links)

test.datagroup(sci_groups)

test.scicleanup()

print(json.dumps(test.output, indent=4, ensure_ascii=False))
