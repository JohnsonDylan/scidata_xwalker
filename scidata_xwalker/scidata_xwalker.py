from itertools import chain, starmap
import ast
import re

dev = False


def crosswalker(sci_dir, sci_input, crosswalks, cw_field, cw_table):
    """Match directory/table name and field to a crosswalk entry.
    If match exists, modify value to a string representation of a
    dictionary containing the value and crosswalk information """
    for x, y in sci_input.items():
        if type(y) in [str, int, float]:
            if cw_table:
                match = next(
                    (cw for cw in crosswalks
                     if cw[cw_field] == x and cw[cw_table] == sci_dir),
                    None)
            else:
                match = next(
                    (cw for cw in crosswalks if cw[cw_field] == x), None)
            if match:
                out = {"scidata_value": y}
                out.update(match)
                out.pop('updated', None)
                sci_input[x] = str(out)
        if isinstance(y, dict):
            crosswalker(x, y, crosswalks, cw_field, cw_table)
        if isinstance(y, list):
            for z in y:
                crosswalker(x, z, crosswalks, cw_field, cw_table)


def flatten_json_iterative_solution(dictionary):
    """Flatten a nested json file"""

    def unpack(parent_key, parent_value):
        """Unpack one level of nesting in json file"""
        # Unpack one level only!!!

        if isinstance(parent_value, dict):
            for key, value in parent_value.items():
                temp1 = parent_key + ';' + key
                yield temp1, value
        elif isinstance(parent_value, list):
            i = 0
            for value in parent_value:
                temp2 = parent_key + ';' + str(i)
                i += 1
                yield temp2, value
        else:
            yield parent_key, parent_value

            # Keep iterating until the termination condition is satisfied

    while True:
        # Keep unpacking the json file until all values are atomic elements
        # (not dictionary or list)
        dictionary = dict(
            chain.from_iterable(
                starmap(
                    unpack,
                    dictionary.items())))
        # Terminate condition: not any value in the json file is dictionary or
        # list
        if not any(
            isinstance(
                value,
                dict) for value in dictionary.values()) and not any(
            isinstance(
                value,
                list) for value in dictionary.values()):
            break

    return dictionary


def cleanup_flattened(sci_input):
    """Converts stringified dictionary values to literal dicts.
    Also populates scidata_dir, scidata_key, scidata_group"""
    output = {}
    for k, v in sci_input.items():
        if isinstance(v, str) and v is not None and v.startswith('{'):
            if re.search(r'(\d+)(?!.*\d)', k):
                reg = re.split(r'(\d+)(?!.*\d)', k)
                scidata_key_path = reg.pop().split(';', 1)[1]
                sci_group = ''.join([str(elem) for elem in reg])
            else:
                sci_group = k.rsplit(';', 1)[0]
                scidata_key_path = k.rsplit(';', 1)[1]
            valx = ast.literal_eval(v)
            val = ({'scidata_dir': k,
                    'scidata_key': k.rsplit(';', 1)[1],
                    'scidata_key_path': scidata_key_path,
                    'scidata_group': sci_group})
            if valx['sdsection'] == 'dataset':
                val.update(
                    {'scidata_group_link': k.rsplit(';', 1)[0] + '/'})
                val.update(
                    {'scidata_group_link_original': k.rsplit(';', 1)[0] + '/'})
            else:
                val.update({'scidata_group_link': k.rsplit(';', 1)
                           [0] + '/' + str(valx['sdsubsection'])})
                val.update({'scidata_group_link_original': k.rsplit(
                    ';', 1)[0] + '/' + str(valx['sdsubsection'])})
            val.update(valx)
            output.update({k: val})
    return output


def group_link_override(sci_input, group_overrides):
    """redefines the scidata_group_link value based on the rules defined in
    group_overrides """
    for k, v in sci_input.items():
        for pattern, replacement in group_overrides.items():
            if re.search(pattern, v['scidata_group_link']):
                try:
                    replace = replacement.replace(
                        '$!@%', re.match(pattern,
                                         v['scidata_group_link']).group(2))
                    v['scidata_group_link'] = replace
                except Exception:
                    v['scidata_group_link'] = replacement
            if re.search(pattern, v['scidata_group_link']):
                try:
                    replace = replacement.replace(
                        '$!@%', re.match(pattern,
                                         v['scidata_group_link']).group(2))
                    v['scidata_group_link'] = replace
                except Exception:
                    v['scidata_group_link'] = replacement


def get_semantics(sci_input, onttermslist, nspaceslist):
    """Scans object_data_json_flat_filtered to create list of all ontterms
    present """
    ontoterms = []
    namespaces = {}
    for x in sci_input.values():
        ontoterms.append(x['ontterm_id'])
        ontsearch = next(
            item for item in onttermslist if item["id"] == x['ontterm_id'])
        nspacesearch = next(
            item for item in nspaceslist
            if item["id"] == ontsearch['nspace_id'])
        namespaces.update({nspacesearch['ns']: nspacesearch['path']})
    return namespaces


def binner(sci_input):
    """Separates into the three SciData Bins System/Methodology/Dataset"""
    bins = {'methodology': [], 'system': [], 'dataset': []}
    for k, v in sci_input.items():
        bins[v['sdsection']].append(v)
    return bins


def bin_grouper(sci_input):
    """Groups bins based on the 'scidata_group_link' key"""
    bins_grouped = {}
    for a, b in sci_input.items():
        bin_groups_set = set()
        for c in b:
            bin_groups_set.add(c['scidata_group_link'])
        subset_list = []
        for d in bin_groups_set:
            sdsubsection_group = {}
            for c in b:
                if c['scidata_group_link'] == d:
                    if c['sdsection'] == 'dataset':
                        if sdsubsection_group.get(c['scidata_key']):
                            if isinstance(
                                    sdsubsection_group.get(
                                        c['scidata_key']),
                                    dict):
                                sdsubsection_group[c['scidata_key']] = [
                                    sdsubsection_group.get(c['scidata_key'])]
                            sdsubsection_group[c['scidata_key']].append(c)
                        else:
                            sdsubsection_group.update(
                                {"@id": 'datum', c['scidata_key']: c})
                    else:
                        if sdsubsection_group.get(c['scidata_key']):
                            if isinstance(
                                    sdsubsection_group.get(
                                        c['scidata_key']),
                                    dict):
                                sdsubsection_group[c['scidata_key']] = [
                                    sdsubsection_group.get(c['scidata_key'])]
                            sdsubsection_group[c['scidata_key']].append(c)
                        else:
                            sdsubsection_group.update(
                                {"@id": c['sdsubsection'],
                                 c['scidata_key']: c})
            subset_list.append(sdsubsection_group)
        if subset_list:
            bins_grouped.update({a: subset_list})
    return bins_grouped


def remove_extra_metadata(sci_input):
    """The dictionary containing the value and metadata is replaced with
    only the value (scidata_value) """
    for entry in sci_input:
        reference = {}
        for k, v in entry.items():
            if isinstance(v, dict):
                entry[k] = v['scidata_value']
                reference.update({'#': v['scidata_group_link']})
                reference.update({'##': v['scidata_group_link_original']})

                if dev:
                    reference.update({k + '#': v['scidata_group_link']})
                    reference.update({k + '##': v['scidata_dir']})
            if isinstance(v, list):
                entry[k] = []
                for vlist in v:
                    entry[k].append(vlist['scidata_value'])
                    reference.update({'#': vlist['scidata_group_link']})
                    reference.update(
                        {'##': vlist['scidata_group_link_original']})
                    if dev:
                        reference.update(
                            {k + '#': vlist['scidata_group_link']})
                        reference.update({k + '##': vlist['scidata_dir']})
                entry[k] = set(entry[k])
                entry[k] = list(entry[k])
                if len(entry[k]) == 1:
                    entry[k] = entry[k][0]
        entry.update(reference)


def datasetmodder(sci_input):
    """Functions for distributing data amongst datum and value sections"""
    datasetmod = []
    for x in sci_input:
        datumtypes = set()
        for k, v in x.items():
            if k != '@id':
                datumtypes.add(v['sdsubsection'])
        datumset = []
        for sdsub in datumtypes:
            value = {}
            for k, v in x.items():
                if k != '@id':
                    if v['sdsubsection'] == sdsub:
                        value.update({"@id": "value",
                                      "@type": "sdo:value",
                                      v['scidata_key']: v['scidata_value'],
                                      # '#': v['scidata_group_link'].split('/')[0],
                                      # '##': v['scidata_group_link_original'].split('/')[0]
                                      })
            datum = ({
                "@id": "datum",
                "@type": "sdo:" + sdsub,
                "value": value})
            datumset.append(datum)

        dataset = {
            "@id": "datapoint",
            "@type": "sdo:datapoint",
            "data": datumset,
            '#': v['scidata_group_link'].split('/')[0],
            '##': v['scidata_group_link_original'].split('/')[0]
        }
        datasetmod.append(dataset)
    return datasetmod
