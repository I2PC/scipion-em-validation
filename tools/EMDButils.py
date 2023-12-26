import os
import requests
import urllib.request
import gzip
from bs4 import BeautifulSoup

def does_map_exist(emdbid):
    """
    Checks if a map exist in EMDB given an id
    """
    url_rest_api = 'https://www.ebi.ac.uk/emdb/api/entry/%s' % emdbid
    print("Checking if EMD-%s exists..." % emdbid)
    response = requests.get(url_rest_api)
    if response.status_code == 200:
        return True
    else:
        return False

def get_map_metadata(emdbid):
    """
    Returns the sampling, threshold and resolution from an EMDB map
    """
    url_rest_api = 'https://www.ebi.ac.uk/emdb/api/entry/%s' % emdbid
    print("Getting some EMD-%s metadata parameters (sampling, threshold and resolution) ..." % emdbid)
    try:
        json_results = requests.get(url_rest_api).json()
        sampling = float(json_results["map"]["pixel_spacing"]["x"]["valueOf_"])
        threshold = float(json_results["map"]["contour_list"]["contour"][0]["level"])
        resolution = float(json_results["structure_determination_list"]["structure_determination"][0]["image_processing"][0]["final_reconstruction"]["resolution"]["valueOf_"])
        return sampling, threshold, resolution
    except:
        return None, None, None

def has_halfmaps(emdbid):
    """
    Checks if an EMDB map has half-maps associated
    """
    url_rest_api = 'https://www.ebi.ac.uk/emdb/api/entry/supplement/%s' % emdbid
    print("Checking if EMD-%s has associated halfmaps in %s..." % (emdbid, url_rest_api))
    json_results = requests.get(url_rest_api).json()
    try:
        json_results['interpretation']['half_map_list']['half_map'][0]['file']
        json_results['interpretation']['half_map_list']['half_map'][1]['file']
    except KeyError:
        return False
    return True

def download_emdb_halfmaps(emdbid, directory):
    """
    Downloads half-maps associated to an EMDB map
    """
    url_supplement_info_rest_api = 'https://www.ebi.ac.uk/emdb/api/entry/supplement/%s' % emdbid
    json_results = requests.get(url_supplement_info_rest_api).json()
    half_maps = [json_results['interpretation']['half_map_list']['half_map'][0]['file'], json_results['interpretation']['half_map_list']['half_map'][1]['file']]

    url_ftp_other = 'ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-%s/other/%s'

    print("Downloading associated halfmaps...")

    for map_gz_name in half_maps:
        map_url = url_ftp_other % (emdbid, map_gz_name)
        map_name = map_gz_name.replace('.gz', '')
        noCompressName = os.path.join(directory, map_name)
        compressName = os.path.join(directory, map_gz_name)
        print('Downloading file from %s' % map_url)
        urllib.request.urlretrieve(map_url, filename=compressName)
        gunzip(compressName, noCompressName)

    return half_maps

def has_atomicmodel(emdbid):
    """
    Returns the associated atomic model of an EMDB map
    """
    url_rest_api = 'https://www.ebi.ac.uk/emdb/api/entry/%s' % emdbid
    print("Checking if EMD-%s has an atomic model associated..." % emdbid)
    try:
        json_results = requests.get(url_rest_api).json()
        return json_results['crossreferences']['pdb_list']['pdb_reference'][0]['pdb_id']
    except:
        return None

def download_atomicmodel(pdbid, path):
    """
    Downloads an atomic model
    """
    url = 'https://files.rcsb.org/download/%s.cif' % pdbid
    print("Downloading atomic model %s..." % pdbid)
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        cifpath = os.path.join(path, pdbid + '.cif')
        with open(cifpath, 'w') as fd:
            fd.write(response.text)
        return cifpath
    else:
        return None

def proper_map_axis_order(emdbid):
    """
    Checks if axis volume are ordered as x, y and z
    """
    url_rest_api = 'https://www.ebi.ac.uk/emdb/api/entry/%s' % emdbid
    print("Checking if EMD-%s axis map is ordered as x, y and z ..." % emdbid)
    try:
        json_results = requests.get(url_rest_api).json()
        fast = json_results['map']['axis_order']['fast'].lower()
        medium = json_results['map']['axis_order']['medium'].lower()
        slow = json_results['map']['axis_order']['slow'].lower()
        if fast == 'x' and medium == 'y' and slow == 'z':
            print("... yes it is.")
            return True
        else:
            print("... no it's not.")
            return False
    except:
        return False
    return True

def gunzip(gzpath, path):
    gzf = gzip.open(gzpath)
    f = open(path, 'wb')
    f.write(gzf.read())
    f.close()
    gzf.close()
    os.remove(gzpath)

def get_emdb_entries(output_file):
    """
    Gets all EMDB entries
    """
    url = 'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/'
    response = requests.get(url)

    if response.status_code == 200:
        content = BeautifulSoup(response.content, 'html.parser')
        emdb_entries = content.find_all('a')

        with open(output_file, 'a') as output:
            output.writelines(('%s\n' % emdb_entry.get_text().replace('/','') for emdb_entry in emdb_entries if 'EMD-' in emdb_entry.get_text()))
    else:
        print('Error getting the page:', response.status_code)

def is_spa_entry(emdbid):
    """
    Checks if an EMDB entry was obtained using Single Particle Analysis method
    """
    url_rest_api = 'https://www.ebi.ac.uk/emdb/api/entry/%s' % emdbid
    print("Checking if EMD-%s was obtained using Single Particle Analysis method..." % (emdbid))
    json_results = requests.get(url_rest_api).json()
    try:
        method = json_results["structure_determination_list"]["structure_determination"][0]["method"]
        if method == 'singleParticle':
            return True
        else:
            return False
    except:
        return False

def check_entry_level(emdb_entries_input_file, level1_output_file, levelA_output_file, failures_output_file):
    """
    Checks which EMDB entry has half-maps (level 1) and atomic model associated (level A)
    """
    with open(emdb_entries_input_file, 'r') as input:
        emdb_entries = [emdb_entry.strip() for emdb_entry in input.readlines()]

    level1_entries = []
    levelA_entries = []
    failures = []

    for emdb_entry in emdb_entries:
        try:
            if has_halfmaps(emdb_entry):
                level1_entries.append(emdb_entry)
            if has_atomicmodel(emdb_entry):
                levelA_entries.append(emdb_entry)
        except:
            failures.append(emdb_entry)

    with open(level1_output_file, 'a') as output:
        output.writelines('%s\n' % level1_entry for level1_entry in level1_entries)
    with open(levelA_output_file, 'a') as output:
        output.writelines('%s\n' % levelA_entry for levelA_entry in levelA_entries)
    with open(failures_output_file, 'a') as output:
        output.writelines('%s\n' % failure for failure in failures)

def get_subsets(levelO_input_file, level1_input_file, levelA_input_file, level0notAnot1_output_file, level0Anot1_output_file, level01AnotA_output_file, level01A_output_file):
    """
    Creates several lists:
        - Level 0 - Level 1 - Level A: Entries with map but without half-maps and atomic model
        - (Level 0 & Level A) - Level 1: Entries with map and atomic model but without half-maps
        - (Level 0 & Level 1) - Level A:: Entries with map and half-maps but without atomic model
        - Level 0 & Level 1 & Level A: Entries with maps, half-maps and atomic model
    """
    with open(levelO_input_file, 'r') as input:
        level0_entries = [emdb_entry.strip() for emdb_entry in input.readlines()]
    with open(level1_input_file, 'r') as input:
        level1_entries = [emdb_entry.strip() for emdb_entry in input.readlines()]
    with open(levelA_input_file, 'r') as input:
        levelA_entries = [emdb_entry.strip() for emdb_entry in input.readlines()]

    level0notAnot1_entries = list(set(level0_entries) - set(levelA_entries) - set(level1_entries))
    level0Anot1_entries = list((set(level0_entries) & set(levelA_entries)) - set(level1_entries))
    level01AnotA_entries = list((set(level0_entries) & set(level1_entries)) - set(levelA_entries))
    level01A_entries = list(set(level0_entries) & set(levelA_entries) & set(level1_entries))

    with open(level0notAnot1_output_file, 'a') as output:
        output.writelines('%s\n' % entry for entry in level0notAnot1_entries)
    with open(level0Anot1_output_file, 'a') as output:
        output.writelines('%s\n' % entry for entry in level0Anot1_entries)
    with open(level01AnotA_output_file, 'a') as output:
        output.writelines('%s\n' % entry for entry in level01AnotA_entries)
    with open(level01A_output_file, 'a') as output:
        output.writelines('%s\n' % entry for entry in level01A_entries)
