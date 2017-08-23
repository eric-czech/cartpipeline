"""
Wrapper API for CGDS (cbioportal) Web Services

See http://www.cbioportal.org/web_api.jsp for more details.
"""
import pandas as pd
import os
import logging
import time
import hashlib

logger = logging.getLogger(__name__)

BASE_URL = 'http://www.cbioportal.org/public-portal/webservice.do?'


def _to_url(cmd, data=None):
    url = '{}cmd={}'.format(BASE_URL, cmd)
    if data:
        url += '&' + '&'.join(['{}={}'.format(k, v) for k, v in data.items()])
    return url


def _to_batches(sequence, batch_size):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(sequence), batch_size):
        yield sequence[i:i+batch_size]


def _get(cmd, data=None):
    url = _to_url(cmd, data)
    logger.debug('Invoking cBioPortal endpoint: {}'.format(url))
    return pd.read_csv(url, sep='\t', comment='#')


def _is_id(idv):
    return isinstance(idv, str) and len(idv) > 0


def _is_iterable(seq, check_empty=True):
    # Ensure sequence is not a string, which is technically iterable
    # but not desired in this context
    if isinstance(seq, str):
        return False
    try:
        _ = (e for e in seq)
        if check_empty and len(seq) == 0:
            return False
        return True
    except TypeError:
        return False


def get_cancer_studies():
    return _get('getCancerStudies')


def get_cancer_types():
    return _get('getTypesOfCancer')


def get_genetic_profiles(cancer_study_id):
    if not _is_id(cancer_study_id):
        raise ValueError('Cancer Study ID must be a non-empty string (e.g. "cellline_ccle_broad")')
    return _get('getGeneticProfiles', {'cancer_study_id': cancer_study_id})


def get_case_lists(cancer_study_id):
    if not _is_id(cancer_study_id):
        raise ValueError('Cancer Study ID must be a non-empty string (e.g. "cellline_ccle_broad")')
    return _get('getCaseLists', {'cancer_study_id': cancer_study_id})


def get_clinical_data(case_list_id):
    if not _is_id(case_list_id):
        raise ValueError('Case list ID must be a non-empty string (e.g. "cellline_ccle_broad_all")')
    return _get('getClinicalData', {'case_set_id': case_list_id})


def _get_batch_batch_result(gene_ids, batch_size, cmd, args, print_progress=True, max_attempts=5, failure_pause_secs=30):
    res = None
    gene_id_batches = list(_to_batches(gene_ids, batch_size))
    n = len(gene_id_batches)
    m = max(int(n / 10.), 1)
    for i, gene_ids in enumerate(gene_id_batches):
        if print_progress and i % m == 0:
            logger.info('Processing batch {} of {}'.format(i + 1, n))
        data = dict(args)
        data['gene_list'] = ','.join(gene_ids)

        attempts = 0
        while True:
            attempts += 1
            try:
                part = _get(cmd, data)
                break
            except:
                if attempts > max_attempts:
                    raise
                logger.warn('An http error occurred.  Will try again in {} seconds ...'.format(failure_pause_secs))
                time.sleep(failure_pause_secs)

        if res is None:
            res = part
        else:
            res = res.append(part)
    return res


def get_genetic_profile_data(case_list_id, genetic_profile_id, gene_ids,
                             batch_size=50, print_progress=True, cache_dir=None):
    if not _is_id(case_list_id):
        raise ValueError('Case list ID must be a non-empty string (e.g. "cellline_ccle_broad_all")')
    if not _is_id(genetic_profile_id):
        raise ValueError('Genetic Profile ID must be a non-empty string (e.g. "cellline_ccle_broad_log2CNA")')
    if not _is_iterable(gene_ids):
        raise ValueError('Gene IDs must be iterable and non-empty')

    if cache_dir is not None:
        cache_key = ':'.join([case_list_id, genetic_profile_id] + list(sorted(set(gene_ids))))
        cache_key = hashlib.md5(cache_key.encode('utf-8')).hexdigest()
        cache_path = os.path.join(cache_dir, cache_key + '.pkl')
        if os.path.exists(cache_path):
            logging.info('Using batch profile data result from cache path "{}"'.format(cache_path))
            return pd.read_pickle(cache_path)

    # Split given gene list into batches and combine results from each batch
    args = {
        'id_type': 'gene_symbol',
        'case_set_id': case_list_id,
        'genetic_profile_id': genetic_profile_id
    }
    res = _get_batch_batch_result(gene_ids, batch_size, 'getProfileData', args, print_progress=print_progress)
    if cache_dir is None:
        return res

    if not os.path.exists(cache_dir):
        os.mkdir(cache_dir)
    cache_key = ':'.join([case_list_id, genetic_profile_id] + list(sorted(set(gene_ids))))
    cache_key = hashlib.md5(cache_key.encode('utf-8')).hexdigest()
    cache_path = os.path.join(cache_dir, cache_key + '.pkl')
    logger.info('Storing CGDS result for query [case_list = {}, profile_id = {}] in cache at "{}"'
                .format(case_list_id, genetic_profile_id, cache_path))
    res.to_pickle(cache_path)
    return res
