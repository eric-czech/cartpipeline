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


# #### Searching CGDS Metadata Examples ####
# pd.set_option('display.max_colwidth', 500)
# dt = cgds.get_genetic_profiles('cellline_ccle_broad')
# dt
# ####

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






# def melt_raw_data(d, rm_null_values=True, value_numeric=True):
#     """
#     Transforms raw CGDS data from wide format (cell lines in columns) to long format.
#
#     Note that the uniqueness of gene + cell line id combinations will also be ensured.
#
#     :param d: Raw CGDS DataFrame
#     :param rm_null_values: If true, rows with a VALUE of null will be removed after transformation to long format
#     :param value_numeric: If true, VALUE will be assumed to be numeric and will be cast as such (and results
#         will be checked for any non-convertible values)
#     :return: Data in long format (ie CELL_LINE, GENE_ID, VALUE)
#     """
#
#     # Currently data is in format where gene ids are in rows and cell line ids are in columns -- to correct this,
#     # first ensure that the two gene identifiers given do not contain any conflicts:
#     assert d.groupby('COMMON')['GENE_ID'].nunique().max() == 1
#     assert d.groupby('GENE_ID')['COMMON'].nunique().max() == 1
#
#     # Next, conform the gene id names to something more informative and then pivot
#     # the frame so that the gene ids are in rows along with cell line ids (pushed down out of columns):
#     d = d.rename(columns={'GENE_ID': 'GENE_ID:CGDS', 'COMMON': 'GENE_ID:HGNC'})
#     d = pd.melt(d, id_vars=['GENE_ID:CGDS', 'GENE_ID:HGNC'], var_name='CELL_LINE_ID', value_name='VALUE')
#
#     # Ensure all non-value fields are non-null
#     assert np.all(d[[c for c in d if c != 'VALUE']].notnull())
#
#     # If specified drop rows with null values
#     if rm_null_values:
#         # Create mask for values to be ignored noting that the CGDS web service sometimes
#         # returns string null for values
#         mask = (d['VALUE'].notnull()) & (d['VALUE'].astype(str) != 'null')
#
#         # Drop any records containing null values for any cell line + gene combination
#         d = subset(d, lambda df: df[mask], subset_op='Remove null values for column "VALUE"')
#
#         # Ensure no values are null
#         assert np.all(d['VALUE'].notnull())
#
#     # Cast value to numeric unless disabled (which should rarely be the case -- CGDS table results
#     # are almost always numeric)
#     if value_numeric:
#         d['VALUE'] = pd.to_numeric(d['VALUE'])
#         assert np.all(d['VALUE'].notnull())
#
#     return d


# DEFAULT_IGNORABLE_MUTATION_COLS =  [
#     'mutation_status',
#     'validation_status',
#     'xvar_link', 'xvar_link_pdb', 'xvar_link_msa',
#     'reference_read_count_normal',
#     'variant_read_count_normal'
# ]
# def prep_mutation_data(d, c_rm=DEFAULT_IGNORABLE_MUTATION_COLS):
#     d_exp = d.drop(c_rm, axis=1).copy()
#
#     # Convert Entrez ID to integer (and make sure this doesn't some how lead to differences)
#     assert np.all(d_exp['entrez_gene_id'] == d_exp['entrez_gene_id'].astype(np.int64))
#     d_exp['entrez_gene_id'] = d_exp['entrez_gene_id'].astype(np.int64)
#
#     # Fill "Functional Impact Score" NAs with unknown (it is always one of 'L', 'M', 'N', or 'H')
#     # d_exp['functional_impact_score'] = d_exp['functional_impact_score'].fillna('Unknown')
#     # d_exp['sequencing_center'] = d_exp['sequencing_center'].fillna('Unknown')
#
#     # Conform field names and convert to ucase
#     d_exp = d_exp.rename(columns={
#         'entrez_gene_id': 'GENE_ID:ENTREZ',
#         'gene_symbol': 'GENE_ID:HGNC',
#         'case_id': 'CELL_LINE_ID'
#     })
#     d_exp = d_exp.rename(columns=lambda c: c.upper())
#
#     # Remove completely duplicated records
#     d_exp = subset(d_exp, lambda df: df[~df.duplicated()], subset_op='Remove duplicate records')
#
#     return d_exp

