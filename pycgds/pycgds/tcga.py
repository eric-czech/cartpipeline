
import pandas as pd
from pycgds import api
import logging
logger = logging.getLogger(__name__)


def get_data(tcga_study_ids, data_type, gene_ids, batch_size=50, cache_dir=None):
    """
    Fetch TCGA study data using the cBioPortal (aka CGDS) web service

    :param tcga_study_ids: List of TCGA-related cancer study identifiers (e.g. "kich_tcga", "ucec_tcga_pub")
    :param data_type: Type of genetic data to fetch (see module constants prefixed by DATA_TYPE.* for possibilities)
    :param gene_ids: List of gene names to collect data for
    :param batch_size: Size of batch requests (in terms of genes) submitted to CGDS
    :param cache_dir: Optional location of directory in which unique calls will be cached; if not present,
        no caching will occur
    :return: Data frame containing study id, gene, and value for each TCGA sample
    """

    # Validate the given list of study identifiers (they should all be something like "kich_tcga", "kirc_tcga",
    # "ucec_tcga_pub", etc where tcga is a substring in the id) and they should all be present within CGDS as
    # well as have data for the given genetic data type
    cancer_studies = api.get_cancer_studies()['cancer_study_id'].unique()
    valid_study_ids = []
    for study_id in tcga_study_ids:
        # Check for tcga substring
        if 'tcga' not in study_id:
            raise ValueError('Study Id "{}" is not valid in this context because it does not pertain '
                             'to TCGA studies'.format(study_id))
        # Check that this study actual exists
        if study_id not in cancer_studies:
            raise ValueError('Study Id "{}" is not a known identifier within cBioPortal'.format(study_id))

        # Check that the requested data actually exists for this study
        genetic_profile_id = study_id + '_' + data_type
        genetic_profiles = api.get_genetic_profiles(study_id)['genetic_profile_id'].unique()
        if genetic_profile_id not in genetic_profiles:
            logging.warning('TCGA study id "{}" does not have data for type "{}" so it will '
                            'be ignored'.format(study_id, data_type))
            continue

        valid_study_ids.append(study_id)

    # Raise if none of the given studies could be validated
    if len(valid_study_ids) == 0:
        raise ValueError('No applicable TCGA study ids found for data type "{}" (study id list = "{}")'
                         .format(data_type, tcga_study_ids))

    # Begin data collection for all remaining, valid TCGA studies
    logger.info('Beginning data collection for TCGA study ids: "{}"'.format(valid_study_ids))
    data = []

    for i, study_id in enumerate(valid_study_ids):

        # Create a CGDS genetic profile (eg "lusc_tcga_mrna")
        genetic_profile_id = study_id + '_' + data_type

        # Create a CGDS case list id (eg "lusc_tcga_all")
        case_list_id = study_id + '_all'

        logger.info(
            'Importing data for study "{}" ({} of {}) -> profile="{}"'
            .format(study_id, i + 1, len(tcga_study_ids), genetic_profile_id)
        )

        # Invoke CGDS web service through batch, possibly cached requests (if cache_dir set)
        d = api.get_genetic_profile_data(
            case_list_id, genetic_profile_id, gene_ids,
            batch_size=batch_size, cache_dir=cache_dir
        )

        if len(d) == 0:
            continue
        data.append(d.assign(STUDY_ID=study_id))

    # Raise on empty results (before transformations/assertions that will fail otherwise)
    if len(data) == 0:
        raise ValueError('No data found for study ids = "{}", data type = "{}"'.format(tcga_study_ids, data_type))

    # Concatenate results from all TCGA studies
    d = pd.concat(data).rename(columns={'COMMON': 'GENE'})

    # Assert that gene ids and names are unique to one another
    assert d.groupby(['GENE'])['GENE_ID'].nunique().max() == 1
    assert d.groupby(['GENE_ID'])['GENE'].nunique().max() == 1

    # Assert that there are no duplicates per study + gene
    assert d.groupby(['STUDY_ID', 'GENE_ID', 'GENE']).size().max() == 1

    # Stack result to transfrom sample ids out of columns (ie convert data to long format)
    d = d.set_index(['STUDY_ID', 'GENE_ID', 'GENE'])
    d.columns.name = 'SAMPLE_ID'
    d = d.stack().rename('VALUE')

    # Rename upper-underscore fields and return result
    return d.reset_index().rename(columns=lambda c: c.title().replace('_', ''))


#-------------------------#
# TCGA Specific Constants #
#-------------------------#

# Methylation (HM450)
# --------------------------------
# Description: Methylation (HM450) beta-values for genes in $X cases. For genes with multiple methylation
# probes, the probe least correlated with expression.
DATA_TYPE_METHYLATION = 'methylation_hm450'

# Mutations
# --------------------------------
# Description: Mutation data from whole exome sequencing.
DATA_TYPE_MUTATIONS = 'mutations'


# ### Copy Number ###
# Putative copy-number alterations from GISTIC
# --------------------------------
# Description: Putative copy-number calls on X cases determined using GISTIC 2.0. Values: -2 = homozygous deletion;
#  -1 = hemizygous deletion; 0 = neutral / no change; 1 = gain; 2 = high level amplification.
DATA_TYPE_CNA_PUTATIVE = 'gistic'
# Relative linear copy-number values
# --------------------------------
# Description: Relative linear copy-number values for each gene (from Affymetrix SNP6).
DATA_TYPE_CNA = 'linear_CNA'


# ### RNA Expression ###
# mRNA Expression z-Scores (microarray)
# --------------------------------
# Description: mRNA z-Scores (Agilent microarray) compared to the expression distribution of each gene
# tumors that are diploid for this gene.
DATA_TYPE_EXPRESSION_ZSCORE = 'mrna_median_Zscores'
# mRNA expression (microarray)
# --------------------------------
# Description: Expression levels (Agilent microarray).
DATA_TYPE_EXPRESSION = 'mrna'


# ### RNA-Seq ###
# mRNA Expression z-Scores (RNA Seq V2 RSEM)
# --------------------------------
# Description: mRNA z-Scores (RNA Seq V2 RSEM) compared to the expression distribution of each
# gene tumors that are diploid for this gene.
DATA_TYPE_RNASEQ_ZSCORE = 'rna_seq_v2_mrna_median_Zscores'
# mRNA expression (RNA Seq V2 RSEM)
# --------------------------------
# Description: Expression levels for $X genes in $Y $CANCER_STUDY cases (RNA Seq V2 RSEM).
DATA_TYPE_RNASEQ = 'rna_seq_v2_mrna'


# ### RPPA ###
# Protein expression Z-scores (RPPA)
# --------------------------------
# Description: Protein expression, measured by reverse-phase protein array, Z-scores
DATA_TYPE_ZSCORE = 'rppa_Zscores'
# Protein expression (RPPA)
# --------------------------------
# Description: Protein expression measured by reverse-phase protein array
DATA_TYPE_RPPA = 'rppa'
