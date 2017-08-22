
import pandas as pd
import numpy as np
import logging
logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument(
        '--gene-meta-path',
        required=True,
        metavar='PATH',
        help='Path to CSV file containing gene metadata to be used for expression data collection'
    )
    parser.add_argument(
        '--gene-exp-paths',
        required=True,
        nargs='+',
        metavar='PATH',
        help='Path(s) to CSV file containing TCGA expression data'
    )
    return parser


def get_gene_meta(path):

    # Read in raw gene meta data frame and subset to relevant fields
    d = pd.read_csv(path)
    d = d[[
        'Gene', 'Gene synonym', 'Ensembl', 'Chromosome',
        'RNA tissue category', 'RNA TS', 'RNA TS TPM',
        'Protein classes'
    ]].rename(columns=lambda c: c.title().replace(' ', ''))

    # Resolve duplicated genes by taking only records for duplicates when Gene Synonym
    # is also present; At TOW, this was an effective way to remove unnecessarily repeated
    # records as the ones with missing synonyms were less informative
    cts = d['Gene'].value_counts()
    dupe_genes = cts[cts>1].index.values
    if len(dupe_genes) > 1:
        logger.warning('Removing duplicated records for the following genes: "{}"'.format(dupe_genes))
        d = d[~d['Gene'].isin(dupe_genes) | d['GeneSynonym'].notnull()]

    # Ensure that an error was not made in the above assumption, that gene synonym non-existence
    # was sufficient to remove all duplicates (if not, this will require further review)
    assert not d['Gene'].duplicated().any()

    return d


def get_exp_stats(path):
    d = pd.read_csv(path)

    # It is assumed that this data never contains duplicates, but verify that here to be sure
    assert not d[['StudyId', 'Gene', 'SampleId']].duplicated().any().any()

    # Group by study name and gene and calculate statistics for expression levels
    d = (
        d.groupby(['StudyId', 'Gene'])['Value']
        .describe(percentiles=list(np.arange(.1, 1, .1)) + [.95, .99])
    )

    # For some reason, percentiles occasionally have .0 suffixed to them due to rounding (so remove that)
    d = d.rename(columns=lambda c: c.replace('.0', '').title())

    return d


def merge(d_gene, d_exp):
    # Combine gene meta data with expression data, merging on gene symbol (not id of some kind)
    # and for now, ignore any matches from either side (inner join)
    # TODO: Analyze unjoined genes from these data sets to see if any are recoverable
    d = pd.merge(
        d_gene.set_index('Gene').add_prefix('Meta:').reset_index(),
        d_exp.add_prefix('Stat:').reset_index(level='StudyId').reset_index(),
        how='inner',
        on='Gene'
    )
    return d


def aggregate_pipeline_results(args):
    d_gene = get_gene_meta(args.gene_meta_path)
    d_exp = pd.concat([get_exp_stats(path) for path in args.gene_exp_paths])
    return merge(d_gene, d_exp)
