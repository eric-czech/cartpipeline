
from pycgds import tcga as tcga_data
import pandas as pd

def add_args(parser):
    parser.add_argument(
        '--gene-meta-path',
        required=True,
        metavar='PATH',
        help='Path to CSV file containing gene metadata to be used for expression data collection'
    )
    parser.add_argument(
        '--cache-dir',
        metavar='DIR',
        help='Path to cache directory for all CGDS web service calls '
             '(ensures that repeat calls to this command are faster)'
    )
    parser.add_argument(
        '--study-id',
        required=True,
        metavar='STUDYID',
        help='Name of TCGA cohort/study id (eg prad_tcga or prad_tcga_pub)'
    )
    parser.add_argument(
        '--use-rna-seq',
        default=True,
        metavar='True|False',
        help='Flag indicating whether to collect RNA-seq or microarray expression data (defaults to RNA-seq)'
    )
    return parser


def get_expression_data(args):
    gene_meta_path = args.gene_meta_path
    use_rna_seq = args.use_rna_seq
    study_id = args.study_id
    cache_dir = args.cache_dir

    data_type = tcga_data.DATA_TYPE_RNASEQ_ZSCORE if use_rna_seq else tcga_data.DATA_TYPE_EXPRESSION_ZSCORE
    gene_list = pd.read_csv(gene_meta_path)['Gene'].unique()

    df = tcga_data.get_data([study_id], data_type, gene_list, cache_dir=cache_dir)

    return df

