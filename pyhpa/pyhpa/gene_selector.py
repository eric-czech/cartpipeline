
from pyhpa import data as hpa_data


def add_args(parser):
    parser.add_argument(
        "--hpa-version",
        default=16,
        type=int,
        help="HPA version number"
    )
    # parser.add_argument(
    #     "protein_classes",
    #     metavar='C',
    #     type=str,
    #     nargs='+',
    #     help="Protein class names to filter on"
    # )
    return parser


def select_genes(args):
    hpa_version = args.hpa_version
    hpa_protein_classes = [
        'FDA approved drug targets', 'Predicted membrane proteins',
        'Cancer-related genes', 'CD markers'
    ]
    # hpa_protein_classes = [
    #     'FDA approved drug targets'
    # ]
    df = hpa_data.get_hpa_data(hpa_version)
    df = hpa_data.prepare_hpa_data(df)
    df = hpa_data.filter_by_protein_class(df, hpa_protein_classes)
    return df
