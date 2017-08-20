import logging
from io import StringIO
from argparse import ArgumentParser
from pycgds.tcga_expression import add_args, get_expression_data

logger = logging.getLogger(__name__)


def make_arg_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--output",
        required=True,
        metavar='PATH',
        help="Name of CSV file to contain resulting TCGA expression data"
    )
    return parser

if __name__ == "__main__":
    # Parse arguments
    parser = add_args(make_arg_parser())
    args = parser.parse_args()
    logger.info('TCGA expression arguments: {}'.format(args))

    # Run selection
    df = get_expression_data(args)

    # Print result info
    info = StringIO()
    df.info(buf=info)
    logging.info('TCGA expression result info:\n{}'.format(info.getvalue()))

    # Write results to file
    df.to_csv(args.output, index=False)