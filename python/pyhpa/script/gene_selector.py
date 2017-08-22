
import logging
from io import StringIO
from argparse import ArgumentParser
from pyhpa.gene_selector import select_genes, add_args

logger = logging.getLogger(__name__)


def make_arg_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--output",
        required=True,
        metavar='PATH',
        help="Name of CSV file to contain resulting selected gene/protein metadata"
    )
    return parser

if __name__ == "__main__":
    # Parse arguments
    parser = add_args(make_arg_parser())
    args = parser.parse_args()
    logger.info('Gene selection arguments: {}'.format(args))

    # Run selection
    df = select_genes(args)

    # Print result info
    info = StringIO()
    df.info(buf=info)
    logging.info('Gene selection result info:\n{}'.format(info.getvalue()))

    # Write results to file
    df.to_csv(args.output, index=False)