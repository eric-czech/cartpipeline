import logging
from io import StringIO
from argparse import ArgumentParser
from pyagg.aggregation import add_args, aggregate_pipeline_results

logger = logging.getLogger(__name__)


def make_arg_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--output",
        required=True,
        metavar='PATH',
        help="Name of CSV file to contain aggregated pipeline data"
    )
    return parser

if __name__ == "__main__":
    # Parse arguments
    parser = add_args(make_arg_parser())
    args = parser.parse_args()
    logger.info('Pipeline aggregation arguments: {}'.format(args))

    # Run selection
    df = aggregate_pipeline_results(args)

    # Print result info
    info = StringIO()
    df.info(buf=info)
    logging.info('Pipeline aggregation result info:\n{}'.format(info.getvalue()))

    # Write results to file
    df.to_csv(args.output, index=False)