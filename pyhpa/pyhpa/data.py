import re
from collections import Iterable
import pandas as pd
import numpy as np

HPA_URL_FMT = 'http://v{}.proteinatlas.org/download/proteinatlas.tab.gz'
REGEX_SUB_PARENS = re.compile('\(.*\)')


def get_hpa_data(versions):
    """
    Download HPA data by version

    Data for this source currently available at "http://v<version_number>.proteinatlas.org/download/proteinatlas.tab.gz"

     * At TOW, the following versions were available: 13, 14, 15, 16

    :param versions: String or sequence of strings or integers indicating version numbers to collect data for
    :return: Data frame with a "Version" field reflecting which version the rows in the frame correspond to
    """
    if isinstance(versions, str):
        versions = [versions]
    if not isinstance(versions, Iterable):
        versions = [versions]

    data = []
    for v in versions:
        url = HPA_URL_FMT.format(v)
        try:
            data.append(pd.read_csv(url, compression='gzip', sep='\t').assign(Version=v))
        except Exception as e:
            raise ValueError('Failed to retrieve data for HPA version {} (URL = "{}")'.format(v, url)) from e

    return pd.concat(data)


def prepare_hpa_data(d):
    """
    Prepare raw HPA data by normalizing (some) differences across and within versions

    :param d: Data frame from `get_hpa_data`
    :return: Prepared data frame
    """
    # Prep protein class lists
    d['Protein classes'] = (
        # Split class string on commas (missing protein class will become empty tuple)
        d['Protein class'].fillna('').str.split(',')
        # Apply clean function to each item is split sequences and deduplicate tuple
        .apply(lambda v: tuple(set(map(clean_protein_class_name, v))))
    )
    return d.drop('Protein class', axis=1)


def filter_by_protein_class(d, protein_classes):
    """
    Restrict data to a desired set of protein classes

    When multiple classes are present, records are returned when that list contains at least one of the classes given.

    At TOW, the allowed list of protein class names available for filtering on are (with these exact spellings):
        - Predicted intracellular proteins
        - Predicted membrane proteins
        - Plasma proteins
        - Protein evidence
        - Disease related genes
        - Enzymes
        - Predicted secreted proteins
        - Cancer-related genes
        - Transcription factors
        - Potential drug targets
        - Transporters
        - G-protein coupled receptors
        - FDA approved drug targets
        - Mitochondrial proteins
        - CD markers
        - Cytoskeleton related proteins
        - RAS pathway related proteins
        - Candidate cardiovascular disease genes
        - Ribosomal proteins
        - Voltage-gated ion channels
        - Nuclear receptors
        - Blood group antigen proteins
        - RNA polymerase related proteins
        - Citric acid cycle related proteins

    If any of the class names given to filter on are not present in the data provided, an error will be thrown.

    :param d: Data frame from `prepare_hpa_data`
    :param protein_classes: List of protein class names for which data should be restricted to
    :return: Filtered data frame
    """
    if 'Protein classes' not in d:
        raise ValueError(
            'Filtering can only be executed on prepared HPA data; Make sure data given is from `prepare_hpa_data` '
            'function and not in raw form [Cause = Missing field "Protein classes"]'
        )

    # Given how protein class lists are generated, it should not be possible for them
    # to be null but run sanity check here anyhow
    assert not d['Protein classes'].isnull().any()

    # Validate that all of the protein class names given are actually real protein classes;
    # otherwise these values were probably misspelled or should be removed
    uniq_classes = pd.Series([v for l in d['Protein classes'] for v in l]).unique()
    diff_classes = np.setdiff1d(protein_classes, uniq_classes)
    if len(diff_classes) > 0:
        raise ValueError('The following protein class filters provided do not exist in ANY HPA records: {}'
                         .format(diff_classes))

    # Apply filter
    mask = d['Protein classes'].apply(lambda v: len(np.intersect1d(v, protein_classes)) > 0)
    return d[mask]


def clean_protein_class_name(name):
    """
    Return processed protein class name

    Common string equivalence found in these names based on removal of leading/trailing whitespace
    and parenthetical enclosures.

    Examples:
        - "Protein evidence (Ezkurdia et al 2014)" -> "Protein evidence"
        - "Disease related genes " -> "Disease related genes"
        - " FDA approved drug targets" -> "FDA approved drug targets"

    :param name: Protein class name as string
    :return: Clean string name
    """
    v = REGEX_SUB_PARENS.sub('', name).strip()
    return v


