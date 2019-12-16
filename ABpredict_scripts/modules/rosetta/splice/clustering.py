from . import splice_util
from .splice_design import SpliceDesign
from flab.cluster.cluster import Cluster

from collections import OrderedDict



def get_cluster_of_design(design, segments: OrderedDict) -> str:
    """

    :param design: path to splice design or a SpliceDesign object
    :param segments: ordered dict of [segment] = MaxCluster object
    :return: clustering, e.g. 14.36 (14 for segment 1, 36 for segment 2)
    """
    if str == type(design):
        design = SpliceDesign(path=design)
    return '.'.join(
        [str(clustering.get_protein_cluster(design.get_source(segment)))
         for segment, clustering in segments.items()])


def cluster_splice_structures(splice: list, segments: OrderedDict,
                              print_=True) -> Cluster:
    """

    :param splice: list of SpliceDesigns or pdb paths
    :param segments: ordered dict of [segment] = MaxCluster object. Assumes that
    the protein names in the Cluster objects are the same as in the segments
    part of the splice design pdb
    :param print_: whether to print the progress or not
    :return: Clustering of designs (pdb opath) by the segments clustering
    """
    clusters = Cluster()
    counter = 1
    for pdb in splice:
        cluster = get_cluster_of_design(pdb, segments)
        clusters.add(cluster, pdb)
        if print_ and 0 == (counter % 1000):
            print(counter)
        counter += 1

    return clusters


