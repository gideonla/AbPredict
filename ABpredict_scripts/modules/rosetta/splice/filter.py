import os
from . import clustering
from ...pdb_ import pdb_cluster
from collections import OrderedDict



def find_designs_with_segment(designs: list, segment_name, source):
    """

    :param designs: list of SpliceDesign objects
    :param segment_name:
    :param source:
    :return: list of SpliceDesign objects with the wanted source protein at
    segment "segment_name"
    """
    return [design for design in designs
            if source == design.get_source(segment_name)]


def filter_n_design(designs: list, n=10):
    """

    :param designs: list of SpliceDesign objects
    :return:
    """
    designs.sort(key=lambda x: x.score)
    result = dict()
    flag = False
    for design in designs:
        min_score = 10000000
        for segment in design.get_segments_names():
            source = design.get_source(segment)
            designs_with_segment = find_designs_with_segment(
                                                        designs,
                                                        segment_name=segment,
                                                        source=source)
            # Debug :)
            # if len(designs_with_segment) >= n:
            #     print(source, len(designs_with_segment), os.path.basename(designs_with_segment[n-1].path), designs_with_segment[n-1].score)
            # else:
            #     print(source, len(designs_with_segment))

            # if the current segment has more than N designs with this source in
            # this segment and its n'th design (ordered in ascending energy) has
            # lower energy than previous segment in the current design
            if len(designs_with_segment) >= n and designs_with_segment[n-1].score < min_score:
                flag = True
                min_score = designs_with_segment[n-1].score
                result['design'] = design
                result['segment'] = segment
                result['source'] = source
                # removes all designs with the segment
                result['designs'] = list(set(designs) -
                                         set(designs_with_segment))
                result['n_designs'] = designs_with_segment[:n]

        if flag: # if there was a segment with more than n designs
            return result
    return None


def filter_packstat(designs: list, percent=0.2, pattern='packstat'):
    """

    :param designs: list of SpliceDesign objects
    :param percent:
    :param pattern: pattern of the packstat line in the design
    :return: a list of the top percent designs (according to packstat values)
    """
    for design in designs:
        if not design.packstat:
            design.packstat = design.find_filter_value(path=design.path,
                                                       pattern=pattern)
    designs.sort(key=lambda x: x.packstat, reverse=True)
    return designs[:int(len(designs) * percent)]


def select_best_per_cluster(designs: list, segments: OrderedDict, print=True) -> list:
    """

    :param designs: list of SpliceDesign objects
    :param segments: ordered dict of [segment] = MaxCluster object. Assumes that
    the protein names in the Cluster objects are the same as in the segments
    part of the splice design pdb
    :return: a reduced list of designs - only the best design from each cluster
    """
    splice_paths = [pdb.path for pdb in designs]
    clusters = clustering.cluster_splice_structures(splice_paths, segments)

    best_energy = pdb_cluster.filter_best_energy(clusters)

    if print:
        clusters.save_each_cluster(split=False)
        open('best_energy_per_cluster', 'w').writelines(
            ['{:<8}{:<10}{}\n'.format(protein.cluster,
                                      protein.score,
                                      protein.name)
             for protein in best_energy])

    best_names = [best.name for best in best_energy]
    return [design for design in designs if design.path in best_names]



