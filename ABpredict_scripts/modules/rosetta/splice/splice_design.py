from collections import OrderedDict
import flab.pdb_.util as pdb_util
import logging 
import pandas as pd


logger = logging.getLogger(__name__)


class SpliceDesign:
    def __init__(self, path: str, frame='@@@@@', packstat='packstat',
                 set_filters=True, silent=False):
        """

        :param str:
        """
        logger.info('SpliceDesign: ' + path)
        self.path = path
        self.silent = silent
        self.segments = self.get_segments(path, frame)
        if set_filters:
            self.score = float(self.find_filter_value(self.path,
                                                      silent=silent)) 
            self.packstat = float(self.find_filter_value(self.path, 
                                                         pattern=packstat,
                                                         silent=silent))
        self.cluster = None
        self.sequence = (self._filter_value(self.path, 'SEQUENCE:') 
                         if self.silent else 
                         str(pdb_util.get_sequence(self.path)))

    @staticmethod
    def get_segments(pdb, frame='@@@@@'):    
        """
        Extracts splice blade information from comments segment of the pdb
        :param pdb: path to pdb file output of splice
        :param frame: the frames prefix in the segment part of splice (to
        exclude them form the segments)
        :return: ordered dictionary of tuples of [segment name] = chimeric blade
        """
        s = 'segment_'
        return OrderedDict((i[i.find(s) + len(s):].split()[0], 
                            i[i.find(s) + len(s):].split()[-1])
                           for i in open(pdb)
                           if ('segment_' in i) and (frame not in i))

    def get_source(self, segment):
        """

        :param segment:
        :return: the source pdb of the wanted segment
        """
        assert segment in list(self.segments.keys()), \
            'segment {} is not part of design {}'.format(segment, self.path)
        return self.segments[segment]
    
    @staticmethod
    def find_filter_value(path, pattern='', score_location=-1, silent=False):
        """"""
        if not pattern:
            pattern = 'score' if silent else 'pose'
        if silent:
            logger.debug('silent filter value')
            return SpliceDesign._filter_value_silent(path, pattern) 
        else:
            logger.debug('normal filter value')
            return SpliceDesign._filter_value(path, pattern,score_location)  

    @staticmethod 
    def _filter_value(path, pattern='pose', score_location=-1):
        """
        :param pattern: how the filter line starts (default: pose)
        :param score_location: where the score is in the line (default: -1)
        :return: the value of the filter
        """
        pdb = open(path).readlines()
        try:
            line = ([x for x in pdb if x.startswith(pattern)][0]).strip()
            return line.split()[score_location]
        except Exception as e:
            logger.warning('no {} in {}'.format(pattern, path), exc_info=True)
            return None
    
    @staticmethod 
    def _filter_value_silent(path, pattern='score'):
        """Locates the wanted pattern in the SCORE lines. 
        """
        lines = [i.split()[1:] 
                 for i in open(path) if i.startswith('SCORE:')]
        silent_filters = pd.DataFrame.from_records(lines[1:], columns=lines[0])
        return silent_filters.loc[0, pattern]

    def get_segments_names(self):
        return list(self.segments.keys())

    def __repr__(self):
        return self.path


# def similar(des1: SpliceDesign, des2: SpliceDesign) -> bool: 
#     """True if the designs are different only by one segment"""
#     count_different = 0
#     for seg in des1.get_segments_names():
#         if des1.get_source(seg) != des2.get_source(seg):
#             count_different += 1
#     return 1 == count_different
# 
# 
# def __extract_neighbors(designs: list):
#     segments = designs[0].get_segments_names()
#     result = {seg: list() for seg in segments}
#     for des in designs:
#         for seg in segments:
#             result[seg].append(des.get_source(seg))
#     return result
#         
# 
# def related(des1: SpliceDesign, des1_neighbors: list, 
#             des2: SpliceDesign, des2_neighbors: list) -> bool:
#     neighbors1 = __extract_neighbors(des1_neighbors) 
#     neighbors2 = __extract_neighbors(des2_neighbors)
#     segments = des1.get_segments_names()
#     for seg in segments:
#         source1 = des1.get_source(seg)
#         source2 = des2.get_source(seg)
#         if ((source1 == source2) or 
#             ((source1 in neighbors2[seg]) and (source2 in neighbors1[seg]))):
#                 continue
#         else:
#             return False
#     return True
# 

def create_graph(des: list):
    import networkx as nx
    g = nx.Graph()
    for d in des: 
        g.add_node(str(d), design=d)
        for node in g.nodes():
            if similar(d, g.node[node]['design']):
                g.add_edge(node, str(d), color='b')
       
    for n1 in g.node:
        print(n1)
        for n2 in g.node:
            if n1 == n2: continue
            if (not g.edge[n1].keys()) or (not g.edge[n2].keys()): continue
            if ((n1, n2) in g.edges()) or ((n2, n1) in g.edges()): continue
            n1_des = g.node[n1]['design']
            n2_des = g.node[n2]['design']
            n1_neighbors = [g.node[neigh]['design'] 
                            for neigh in g.edge[n1].keys()]
            n2_neighbors = [g.node[neigh]['design'] 
                            for neigh in g.edge[n2].keys()]
            if related(n1_des, n1_neighbors, n2_des, n2_neighbors):
                g.add_edge(n1, n2, color='r')
    return g

