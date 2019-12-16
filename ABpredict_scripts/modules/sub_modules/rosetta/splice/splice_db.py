from ...pdb_.util import *
from Bio.SeqUtils import seq1, seq3
import os



class SpliceDB:
    def __init__(self, db):
        self.original_path = db if os.path.isfile(db) else ''
        # self.new_path = new_path if new_path else self.original_path

        self.db = self.parse_db(db)
        self.start_res = int(self.db[-1][0])
        self.end_res = int(self.db[-1][1])
        self.chainbreak = int(self.db[-1][2])
        self.protein_name = self.db[-1][3]

    def parse_db(self, db):
        """

        :param db:
        :return: list of lists [phi, psi, omega, res]. The last list is
        [start res, end res, chainbreak, protein name]
        """
        if os.path.isfile(db):
            db = open(db).read().split()
        else:
            db = db.split()
        db = [db[i:i + 4] for i in range(0, len(db), 4)]
        return db

    def update_db(self, design):
        design = parse_pdb_file(design)
        chain = list(design.get_chains())[0].id
        design_residues = design[0][chain]
        new_db = list()
        index = self.start_res
        print(self.protein_name)
        for phi, psi, omega, res in self.db[:-1]:
            new_res = design_residues[index].resname
            new_db.append([phi, psi, omega, new_res])
            index += 1
            if res != new_res:
                print(index - 1, res, new_res)
        new_db.append(self.db[-1])  # adding the start, end, chainbreak and name line
        self.db = new_db

    def write_db(self, path=''):
        """

        :param path: path to write the new database to
        :return:
        """
        if not path:
            path = self.protein_name + '_new.db'
        open(path, 'w').write(' '.join([' '.join(res) for res in self.db]) + '\n')

    def get_sequence(self):
        """Returns the sequence of the database (the forth element of each
        residue"""
        seq = ''.join([seq1(pos[-1]) for pos in self.db[:-1]])
        return seq

