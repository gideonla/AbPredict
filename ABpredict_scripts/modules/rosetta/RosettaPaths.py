"""
a class to keep Rosetta related paths
"""
import os


class RosettaPaths():
    """RosettaPaths"""

    def __init__(self,
                 rosetta_dir: str=None,
                 static_database: str=None,
                 use_static: bool=False):
        """__init__

        :param rosetta_dir: a dir for the Rosetta repo
        :type rosetta_dir: str
        """
        self.use_static = use_static
        if rosetta_dir:
            self.rosetta_dir = os.path.abspath(os.path.expanduser(rosetta_dir))
        else:
            self.rosetta_dir = os.path.abspath(os.path.expanduser('~/Rosetta'))
        self.source = '%s/main/source/' % self.rosetta_dir
        self.bin = '%s/bin/' % self.source
        self.rosetta_scripts_default = '%s/rosetta_scripts.default.linuxgccrelease' % self.bin
        self.rosetta_static_default = '%s/rosetta_scripts.static.linuxgccrelease' % rosetta_dir
        self.database = '%s/main/database' % self.rosetta_dir
        self.fragment_picker = '%s/fragment_picker.default.linuxgccrelease' % self.bin
        self.symm_dir = '%s/src/apps/public/symmetry/' % self.source
        self.make_symm_denovo = '%s/make_symmdef_file_denovo.py' % self.symm_dir
        self.make_symmdef = '%s/make_symmdef_file.pl' % self.symm_dir
        self.static_databse = static_database
        self.protein_tools = '%s/tools/protein_tools/' % self.rosetta_dir

    def r_scripts(self) -> str:
        print (self.rosetta_static_default)
        return self.rosetta_scripts_default if not self.use_static else self.rosetta_static_default

    def db(self) -> str:
        return self.database if not self.use_static else self.static_databse
