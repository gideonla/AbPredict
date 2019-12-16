#!/bin/env python3
import logging
import re


logger = logging.getLogger(__name__)


class RosettaXML:
    def __init__(self, path: str):
        """
        TODO: probably some parsing.... 
        """
        self.path = path

    def get_script_vars(self) -> list:
        return list(set([i.replace('%', '') 
                for i in re.findall('%%.+?%%', open(self.path).read())]))

