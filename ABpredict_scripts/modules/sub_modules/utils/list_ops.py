#! /usr/bin/env python3
"""
commonly used list related operations
"""
from typing import Generator, List


def even_split(lst: List, max_size: int) -> Generator[List, None, None]:
    """Yield successive max_size-sized chunks from lst."""
    for ind in range(0, len(lst), max_size):
        yield lst[ind:ind + max_size]
