from itertools import groupby
from pathlib import Path
from sys import argv
from time import perf_counter
from statistics import median, stdev
from itaxotools.taxi3.align import *
from itaxotools.taxi3.distances import *
from itaxotools.taxi3.pairs import *
from itaxotools.taxi3.sequences import *
from itaxotools.taxi3.partitions import *
from itaxotools.taxi3.handlers import *
from itaxotools.taxi3.handlers import *
from itaxotools.taxi3.statistics import *
from typing import NamedTuple
import numpy as np
import os

class SequenceInfo(NamedTuple):
    type: str
    file_size: int


class TabularSequenceInfo(SequenceInfo):
    headers: list[str]


class PartitionInfo(NamedTuple):
    type: str
    file_size: int


class TabularPartitionInfo(PartitionInfo):
    headers: list[str]


class SpartitionInfo(PartitionInfo):
    spartitions: list[str]


class Results(NamedTuple):
    stats_all: Path | None
    stats_species: Path | None
    stats_genus: Path | None

    aligned_pairs: Path | None

    summary: Path

    ALOTOFFILES: ...

    matrices: dict[str, Path]

    time_taken: ...

    number_of_files_created: int


class VersusAll:

    def __init__(self):

        self.work_dir: Path
        self.progress_handler: Callable

        self.input_sequences: Path
        self.input_columns: tuple[str]
        self.input_columns_res: tuple[str]

        self.perform_species_analysis: bool
        self.species_path: Path
        self.species_columns: tuple[str]
        self.species_columns_res: tuple[str]
        self.species_spartition: str

        self.perform_genus_analysis: bool
        self.genera_path: Path
        self.genera_columns: tuple[str]
        self.genera_columns_res: tuple[str]
        self.genera_spartition: str

        self.do_pairwise_alignment: bool
        self.pairwise_scores: Scores

        self.distance_metrics: list[DistanceMetric]
        self.distance_linear_output: bool
        self.distance_matrix_output: bool
        self.distances_percentile: bool
        self.distances_formatter: str = '{:.4f}'
        self.distances_missing: str = 'NA'

        self.stats_all: bool
        self.stats_species: bool
        self.stats_genus: bool

    def set_sequence_path(self, path: Path) -> SequenceInfo:
        pass

    def set_species_path(self, path: Path) -> PartitionInfo:
        pass

    def set_genus_path(self, path: Path) -> PartitionInfo:
        pass

    def createDir(self):
        fileDir = ['Statistics', 'Pairs', 'Distances', 'Summary', 'Species', 'Genus']
        try:
            for directory in fileDir:
                dirPath = os.path.join(self.work_dir, directory)
                os.mkdir(dirPath)
        except OSError as error:
            raise error

    def start(self, progress_handler, warnings_handler) -> Results:
        # Construct conditional pipeline
        # execute it
        # return

        #Create all Dir
        self.createDir()



        statsCalculator = StatisticsCalculator()
        if self.stats_all:

            if self.stats_species:
                pass
            if self.stats_genus:
                pass


        if self.perform_species_analysis:
            pass
        if self.perform_genus_analysis:
            pass

        if self.do_pairwise_alignment:
            pass

        #Distance Metrics

        if self.distance_linear_output:
            if self.distances_percentile:
                pass
            else:
                USE = self.distances_formatter
            pass

        if self.distance_matrix_output:
            if self.distances_percentile:
                pass
            else:
                USE = self.distances_formatter
            pass



        pass
