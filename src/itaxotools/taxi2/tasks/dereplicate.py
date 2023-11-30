from __future__ import annotations

from itertools import chain, groupby
from pathlib import Path
from time import perf_counter
from typing import Callable, Iterator, Literal, NamedTuple

from itaxotools.common.utility import AttrDict

from ..align import PairwiseAligner, Scores
from ..distances import Distance, DistanceHandler, DistanceMetric
from ..file_types import FileFormat
from ..files import identify_format
from ..handlers import FileHandler, ReadHandle, WriteHandle
from ..pairs import SequencePair, SequencePairHandler, SequencePairs
from ..sequences import Sequence, SequenceHandler, Sequences


def multiply(iterator: iter, n: int):
    return (item for item in iterator for i in range(n))


def split(source: iter, *funcs: list[Callable]):
    source = multiply(source, len(funcs))
    return [map(func, source) for func in funcs]


def console_report(caption, index, total):
    if caption == "Finalizing...":
        print(f"\rCalculating... {total}/{total} = {100:.2f}%", end="")
        print("\nFinalizing...")
    else:
        print(f"\rCalculating... {index}/{total} = {100*index/total:.2f}%", end="")


class AllInfo(NamedTuple):
    query: Sequence
    id_x: str
    id_y: str
    len_x: int
    len_y: int
    distance: float
    similar: bool


class SummaryLine(NamedTuple):
    query_id: str
    query_length: str
    included_id: int
    included_length: int
    included_distance: float
    excluded_id: int
    excluded_length: int
    excluded_distance: float


class Results(NamedTuple):
    output_directory: Path
    seconds_taken: float


class SummaryHandle(FileHandler[SummaryLine]):
    def _open(
        self,
        path: Path,
        mode: Literal["r", "w"] = "w",
        missing: str = "NA",
        formatter: str = "{:f}",
        *args,
        **kwargs,
    ):
        self.missing = missing
        self.formatter = formatter
        super()._open(path, mode, *args, **kwargs)

    def distance_to_text(self, d: float | None) -> str:
        if d is None:
            return self.missing
        return self.formatter.format(d)

    def _iter_read(self, *args, **kwargs) -> ReadHandle[SummaryLine]:
        raise NotImplementedError()

    def format_line(self, line: SummaryLine):
        return (
            line.query_id,
            str(line.query_length),
            line.included_id,
            str(line.included_length),
            self.distance_to_text(line.included_distance),
            line.excluded_id,
            str(line.excluded_length),
            self.distance_to_text(line.excluded_distance),
        )

    def _iter_write(self, *args, **kwargs) -> WriteHandle[SummaryLine]:
        try:
            headers = SummaryLine._fields
            with FileHandler.Tabfile(self.path, "w", columns=headers) as file:
                while True:
                    line = yield
                    file.write(self.format_line(line))
        except GeneratorExit:
            return


class Dereplicate:
    def __init__(self):
        self.work_dir: Path = None
        self.paths = AttrDict()

        self.progress_handler: Callable = console_report
        self.progress_interval: float = 0.015

        self.input: Sequences = None
        self.output_format: FileFormat = None
        self.excluded: set[str] = set()

        self.params = AttrDict()

        self.params.thresholds = AttrDict()
        self.params.thresholds.similarity: float = 0.07
        self.params.thresholds.length: int = 10

        self.params.pairs = AttrDict()
        self.params.pairs.align: bool = True
        self.params.pairs.write: bool = True
        self.params.pairs.scores: Scores = None

        self.params.distances = AttrDict()
        self.params.distances.metric: DistanceMetric = None
        self.params.distances.write_linear: bool = True
        self.params.distances.write_matricial: bool = True

        self.params.format = AttrDict()
        self.params.format.float: str = "{:.4f}"
        self.params.format.missing: str = "NA"
        self.params.format.percentage_multiply: bool = False

    def set_output_format_from_path(self, path: Path):
        self.output_format = identify_format(path)

    def get_output_handler(self, path: Path):
        if self.output_format == FileFormat.Fasta:
            return SequenceHandler.Fasta(path, "w", write_organism=True)
        if self.output_format == FileFormat.Tabfile:
            return SequenceHandler.Tabfile(
                path, "w", idHeader="seqid", seqHeader="sequence"
            )
        raise Exception("Unknown file format")

    def check_params(self):
        self.output_format = self.output_format or FileFormat.Tabfile
        self.params.distances.metric = (
            self.params.distances.metric or DistanceMetric.Uncorrected()
        )

    def generate_paths(self):
        assert self.work_dir
        self.create_parents(self.work_dir)
        metric = str(self.params.distances.metric)
        extension = self.output_format.extension

        self.paths.summary = self.work_dir / "summary.tsv"
        self.paths.dereplicated = self.work_dir / f"dereplicated{extension}"
        self.paths.excluded = self.work_dir / f"excluded{extension}"
        self.paths.aligned_pairs = self.work_dir / "aligned_pairs.txt"
        self.paths.distances_linear = (
            self.work_dir / "distances" / f"{metric}.linear.tsv"
        )
        self.paths.distances_matricial = (
            self.work_dir / "distances" / f"{metric}.matricial.tsv"
        )

    def create_parents(self, path: Path):
        if path.suffix:
            path = path.parent
        path.mkdir(parents=True, exist_ok=True)

    def drop_short_sequences(self, sequences: iter[Sequence]) -> Iterator[Sequence]:
        for sequence in sequences:
            if len(sequence.seq) >= self.params.thresholds.length:
                yield sequence

    def drop_identical_pairs(self, pairs: iter[SequencePair]) -> Iterator[SequencePair]:
        for pair in pairs:
            if pair.x.id != pair.y.id:
                yield pair

    def drop_excluded_pairs(self, pairs: iter[SequencePair]) -> Iterator[SequencePair]:
        for pair in pairs:
            if all(
                (
                    pair.x.id not in self.excluded,
                    pair.y.id not in self.excluded,
                )
            ):
                yield pair

    def normalize_pairs(self, pairs: iter[SequencePair]) -> Iterator[SequencePair]:
        if not self.params.pairs.align:
            yield from pairs
            return

        for pair in pairs:
            yield SequencePair(pair.x.normalize(), pair.y.normalize())

    def align_pairs(self, pairs: iter[SequencePair]) -> Iterator[SequencePair]:
        if not self.params.pairs.align:
            yield from pairs
            return

        aligner = PairwiseAligner.Biopython(self.params.pairs.scores)
        yield from aligner.align_pairs(pairs)

    def write_pairs(self, pairs: iter[SequencePair]) -> Iterator[SequencePair]:
        if not self.params.pairs.write:
            yield from pairs
            return

        self.create_parents(self.paths.aligned_pairs)
        with SequencePairHandler.Formatted(self.paths.aligned_pairs, "w") as file:
            for pair in pairs:
                file.write(pair)
                yield pair

    def calculate_distances(self, pairs: iter[SequencePair]) -> Iterator[Distance]:
        metric = self.params.distances.metric
        for x, y in pairs:
            yield metric.calculate(x, y)

    def adjust_distances(self, distances: iter[Distance]) -> Iterator[Distance]:
        if not self.params.format.percentage_multiply:
            yield from distances
            return

        for distance in distances:
            if distance.d is not None:
                distance = distance._replace(d=distance.d * 100)
            yield distance

    def write_distances_linear(self, distances: iter[Distance]) -> Iterator[Distance]:
        if not self.params.distances.write_linear:
            yield from distances
            return

        self.create_parents(self.paths.distances_linear)
        with DistanceHandler.Linear.WithExtras(
            self.paths.distances_linear,
            "w",
            missing=self.params.format.missing,
            formatter=self.params.format.float,
        ) as file:
            for distance in distances:
                file.write(distance)
                yield distance

    def write_distances_matrix(self, distances: iter[Distance]) -> Iterator[Distance]:
        if not self.params.distances.write_matricial:
            yield from distances
            return

        self.create_parents(self.paths.distances_matricial)
        with DistanceHandler.Matrix(
            self.paths.distances_matricial,
            "w",
            missing=self.params.format.missing,
            formatter=self.params.format.float,
        ) as file:
            for distance in distances:
                file.write(distance)
                yield distance

    def check_similar_distances(self, distances: iter[Distance]) -> Iterator[bool]:
        similarity = self.params.thresholds.similarity
        for distance in distances:
            if distance.d is None:
                yield False
            else:
                yield bool(distance.d <= similarity)

    def combine_all_info(
        self,
        pairs: iter[SequencePair],
        distances: iter[Distance],
        are_similar: iter[bool],
    ) -> Iterator[AllInfo]:
        all = zip(pairs, distances, are_similar)
        return (
            AllInfo(
                query=pair.x,
                id_x=pair.x.id,
                id_y=pair.y.id,
                len_x=len(pair.x.seq),
                len_y=len(pair.y.seq),
                distance=distance.d,
                similar=is_similar,
            )
            for pair, distance, is_similar in all
        )

    def group_all_left(self, infos: iter[AllInfo]) -> Iterator[iter[AllInfo]]:
        for _, group in groupby(infos, lambda info: info.id_x):
            yield group

    def find_replicates(self, groups: iter[iter[AllInfo]]) -> Iterator[SummaryLine]:
        """From each group of replicate sequences, we only keep the longest"""

        for infos in groups:
            first = next(infos)
            query_id = first.id_x
            query_length = first.len_x
            max_id = first.id_x
            max_length = first.len_x
            max_distance = first.distance

            for _, _, id_y, _, len_y, distance, similar in chain([first], infos):
                if not similar:
                    continue

                if len_y > max_length:
                    included_id = id_y
                    included_length = len_y
                    included_distance = distance
                    excluded_id = max_id
                    excluded_length = max_length
                    excluded_distance = max_distance
                else:
                    included_id = max_id
                    included_length = max_length
                    included_distance = max_distance
                    excluded_id = id_y
                    excluded_length = len_y
                    excluded_distance = distance

                self.excluded.add(excluded_id)
                yield SummaryLine(
                    query_id=query_id,
                    query_length=query_length,
                    included_id=included_id,
                    included_length=included_length,
                    included_distance=included_distance,
                    excluded_id=excluded_id,
                    excluded_length=excluded_length,
                    excluded_distance=excluded_distance,
                )

                if len_y > max_length:
                    max_id = id_y
                    max_length = len_y
                    max_distance = distance

    def write_summary(self, lines: iter[SummaryLine]) -> Iterator[SummaryLine]:
        with SummaryHandle(
            self.paths.summary,
            "w",
            missing=self.params.format.missing,
            formatter=self.params.format.float,
        ) as file:
            for line in lines:
                if line is not None:
                    file.write(line)
                yield line

    def write_file_dereplicated(self, sequences: iter[Sequence]) -> Iterator[Sequence]:
        with self.get_output_handler(self.paths.dereplicated) as file:
            for sequence in sequences:
                if sequence.id not in self.excluded:
                    file.write(sequence)
                yield sequence

    def write_file_excluded(self, sequences: iter[Sequence]) -> Iterator[Sequence]:
        with self.get_output_handler(self.paths.excluded) as file:
            for sequence in sequences:
                if sequence.id in self.excluded:
                    file.write(sequence)
                yield sequence

    def report_progress(self, distances: iter[Distance], data):
        section = len(data)
        total = len(data) ** 2
        last_time = perf_counter()
        for index, distance in enumerate(distances, 1):
            new_time = perf_counter()
            if new_time - last_time >= self.progress_interval:
                self.progress_handler(
                    "distance.x.id", index, total - len(self.excluded) * section
                )
                last_time = new_time
            yield distance
        self.progress_handler("Finalizing...", total, total)

    def start(self) -> None:
        ts = perf_counter()

        self.excluded = set()
        self.check_params()
        self.generate_paths()

        data = Sequences(self.drop_short_sequences, self.input)

        pairs = SequencePairs.fromProduct(data, data)
        pairs = self.drop_identical_pairs(pairs)
        pairs = self.drop_excluded_pairs(pairs)

        pairs = multiply(pairs, 2)
        unaligned_pairs = pairs

        aligned_pairs = self.normalize_pairs(pairs)
        aligned_pairs = self.align_pairs(aligned_pairs)
        aligned_pairs = self.write_pairs(aligned_pairs)

        distances = self.calculate_distances(aligned_pairs)
        distances = self.report_progress(distances, data)

        distances = self.adjust_distances(distances)
        distances = self.write_distances_linear(distances)
        distances = self.write_distances_matrix(distances)

        distances = multiply(distances, 2)
        are_similar = self.check_similar_distances(distances)

        all = self.combine_all_info(unaligned_pairs, distances, are_similar)

        groups = self.group_all_left(all)
        lines = self.find_replicates(groups)
        lines = self.write_summary(lines)

        for _ in lines:
            pass

        data = self.write_file_dereplicated(data)
        data = self.write_file_excluded(data)

        for _ in data:
            pass

        tf = perf_counter()

        return Results(self.work_dir, tf - ts)
