from __future__ import annotations

from itertools import chain, product
from math import inf
from pathlib import Path
from time import perf_counter
from typing import Callable, Generator, Iterator, Literal, NamedTuple, TextIO

from itaxotools.common.utility import AttrDict

from ..align import PairwiseAligner, Scores
from ..distances import Distance, DistanceHandler, DistanceMetric, Distances
from ..handlers import FileHandler, WriteHandle
from ..pairs import SequencePairHandler, SequencePairs
from ..partitions import Partition
from ..plot import ComparisonType, HistogramPlotter
from ..sequences import Sequences
from ..statistics import StatisticsCalculator, StatisticsHandler


def multiply(iterator: iter, n: int):
    return (item for item in iterator for i in range(n))


def console_report(caption, index, total):
    if caption == "Finalizing...":
        print(f"\rCalculating... {total}/{total} = {100:.2f}%", end="")
        print("\nFinalizing...")
    else:
        print(f"\rCalculating... {index}/{total} = {100*index/total:.2f}%", end="")


class GenericDistance(NamedTuple):
    metric: DistanceMetric
    idx: str
    idy: str
    d: float


class SimpleStatistics(NamedTuple):
    min: float
    max: float
    mean: float
    count: int


class DistanceStatistics(NamedTuple):
    metric: DistanceMetric
    idx: str
    idy: str
    min: float
    max: float
    mean: float
    count: int


class SimpleAggregator:
    def __init__(self):
        self.sum = 0.0
        self.min = inf
        self.max = 0.0
        self.n = 0

    def add(self, value: float | None):
        if value is None:
            return
        self.sum += value
        if value < self.min:
            self.min = value
        if value > self.max:
            self.max = value
        self.n += 1

    def calculate(self):
        if not self.n:
            return SimpleStatistics(None, None, None, 0)
        return SimpleStatistics(self.min, self.max, self.sum / self.n, self.n)


class DistanceAggregator:
    def __init__(self, metric: DistanceMetric):
        self.metric = metric
        self.aggs = dict()

    def add(self, d: GenericDistance):
        if (d.idx, d.idy) not in self.aggs:
            self.aggs[(d.idx, d.idy)] = SimpleAggregator()
        self.aggs[(d.idx, d.idy)].add(d.d)

    def __iter__(self) -> Iterator[DistanceStatistics]:
        for (idx, idy), agg in self.aggs.items():
            stats = agg.calculate()
            yield DistanceStatistics(
                self.metric, idx, idy, stats.min, stats.max, stats.mean, stats.count
            )


class SubsetStatisticsHandler(FileHandler[tuple[DistanceStatistics]]):
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

    def distanceToText(self, d: float | None) -> str:
        if d is None:
            return self.missing
        return self.formatter.format(d)

    def _iter_write(self, *args, **kwargs) -> WriteHandle[tuple[DistanceStatistics]]:
        buffer = None
        with FileHandler.Tabfile(self.path, "w") as file:
            try:
                bunch = yield
                self._write_headers(file, bunch)
                self._write_stats(file, bunch)
                while True:
                    bunch = yield
                    self._write_stats(file, bunch)
            except GeneratorExit:
                if not buffer:
                    return
                self._write_headers(file, buffer)
                self._write_stats(file, buffer)

    def _write_headers(self, file: TextIO, bunch: tuple[DistanceStatistics]):
        raise NotImplementedError()

    def _write_stats(self, file: TextIO, bunch: tuple[DistanceStatistics]):
        raise NotImplementedError()

    def _iter_read(self, *args, **kwargs) -> None:
        raise NotImplementedError()


class SubsetPairsStatisticsHandler(SubsetStatisticsHandler):
    def _write_headers(self, file: TextIO, bunch: tuple[DistanceStatistics]):
        metrics = (str(stats.metric) for stats in bunch)
        combinations = product(metrics, ["mean", "min", "max"])
        headers = (f"{metric} {stat}" for metric, stat in combinations)
        out = ("target", "query", *headers)
        file.write(out)

    def _write_stats(self, file: TextIO, bunch: tuple[DistanceStatistics]):
        idx = bunch[0].idx
        idy = bunch[0].idy
        if idx is None:
            idx = "?"
        if idy is None:
            idy = "?"
        stats = ((stats.mean, stats.min, stats.max) for stats in bunch)
        stats = (self.distanceToText(stat) for stat in chain(*stats))
        out = (idx, idy, *stats)
        file.write(out)


class SubsetIdentityStatisticsHandler(SubsetStatisticsHandler):
    def _write_headers(self, file: TextIO, bunch: tuple[DistanceStatistics]):
        metrics = (str(stats.metric) for stats in bunch)
        combinations = product(metrics, ["mean", "min", "max"])
        headers = (f"{metric} {stat}" for metric, stat in combinations)
        out = ("target", *headers)
        file.write(out)

    def _write_stats(self, file: TextIO, bunch: tuple[DistanceStatistics]):
        idx = bunch[0].idx
        if idx is None:
            idx = "?"
        stats = ((stats.mean, stats.min, stats.max) for stats in bunch)
        stats = (self.distanceToText(stat) for stat in chain(*stats))
        out = (idx, *stats)
        file.write(out)


class SubsetMatrixStatisticsHandler(SubsetStatisticsHandler):
    def _open(
        self,
        path: Path,
        mode: Literal["r", "w"] = "w",
        template: str = "{mean} ({min}-{max})",
        *args,
        **kwargs,
    ):
        self.template = template
        super()._open(path, mode, *args, **kwargs)

    def statsToText(self, stats: DistanceStatistics):
        if not stats.count:
            return self.missing
        mean = self.distanceToText(stats.mean)
        min = self.distanceToText(stats.min)
        max = self.distanceToText(stats.max)
        return self.template.format(mean=mean, min=min, max=max)

    def _iter_write(self) -> WriteHandle[DistanceStatistics]:
        self.buffer: list[DistanceStatistics] = []
        self.wrote_headers = False

        with FileHandler.Tabfile(self.path, "w") as file:
            try:
                line = yield from self._assemble_line()
                self._write_headers(file, line)
                self._write_scores(file, line)
                while True:
                    line = yield from self._assemble_line()
                    self._write_scores(file, line)
            except GeneratorExit:
                line = self.buffer
                self._write_headers(file, line)
                self._write_scores(file, line)
                return

    def _assemble_line(
        self
    ) -> Generator[None, DistanceStatistics, list[DistanceStatistics]]:
        buffer = self.buffer
        try:
            while True:
                distance = yield
                buffer.append(distance)
                if buffer[0].idx != buffer[-1].idx:
                    self.buffer = buffer[-1:]
                    return buffer[:-1]
        except GeneratorExit:
            return

    def _write_headers(self, file: FileHandler.Tabfile, line: list[DistanceStatistics]):
        if self.wrote_headers:
            return
        idys = [distance.idy for distance in line]
        idys = [idy if idy is not None else "?" for idy in idys]
        out = ("", *idys)
        file.write(out)
        self.wrote_headers = True

    def _write_scores(self, file: FileHandler.Tabfile, line: list[DistanceStatistics]):
        scores = [self.statsToText(stats) for stats in line]
        idx = line[0].idx
        if idx is None:
            idx = "?"
        out = (idx, *scores)
        file.write(out)


class SubsetPair(NamedTuple):
    x: str
    y: str


class SubsetDistance(NamedTuple):
    distance: Distance
    genera: SubsetPair | None
    species: SubsetPair | None

    def get_comparison_type(self) -> str:
        same_genera = bool(self.genera.x == self.genera.y) if self.genera else None
        same_species = bool(self.species.x == self.species.y) if self.species else None
        return {
            (None, None): ComparisonType.Unknown,
            (None, True): ComparisonType.IntraSpecies,
            (None, False): ComparisonType.InterSpecies,
            (False, None): ComparisonType.InterGenus,
            (False, True): ComparisonType.InterGenus,
            (False, False): ComparisonType.InterGenus,
            (True, None): ComparisonType.IntraGenus,
            (True, True): ComparisonType.IntraSpecies,
            (True, False): ComparisonType.InterSpecies,
        }[(same_genera, same_species)]


class SummaryHandler(DistanceHandler.Linear.WithExtras):
    def _open(self, path, mode, *args, **kwargs):
        super()._open(path, mode, tagX=" (query 1)", tagY=" (query 2)", *args, **kwargs)

    def _assemble_line(self) -> Generator[None, SubsetDistance, list[SubsetDistance]]:
        buffer = self.buffer
        try:
            while True:
                distance = yield
                buffer.append(distance)
                if any(
                    (
                        buffer[0].distance.x.id != buffer[-1].distance.x.id,
                        buffer[0].distance.y.id != buffer[-1].distance.y.id,
                    )
                ):
                    self.buffer = buffer[-1:]
                    return buffer[:-1]
        except GeneratorExit:
            return

    def _write_headers(self, file: FileHandler.Tabfile, line: list[SubsetDistance]):
        if self.wrote_headers:
            return
        idxHeader = self.idxHeader + self.tagX
        idyHeader = self.idyHeader + self.tagY
        extrasX = [key + self.tagX for key in line[0].distance.x.extras.keys()]
        extrasY = [key + self.tagY for key in line[0].distance.y.extras.keys()]
        metrics = [str(subset_distance.distance.metric) for subset_distance in line]
        infoX = ("genus" + self.tagX, "species" + self.tagX)
        infoY = ("genus" + self.tagY, "species" + self.tagY)
        out = (
            idxHeader,
            idyHeader,
            *metrics,
            *extrasX,
            *extrasY,
            *infoX,
            *infoY,
            "comparison_type",
        )
        file.write(out)
        self.wrote_headers = True

    def _write_scores(self, file: FileHandler.Tabfile, line: list[SubsetDistance]):
        first = line[0]
        idx = first.distance.x.id
        idy = first.distance.y.id
        extrasX = first.distance.x.extras.values()
        extrasY = first.distance.y.extras.values()
        extrasX = [x if x is not None else self.missing for x in extrasX]
        extrasY = [y if y is not None else self.missing for y in extrasY]
        scores = [
            self.distanceToText(subset_distance.distance.d) for subset_distance in line
        ]
        genusX = first.genera.x if first.genera else "-"
        genusY = first.genera.y if first.genera else "-"
        speciesX = first.species.x if first.species else "-"
        speciesY = first.species.y if first.species else "-"
        comparison_type = first.get_comparison_type()
        out = (
            idx,
            idy,
            *scores,
            *extrasX,
            *extrasY,
            genusX or "-",
            speciesX or "-",
            genusY or "-",
            speciesY or "-",
            comparison_type.label,
        )
        file.write(out)


class Results(NamedTuple):
    output_directory: Path
    seconds_taken: float


#     stats_all: Path | None
#     stats_species: Path | None
#     stats_genus: Path | None
#
#     aligned_pairs: Path | None
#
#     summary: Path
#
#     ALOTOFFILES: ...
#
#     matrices: dict[str, Path]
#
#
#     number_of_files_created: int


class VersusAll:
    def __init__(self):
        self.work_dir: Path = None
        self.paths = AttrDict()

        self.progress_handler: Callable = console_report

        self.input = AttrDict()
        self.input.sequences: Sequences = None
        self.input.species: Partition = None
        self.input.genera: Partition = None

        self.params = AttrDict()
        self.progress_interval: float = 0.015

        self.params.pairs = AttrDict()
        self.params.pairs.align: bool = True
        self.params.pairs.write: bool = True
        self.params.pairs.scores: Scores = None

        self.params.distances = AttrDict()
        self.params.distances.metrics: list[DistanceMetric] = None
        self.params.distances.write_linear: bool = True
        self.params.distances.write_matricial: bool = True

        self.params.plot = AttrDict()
        self.params.plot.histograms: bool = True
        self.params.plot.binwidth: float = 0.05
        self.params.plot.formats: list[str] = None
        self.params.plot.palette: list[tuple] = None

        self.params.format = AttrDict()
        self.params.format.float: str = "{:.4f}"
        self.params.format.percentage: str = "{:.2f}"
        self.params.format.missing: str = "NA"
        self.params.format.stats_template: str = "{mean} ({min}-{max})"
        self.params.format.percentage_multiply: bool = False

        self.params.stats = AttrDict()
        self.params.stats.all: bool = True
        self.params.stats.species: bool = True
        self.params.stats.genera: bool = True

    def generate_paths(self):
        assert self.work_dir

        self.paths.summary = self.work_dir / "summary.tsv"
        self.paths.stats_all = self.work_dir / "stats" / "all.tsv"
        self.paths.stats_species = self.work_dir / "stats" / "species.tsv"
        self.paths.stats_genera = self.work_dir / "stats" / "genera.tsv"
        self.paths.aligned_pairs = self.work_dir / "align" / "aligned_pairs.txt"
        self.paths.distances_linear = self.work_dir / "distances" / "linear.tsv"
        self.paths.distances_matricial = self.work_dir / "distances" / "matricial"
        self.paths.subsets = self.work_dir / "subsets"
        self.paths.plots = self.work_dir / "plots"

        for path in [
            self.paths.summary,
        ]:
            self.create_parents(path)

    def create_parents(self, path: Path):
        if path.suffix:
            path = path.parent
        path.mkdir(parents=True, exist_ok=True)

    def check_metrics(self):
        self.params.distances.metrics = self.params.distances.metrics or [
            DistanceMetric.Uncorrected(),
            DistanceMetric.UncorrectedWithGaps(),
            DistanceMetric.JukesCantor(),
            DistanceMetric.Kimura2P(),
        ]

    def calculate_statistics_all(self, sequences: Sequences):
        if not self.params.stats.all:
            yield from sequences
            return

        allStats = StatisticsCalculator()

        for sequence in sequences:
            allStats.add(sequence.seq.upper())
            yield sequence

        self.create_parents(self.paths.stats_all)
        with StatisticsHandler.Single(
            self.paths.stats_all,
            "w",
            float_formatter=self.params.format.float,
            percentage_formatter=self.params.format.percentage,
            percentage_multiply=self.params.format.percentage_multiply,
        ) as file:
            stats = allStats.calculate()
            file.write(stats)

    def calculate_statistics_species(self, sequences: Sequences):
        if not self.input.species:
            return sequences
        if not self.params.stats.species:
            return sequences

        return self._calculate_statistics_partition(
            sequences, self.input.species, "species", self.paths.stats_species
        )

    def calculate_statistics_genera(self, sequences: Sequences):
        if not self.input.genera:
            return sequences
        if not self.params.stats.genera:
            return sequences

        return self._calculate_statistics_partition(
            sequences, self.input.genera, "genera", self.paths.stats_genera
        )

    def _calculate_statistics_partition(
        self, sequences: Sequences, partition: Partition, group_name: str, path: Path
    ):
        try:
            calculators = dict()
            for subset in partition.values():
                if subset not in calculators:
                    calculators[subset] = StatisticsCalculator(group=subset)

            for sequence in sequences:
                subset = partition.get(sequence.id, None)
                if subset is not None:
                    calculators[subset].add(sequence.seq.upper())
                yield sequence

        except GeneratorExit:
            pass

        finally:
            self.create_parents(path)
            with StatisticsHandler.Groups(
                path,
                "w",
                group_name=group_name,
                float_formatter=self.params.format.float,
                percentage_formatter=self.params.format.percentage,
                percentage_multiply=self.params.format.percentage_multiply,
            ) as file:
                for calc in calculators.values():
                    stats = calc.calculate()
                    file.write(stats)

    def normalize_sequences(self, sequences: Sequences) -> Sequences:
        if not self.params.pairs.align:
            return sequences
        return sequences.normalize()

    def align_pairs(self, pairs: SequencePairs):
        if not self.params.pairs.align:
            yield from pairs
            return

        aligner = PairwiseAligner.Biopython(self.params.pairs.scores)
        yield from aligner.align_pairs(pairs)

    def write_pairs(self, pairs: SequencePairs):
        if not self.params.pairs.write:
            yield from pairs
            return

        self.create_parents(self.paths.aligned_pairs)
        with SequencePairHandler.Formatted(self.paths.aligned_pairs, "w") as file:
            for pair in pairs:
                file.write(pair)
                yield pair

    def calculate_distances(self, pairs: SequencePairs):
        for x, y in pairs:
            for metric in self.params.distances.metrics:
                if x != y:
                    yield metric.calculate(x, y)
                else:
                    yield Distance(metric, x, y, None)

    def adjust_distances(self, distances: Distances):
        if not self.params.format.percentage_multiply:
            yield from distances
            return

        for distance in distances:
            if distance.d is not None:
                distance = distance._replace(d=distance.d * 100)
            yield distance

    def write_distances_linear(self, distances: Distances):
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

    def write_distances_multimatrix(self, distances: Distances):
        if not self.params.distances.write_matricial:
            return distances

        self.create_parents(self.paths.distances_matricial)
        for metric in self.params.distances.metrics:
            distances = self._write_distances_matrix(
                distances, metric, self.paths.distances_matricial / f"{str(metric)}.tsv"
            )
        return distances

    def _write_distances_matrix(
        self, distances: Distances, metric: DistanceMetric, path: Path
    ):
        with DistanceHandler.Matrix(
            path,
            "w",
            missing=self.params.format.missing,
            formatter=self.params.format.float,
        ) as file:
            for distance in distances:
                if distance.metric.type == metric.type:
                    file.write(distance)
                yield distance

    def aggregate_distances_species(
        self, distances: Distances
    ) -> Iterator[SubsetPair | None]:
        if not self.input.species:
            return (None for _ in distances)
        return self._aggregate_distances(
            distances, self.input.species, self.paths.subsets / "species"
        )

    def aggregate_distances_genera(
        self, distances: Distances
    ) -> Iterator[SubsetPair | None]:
        if not self.input.genera:
            return (None for _ in distances)
        return self._aggregate_distances(
            distances, self.input.genera, self.paths.subsets / "genera"
        )

    def _aggregate_distances(
        self, distances: Distances, partition: Partition, path: Path
    ) -> Iterator[SubsetPair]:
        try:
            aggregators = dict()
            for metric in self.params.distances.metrics:
                aggregators[str(metric)] = DistanceAggregator(metric)

            for distance in distances:
                subset_x = partition.get(distance.x.id, None)
                subset_y = partition.get(distance.y.id, None)
                generic = GenericDistance(
                    distance.metric, subset_x, subset_y, distance.d
                )
                aggregators[str(generic.metric)].add(generic)
                yield SubsetPair(subset_x, subset_y)

        except GeneratorExit:
            pass

        finally:
            self.write_subset_statistics_linear(aggregators, path / "linear")
            self.write_subset_statistics_matricial(aggregators, path / "matricial")

    def write_subset_statistics_linear(
        self, aggregators: dict[str, DistanceAggregator], path: Path
    ):
        self.create_parents(path)
        with (
            SubsetPairsStatisticsHandler(
                path / "pairs.tsv",
                "w",
                formatter=self.params.format.float,
            ) as pairs_file,
            SubsetIdentityStatisticsHandler(
                path / "identity.tsv",
                "w",
                formatter=self.params.format.float,
            ) as identity_file,
        ):
            aggs = aggregators.values()
            iterators = (iter(agg) for agg in aggs)
            bunches = zip(*iterators)
            for bunch in bunches:
                if bunch[0].idx == bunch[0].idy:
                    identity_file.write(bunch)
                else:
                    pairs_file.write(bunch)

    def write_subset_statistics_matricial(
        self, aggregators: dict[str, DistanceAggregator], path: Path
    ):
        self.create_parents(path)
        for metric, aggregator in aggregators.items():
            with SubsetMatrixStatisticsHandler(
                path / f"{metric}.tsv",
                "w",
                formatter=self.params.format.float,
                template=self.params.format.stats_template,
            ) as file:
                for stats in aggregator:
                    file.write(stats)

    def plot_histograms(self, distances: iter[SubsetDistance]):
        if not self.params.plot.histograms:
            yield from distances

        plotter = HistogramPlotter(
            formats=self.params.plot.formats,
            palette=self.params.plot.palette,
            binwidth=self.params.plot.binwidth,
            binfactor=100.0 if self.params.format.percentage_multiply else 1.0,
        )
        for subset_distance in distances:
            plotter.add(
                str(subset_distance.distance.metric),
                subset_distance.distance.d,
                subset_distance.get_comparison_type(),
            )
            yield subset_distance
        self.create_parents(self.paths.plots)
        plotter.plot(self.paths.plots)

    def write_summary(self, distances: iter[SubsetDistance]):
        with SummaryHandler(
            self.paths.summary,
            "w",
            missing=self.params.format.missing,
            formatter=self.params.format.float,
        ) as file:
            for distance in distances:
                # if 'organism' in distance.x.extras:
                #     del distance.x.extras['organism']
                # if 'organism' in distance.y.extras:
                #     del distance.y.extras['organism']
                file.write(distance)
                yield distance

    def report_progress(self, distances: iter[SubsetDistance]):
        total = len(self.params.distances.metrics) * len(self.input.sequences) ** 2
        last_time = perf_counter()
        for index, distance in enumerate(distances, 1):
            new_time = perf_counter()
            if new_time - last_time >= self.progress_interval:
                self.progress_handler("distance.x.id", index, total)
                last_time = new_time
            yield distance
        self.progress_handler("Finalizing...", total, total)

    def start(self) -> None:
        ts = perf_counter()

        self.generate_paths()
        self.check_metrics()

        sequences = self.input.sequences
        sequences = self.normalize_sequences(sequences)

        sequences_left = sequences
        sequences_left = self.calculate_statistics_all(sequences_left)
        sequences_left = self.calculate_statistics_species(sequences_left)
        sequences_left = self.calculate_statistics_genera(sequences_left)

        pairs = SequencePairs.fromProduct(sequences_left, sequences)
        pairs = self.align_pairs(pairs)
        pairs = self.write_pairs(pairs)

        distances = self.calculate_distances(pairs)
        distances = self.adjust_distances(distances)
        distances = self.write_distances_linear(distances)
        distances = self.write_distances_multimatrix(distances)

        distances = multiply(distances, 3)
        genera_pair = self.aggregate_distances_genera(distances)
        species_pair = self.aggregate_distances_species(distances)
        subset_distances = (
            SubsetDistance(d, g, s)
            for d, g, s in zip(distances, genera_pair, species_pair)
        )

        subset_distances = self.plot_histograms(subset_distances)
        subset_distances = self.write_summary(subset_distances)

        subset_distances = self.report_progress(subset_distances)

        for _ in subset_distances:
            pass

        tf = perf_counter()

        return Results(self.work_dir, tf - ts)
