from pathlib import Path
from sys import argv

from itaxotools.taxi2.distances import DistanceMetric
from itaxotools.taxi2.partitions import Partition, PartitionHandler
from itaxotools.taxi2.sequences import SequenceHandler, Sequences
from itaxotools.taxi2.tasks.versus_all import VersusAll


def main(input_path: Path, output_path: Path):
    task = VersusAll()
    task.work_dir = Path(output_path)
    # task.input.sequences = Sequences.fromPath(input_path, SequenceHandler.Fasta)
    task.input.sequences = Sequences.fromPath(
        input_path, SequenceHandler.Tabfile, idHeader="seqid", seqHeader="sequence"
    )
    task.input.species = Partition.fromPath(
        input_path, PartitionHandler.Tabfile, idHeader="seqid", subHeader="organism"
    )
    task.input.genera = Partition.fromPath(
        input_path,
        PartitionHandler.Tabfile,
        idHeader="seqid",
        subHeader="organism",
        filter=PartitionHandler.subset_first_word,
    )
    task.params.pairs.align = False
    task.params.distances.metrics = [DistanceMetric.NCD()]
    task.params.plot.formats = ["pdf"]
    # task.params.plot.binwidth = 0.02
    # task.params.pairs.align = False
    results = task.start()
    print("")
    print(f"Output directory: {results.output_directory}")
    print(f"Time taken: {results.seconds_taken:.4f}s")


if __name__ == "__main__":
    input_path = Path(argv[1])
    output_path = Path(argv[2])
    main(input_path, output_path)
