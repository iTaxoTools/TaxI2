from pathlib import Path
from sys import argv

from itaxotools.taxi2.sequences import SequenceHandler, Sequences
from itaxotools.taxi2.tasks.decontaminate import Decontaminate


def main(data_path: Path, outgroup_path: Path, output_path: Path):
    task = Decontaminate()
    task.work_dir = Path(output_path)
    # task.input = Sequences.fromPath(data_path, SequenceHandler.Fasta)
    task.input = Sequences.fromPath(
        data_path, SequenceHandler.Tabfile, idHeader="seqid", seqHeader="sequence"
    )
    task.outgroup = Sequences.fromPath(
        outgroup_path, SequenceHandler.Tabfile, idHeader="seqid", seqHeader="sequence"
    )
    task.params.thresholds.similarity = 0.07
    task.params.pairs.align = False
    task.set_output_format_from_path(data_path)
    results = task.start()
    print("")
    print(f"Output directory: {results.output_directory}")
    print(f"Time taken: {results.seconds_taken:.4f}s")


if __name__ == "__main__":
    data_path = Path(argv[1])
    outgroup_path = Path(argv[2])
    output_path = Path(argv[3])
    main(data_path, outgroup_path, output_path)
