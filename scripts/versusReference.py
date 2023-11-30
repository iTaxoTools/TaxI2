from pathlib import Path
from sys import argv

from itaxotools.taxi2.sequences import SequenceHandler, Sequences
from itaxotools.taxi2.tasks.versus_reference import VersusReference


def main(data_path: Path, reference_path: Path, output_path: Path):
    task = VersusReference()
    task.work_dir = Path(output_path)
    task.input.data = Sequences.fromPath(
        data_path, SequenceHandler.Tabfile, idHeader="seqid", seqHeader="sequence"
    )
    task.input.reference = Sequences.fromPath(
        reference_path, SequenceHandler.Tabfile, idHeader="seqid", seqHeader="sequence"
    )
    task.params.pairs.align = False
    results = task.start()
    print("")
    print(f"Output directory: {results.output_directory}")
    print(f"Time taken: {results.seconds_taken:.4f}s")


if __name__ == "__main__":
    data_path = Path(argv[1])
    reference_path = Path(argv[2])
    output_path = Path(argv[3])
    main(data_path, reference_path, output_path)
