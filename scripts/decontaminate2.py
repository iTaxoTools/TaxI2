from pathlib import Path
from sys import argv

from itaxotools.taxi2.sequences import SequenceHandler, Sequences
from itaxotools.taxi2.tasks.decontaminate2 import Decontaminate2


def main(data_path: Path, outgroup_path: Path, ingroup_path: Path, output_path: Path):
    task = Decontaminate2()
    task.work_dir = Path(output_path)
    # task.input = Sequences.fromPath(data_path, SequenceHandler.Fasta)
    task.input = Sequences.fromPath(
        data_path, SequenceHandler.Tabfile, idHeader="seqid", seqHeader="sequence"
    )
    task.outgroup = Sequences.fromPath(
        outgroup_path, SequenceHandler.Tabfile, idHeader="seqid", seqHeader="sequence"
    )
    task.ingroup = Sequences.fromPath(
        ingroup_path, SequenceHandler.Tabfile, idHeader="seqid", seqHeader="sequence"
    )
    task.params.weights.outgroup = 1
    task.params.weights.ingroup = 1
    task.params.pairs.align = False
    task.set_output_format_from_path(data_path)
    results = task.start()
    print("")
    print(f"Output directory: {results.output_directory}")
    print(f"Time taken: {results.seconds_taken:.4f}s")


if __name__ == "__main__":
    data_path = Path(argv[1])
    outgroup_path = Path(argv[2])
    ingroup_path = Path(argv[3])
    output_path = Path(argv[4])
    main(data_path, outgroup_path, ingroup_path, output_path)
