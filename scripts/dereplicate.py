from pathlib import Path
from sys import argv

from itaxotools.taxi3.tasks.dereplicate import Dereplicate
from itaxotools.taxi3.sequences import Sequences, SequenceHandler


def main(data_path: Path, output_path: Path):
    task = Dereplicate()
    task.work_dir = Path(output_path)
    task.input = Sequences.fromPath(data_path, SequenceHandler.Tabfile, idHeader='seqid', seqHeader='sequence')
    task.params.thresholds.length = 20
    task.params.thresholds.similarity = 0.05
    results = task.start()
    print('')
    print(f'Output directory: {results.output_directory}')
    print(f'Time taken: {results.seconds_taken:.4f}s')


if __name__ == '__main__':
    data_path = Path(argv[1])
    output_path = Path(argv[2])
    main(data_path, output_path)
