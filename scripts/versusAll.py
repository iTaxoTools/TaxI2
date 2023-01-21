from pathlib import Path
from sys import argv

from itaxotools.taxi3.tasks.versus_all import VersusAll
from itaxotools.taxi3.sequences import Sequences, SequenceHandler
from itaxotools.taxi3.partitions import Partition, PartitionHandler

def main(input_path: Path, output_path: Path):
    task = VersusAll()
    task.input.sequences = Sequences.fromPath(input_path, SequenceHandler.Tabfile, idHeader='seqid', seqHeader='sequence')
    task.input.species = Partition.fromPath(input_path, PartitionHandler.Tabfile, idHeader='seqid', subHeader='organism')
    task.input.genera = Partition.fromPath(input_path, PartitionHandler.Tabfile, idHeader='seqid', subHeader='organism', filter=PartitionHandler.subset_first_word)
    task.work_dir = Path(output_path)
    results = task.start()
    print('')
    print(f'Output directory: {results.output_directory}')
    print(f'Time taken: {results.seconds_taken:.4f}s')


if __name__ == '__main__':
    input_path = Path(argv[1])
    output_path = Path(argv[2])
    main(input_path, output_path)
