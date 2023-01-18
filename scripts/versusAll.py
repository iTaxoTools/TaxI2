from pathlib import Path
from sys import argv

from itaxotools.taxi3.tasks.versus_all import VersusAll

def main():
    task = VersusAll()
    task.set_input_sequences_from_path(Path(argv[1]))
    task.set_input_species_from_path(Path(argv[1]))
    task.set_input_genera_from_path(Path(argv[1]))
    task.work_dir = Path(argv[2])
    task.start()

if __name__ == '__main__':
    main()
