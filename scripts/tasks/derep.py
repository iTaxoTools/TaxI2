from itaxotools.taxi3.library.datatypes import CompleteData, ValidFilePath, TabfileReader
from itaxotools.taxi3.library.task import Dereplicate
from pathlib import Path

sequences = CompleteData.from_path(ValidFilePath(Path("input.tab")), TabfileReader)
task = Dereplicate(warn=print)

task.similarity = 0.05
# task.length_threshold = 100
# task.keep_most_complete = False
task.data = sequences
task.start()
for output in task.result:
    output.included.append_to_file(Path("dereplicated.txt"))
    output.excluded.append_to_file(Path("excluded.txt"))
