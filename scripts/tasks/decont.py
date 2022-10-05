from itaxotools.taxi3.library.datatypes import CompleteData, SequenceData, ValidFilePath, TabfileReader
from itaxotools.taxi3.library.task import Decontaminate, Alignment
from pathlib import Path

data = CompleteData.from_path(ValidFilePath(Path("input.tab")), TabfileReader)
reference = SequenceData.from_path(ValidFilePath(Path("excluded.txt")), TabfileReader)
task = Decontaminate(warn=print)

task.alignment = Alignment.Pairwise
task.data = data
task.reference = reference
task.start()
for output in task.result:
    print(output)
    output.decontaminated.append_to_file(Path("decontaminated.txt"))
    output.contaminates.append_to_file(Path("contaminates.txt"))
    output.summary.append_to_file(Path("summary.txt"))
