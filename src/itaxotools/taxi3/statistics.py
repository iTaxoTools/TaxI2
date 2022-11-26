
# Abstract!

class Sequence_questionmark(NamedTuple):
    idx: str
    idy: str
    seq: str


class SequenceStatistics(NamedTuple):
    subset: str
    numberOfSequence: int
    lessThan100BP: int
    between100To300BP: int
    between301To1000BP: int
    greaterThan1000BP: int
    minimumLenght: int
    maximumLenght: int
    meanLenght: float
    medianLenght: float
    standardDeviation: float
    N50: float
    L50: float
    N90: float
    L90: float
    totalLength: int
    percentageOfA: float
    percentageOfC: float
    percentageOfG: float
    percentageOfT: float
    GCcontent: float
    percentageOfAmbiguity: float
    percentageOfMissingData: float
    percentageOfMissingDataWithGaps: float


class StatisticsCalculator:
    def __init__(self):
        self.state = ...

    def add(self, sequence):
        self.state.append(len(sequence))

    def calculate(self) -> SequenceStatistics:
        pass



allStats = StatisticsCalculator()

statsDict: dict[str, StatisticsCalculator] = dict()
for subset in subsets:
    statsDict[subset] = StatisticsCalculator()


for sequence in sequences:
    allStats.add(sequence)
    subset = get_subset(sequence)  # subset = 'Boophis piperatus'
    statsDict[subset].add(sequence)

results = allStats.calculate()
