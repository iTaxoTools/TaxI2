from itertools import groupby
from pathlib import Path
from sys import argv
from time import perf_counter
import os
from itaxotools.taxi3.align import *
from itaxotools.taxi3.distances import *
from itaxotools.taxi3.pairs import *
from itaxotools.taxi3.sequences import *


def calc(aligned_pairs, metric=DistanceMetric.Uncorrected()):
    for x, y in aligned_pairs:
        yield metric.calculate(x, y)


def progress(distances, total):

    for index, distance in enumerate(distances):
        index+=1
        print(f"\r Loading... {int(index)}/{total} = {(round(float(index)/float(total)  * 100, 2))} ", end="")
        yield distance


def decontaminateSeq(allpairs, contaminatePath, decontaminatePath, summaryPath, headers, outgroup_weight=1.0):

    with open(decontaminatePath, 'w') as decontaminatedFile, open(contaminatePath, 'w') as contaminatedFile, open(summaryPath, 'w') as summryFile:
        summryFile.write("seqid_query\tclosest possible contaminant\tdistance\tis_contaminant")
        decontaminatedFile.write('\t'.join(headers))
        contaminatedFile.write('\t'.join(headers))

        for ingroup_pair, ingroup_disctance, outgroup_pair, outgroup_disctance in allpairs:
            isContaminate = ingroup_disctance.d > (outgroup_weight * outgroup_disctance.d)
            extras = [v for k, v in ingroup_pair.x.extras.items()]
            extrasString = '\t'.join(extras)
            print(f'\nthis is distance ingroup: {ingroup_disctance.x.id} {ingroup_disctance.y.id}')
            print(f'this is distance outgroup: {outgroup_disctance.x.id} {outgroup_disctance.y.id}\n')

            if isContaminate:
                decontaminatedFile.write(f"\n{ingroup_disctance.x.id}\t{extrasString}\t{ingroup_pair.x.seq}")
            else:
                contaminatedFile.write(f"\n{ingroup_disctance.x.id}\t{extrasString}\t{ingroup_pair.x.seq}")

            summryFile.write(f'\n{ingroup_disctance.x.id}\t{ingroup_disctance.y.id}\t{ingroup_disctance.d}\t{isContaminate}')

            yield ingroup_disctance


def get_minimum(distances):
    for k, g in groupby(distances, lambda d: d.x.id):
        g = (d for d in g if d.d is not None)
        d = min(g, key = lambda x: x.d)
        yield d


def multiply(g, n):
    return (x for x in g for i in range(n))


def main():

    path_data = Path(argv[1])
    path_ingroup_reference = Path(argv[2])
    path_outgroup_reference = Path(argv[3])
    contaminatePath = 'contaminates.txt'
    decontaminatePath = 'decontaminated.txt'
    summaryPath = 'summary.txt'

    ts = perf_counter()

    file_data = SequenceFile.Tabfile(path_data)
    file_ingroup_reference = SequenceFile.Tabfile(path_ingroup_reference)
    file_outgroup_reference = SequenceFile.Tabfile(path_outgroup_reference)

    data = Sequences.fromFile(file_data, idHeader='seqid', seqHeader='sequence')
    ingroup_reference = Sequences.fromFile(file_ingroup_reference, idHeader='seqid', seqHeader='sequence')
    outgroup_reference = Sequences.fromFile(file_outgroup_reference, idHeader='seqid', seqHeader='sequence')

    total = len(data)

    data = data.normalize()

    ingroup_refe = ingroup_reference.normalize()
    outgroup_refe = outgroup_reference.normalize()

    in_pairs = SequencePairs.fromProduct(data, ingroup_refe)
    out_pairs = SequencePairs.fromProduct(data, outgroup_refe)

    #Ingroup
    aligner = PairwiseAligner.Biopython()
    aligned_in_pairs = aligner.align_pairs(in_pairs)

    in_distances = calc(aligned_in_pairs)

    in_minimums = get_minimum(in_distances)

    #outGroup
    aligned_out_pairs = aligner.align_pairs(out_pairs)

    out_distances = calc(aligned_out_pairs)

    out_minimums = get_minimum(out_distances)

    allPairs = zip(in_pairs, in_minimums, out_pairs, out_minimums)

    headers = file_data.getHeader()

    d = decontaminateSeq(allPairs, contaminatePath, decontaminatePath, summaryPath, headers)

    distance = progress(d, total)

    for _ in distance:
        pass

    tf = perf_counter()

    print(f'Time taken: {tf-ts:.4f}s')


if __name__ == '__main__':
    main()
