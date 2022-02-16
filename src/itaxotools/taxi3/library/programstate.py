import tkinter as tk
import pandas as pd
import numpy as np
import networkx as nx
import os
from pathlib import Path
import sys
import math
import re
import logging
import datetime
import time
import gc
import itertools
from typing import Union, TextIO, Iterator, Tuple, Any, Dict, Optional, List
from .fasta import Fastafile
from .genbank import GenbankFile
from .sequence_statistics import (
    sequence_statistics,
    sequence_statistics_with_gaps,
)
from .alfpy_distance import make_alfpy_distance_table, make_alfpy_distance_table2
from .resources import get_resource


with open(get_resource("options.tab")) as options_file:
    option, _, val = options_file.readline().rstrip().partition("\t")
    if option == "distance_calculation":
        if val not in {"Rust", "Python"}:
            raise ValueError(
                'distance_calculation value in data/options.tab should be either "Rust" or "Python"'
            )
        BACKEND = val
    else:
        raise ValueError("distance_calculation value is missing in data/options.tab")
if BACKEND == "Python":
    from .python_backend import (
        make_distance_table,
        make_distance_table2,
        make_aligner,
        PDISTANCE,
        NDISTANCES,
        distances_short_names,
        show_alignment,
    )
elif BACKEND == "Rust":
    from .rust_backend import (
        make_distance_table,
        make_distance_table2,
        make_aligner,
        PDISTANCE,
        NDISTANCES,
        distances_short_names,
        show_alignment,
    )
else:
    raise RuntimeError("Unexpected BACKEND")


distances_names = [
    "pairwise uncorrected distance",
    "Jukes-Cantor distance",
    "Kimura-2-Parameter distance",
    "pairwise uncorrected distance counting gaps",
]


class FileFormat:
    """
    Interface for file formats supported by the program
    """

    rename_dict = {
        "organism": "species",
        "specimen_identifier": "specimen_voucher",
        "specimen identifier": "specimen_voucher",
        "specimen voucher": "specimen_voucher",
    }

    chunk_size = 100

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        raise NotImplementedError

    @staticmethod
    def rename_columns(name: str) -> str:
        regex = "([a-zA-Z0-9_-]+)"
        name = "".join(re.findall(regex, name))
        if "sequence" in name:
            return "sequence"
        else:
            try:
                return FileFormat.rename_dict[name]
            except KeyError:
                return name

    def load_chunks(
            self, filepath_or_buffer: Union[str, TextIO], chunk_size: int
    ) -> Iterator[pd.DataFrame]:
        raise NotImplementedError


class TabFormat(FileFormat):
    """
    Format for delimiter-separated table
    """

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        with open(filepath_or_buffer, errors="replace") as infile:
            return self.process_table(pd.read_csv(infile, sep="\t", dtype=str))

    @staticmethod
    def process_table(table: pd.DataFrame) -> pd.DataFrame:
        try:
            table = table.rename(columns=str.casefold).rename(
                columns=FileFormat.rename_columns
            )
            if "subspecies" in table.columns:
                return table[
                    ["seqid", "specimen_voucher", "species", "subspecies", "sequence"]
                ].drop_duplicates()
            else:
                return table[
                    ["seqid", "specimen_voucher", "species", "sequence"]
                ].drop_duplicates()
        except KeyError as ex:
            raise ValueError(
                "'seqid', 'specimen_voucher', 'species' or 'organism', or 'sequence' column is missing"
            ) from ex

    def load_chunks(
            self, filepath_or_buffer: Union[str, TextIO], chunk_size: int
    ) -> Iterator[pd.DataFrame]:
        with open(filepath_or_buffer, errors="replace") as infile:
            tables = pd.read_csv(
                infile, sep="\t", dtype=str, chunksize=chunk_size
            )
            for table in tables:
                yield self.process_table(table)


class FastaFormat(FileFormat):
    """
    Format for fasta files
    """

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        if isinstance(filepath_or_buffer, str):
            with open(filepath_or_buffer, errors="replace") as infile:
                return self._load_table(infile)
        else:
            return self._load_table(filepath_or_buffer)

    def _load_table(self, file: TextIO) -> pd.DataFrame:
        _, records = Fastafile.read(file)
        return pd.DataFrame(
            ([record["seqid"], record["sequence"]] for record in records()),
            columns=["seqid", "sequence"],
        ).drop_duplicates(subset="seqid")

    def load_chunks(
            self, filepath_or_buffer: Union[str, TextIO], chunk_size: int
    ) -> Iterator[pd.DataFrame]:
        if isinstance(filepath_or_buffer, str):
            with open(filepath_or_buffer, errors="replace") as infile:
                for chunk in self._load_chunks(infile, chunk_size):
                    yield chunk
        else:
            for chunk in self._load_chunks(filepath_or_buffer, chunk_size):
                yield chunk

    def _load_chunks(self, file: TextIO, chunk_size: int) -> Iterator[pd.DataFrame]:
        _, records_gen = Fastafile.read(file)
        records = records_gen()
        while True:
            chunk = itertools.islice(records, chunk_size)
            table = pd.DataFrame(
                ([record["seqid"], record["sequence"]] for record in chunk),
                columns=["seqid", "sequence"],
            )
            if len(table) == 0:
                break
            yield table.drop_duplicates(subset="seqid").copy()


class GenbankFormat(FileFormat):
    """
    Format for fasta files
    """

    def load_table(self, filepath_or_buffer: Union[str, TextIO]) -> pd.DataFrame:
        if isinstance(filepath_or_buffer, str):
            with open(filepath_or_buffer, errors="replace") as infile:
                return self._load_table(infile)
        else:
            return self._load_table(filepath_or_buffer)

    def _load_table(self, file: TextIO) -> pd.DataFrame:
        _, records = GenbankFile.read(file)
        try:
            return (
                pd.DataFrame((record._fields for record in records()))
                .rename(columns=str.casefold)
                .rename(columns=FileFormat.rename_columns)[
                    ["seqid", "specimen_voucher", "species", "sequence"]
                ]
                .drop_duplicates()
            )
        except KeyError as ex:
            raise ValueError(f"{str(ex)} is missing") from ex


class XLSXFormat(FileFormat):
    """
    Format for xlsx files
    """

    def load_table(self, filepath: str) -> pd.DataFrame:
        return self.process_table(pd.read_excel(filepath, dtype=str, engine="openpyxl"))

    @staticmethod
    def process_table(table: pd.DataFrame) -> pd.DataFrame:
        try:
            table = table.rename(columns=str.casefold).rename(
                columns=FileFormat.rename_columns
            )
            if "subspecies" in table.columns:
                return table[
                    ["seqid", "specimen_voucher", "species", "subspecies", "sequence"]
                ].drop_duplicates()
            else:
                return table[
                    ["seqid", "specimen_voucher", "species", "sequence"]
                ].drop_duplicates()
        except KeyError as ex:
            raise ValueError(
                "'seqid', 'specimen_voucher', 'species' or 'organism', or 'sequence' column is missing"
            ) from ex


class DereplicateSettings:
    """
    Settings for ProgramState.dereplicate
    """

    def __init__(self, root: tk.Misc) -> None:
        self.root = root
        self.similarity = tk.StringVar(root, value="0.07")
        self.length_threshold = tk.StringVar(root)
        self.keep_most_complete = tk.BooleanVar(root, value=False)
        self.save_excluded_replicates = tk.BooleanVar(root, value=False)


class DecontaminateSetting():
    """
    Settings for ProgramState.decontaminate
    """

    def __init__(self, root: tk.Misc) -> None:
        self.root = root
        self.similarity = tk.StringVar(root, value="0.07")


class Decontaminate2Setting():
    """
    Settings for ProgramState.decontaminate2
    """

    def __init__(self, root: tk.Misc) -> None:
        self.root = root
        self.alignment_free = tk.BooleanVar(root, value=True)


class ProgramState:
    """
    Encapsulates the state of TaxI2
    """

    SUMMARY_STATISTICS_NAME = "Summary_statistics.txt"

    COMPARE_REFERENCE = 0
    COMPARE_ALL = 1
    DEREPLICATE = 2
    DECONTAMINATE = 3
    DECONT2 = 4

    formats = dict(
        Tabfile=TabFormat, Fasta=FastaFormat, Genbank=GenbankFormat, XLSX=XLSXFormat
    )

    def __init__(self, root: tk.Misc, output_dir: str) -> None:
        self.root = root
        self.input_format_name = tk.StringVar(root, value="Tabfile")
        self.already_aligned = tk.BooleanVar(root, value=False)
        self.alignment_free = tk.BooleanVar(root, value=False)
        self.distance_options = tuple(
            tk.BooleanVar(root, value=False) for _ in range(NDISTANCES)
        )
        self.distance_options[PDISTANCE].set(True)
        self.mode = tk.IntVar(root, value=0)
        self.print_alignments = tk.BooleanVar(root, value=False)
        self.intra_species_lineages = tk.BooleanVar(root, value=False)
        self.perform_clustering = tk.BooleanVar(root, value=False)
        self.cluster_distance = tk.StringVar(root, value=distances_names[PDISTANCE])
        self.cluster_size = tk.StringVar(root, value="0.05")
        self.dereplicate_settings = DereplicateSettings(root)
        self.decontaminate_settings = DecontaminateSetting(root)
        self.decontaminate2_settings = Decontaminate2Setting(root)
        self.output_dir = output_dir

    def show_progress(self, message: str) -> None:
        full_message = f"{time.monotonic() - self.start_time:.1f}s: {message}\n"
        print(full_message)
        self.root.show_progress(full_message)

    @property
    def input_format(self) -> FileFormat:
        return ProgramState.formats[self.input_format_name.get()]()

    def output_name(self, description: str) -> str:
        return os.path.join(self.output_dir, description.replace(" ", "_") + ".txt")

    def output(self, description: str, table: pd.DataFrame, **kwargs) -> None:
        out_name = self.output_name(description)
        with open(out_name, mode="w") as outfile:
            print(description, file=outfile)
            table.to_csv(
                outfile, sep="\t", line_terminator="\n", float_format="%.4g", **kwargs
            )
        if 'distance' in description:
            table = table.copy()
            convert_percent(table)
            out_name = self.output_name(description + " in percent")
            with open(out_name, mode="w") as outfile:
                print(description, file=outfile)
                table.to_csv(
                    outfile, sep="\t", line_terminator="\n", float_format="%.4g", **kwargs
                )

    def output_square(self, description: str, table: pd.DataFrame, **kwargs) -> None:
        out_name = self.output_name(description)
        table = table.copy()
        table.columns.names = [None] * len(table.columns.names)
        table.index.names = [None] * len(table.index.names)
        with open(out_name, mode="w") as outfile:
            print(description, file=outfile)
            table.to_csv(
                outfile, sep="\t", line_terminator="\n", float_format="%.4g", **kwargs
            )
        if 'distance' in description:
            convert_percent(table)
            out_name = self.output_name(description + " in percent")
            with open(out_name, mode="w") as outfile:
                print(description, file=outfile)
                table.to_csv(
                    outfile, sep="\t", line_terminator="\n", float_format="%.4g", **kwargs
                )

    def process(self, input_file: str) -> None:
        self.start_time = time.monotonic()
        if self.input_format_name.get() == "Genbank" and self.already_aligned.get():
            raise ValueError(
                "'Already aligned' option is not allowed for the Genbank format."
            )
        try:
            table = self.input_format.load_table(input_file)
        except ImportError:
            raise
        self.species_analysis = "species" in table.columns
        self.subspecies_analysis = (
            "subspecies" in table.columns
        ) and self.intra_species_lineages.get()
        table.set_index("seqid", inplace=True)
        if not self.already_aligned.get():
            table["sequence"] = normalize_sequences(table["sequence"])

        if self.print_alignments.get():
            with open(
                os.path.join(self.output_dir, "taxi2_alignments.txt"), "w"
            ) as alignment_file:
                print_alignments(table["sequence"], alignment_file)
            self.show_progress("Alignment is printed")

        self.simple_sequence_statistics(table)

        if not self.alignment_free.get():
            distance_table = make_distance_table(table, self.already_aligned.get())
        else:
            distance_table = make_alfpy_distance_table(table)

        if self.species_analysis:
            species_table = pd.DataFrame(table["species"])
        else:
            species_table = None
        del table
        self.show_progress("Distance calculation")

        # The table of most similar sequences
        self.output(f"Most similar sequences", table_closest(distance_table))
        self.show_progress("Table of closest sequence")

        # The matrix of distances between seqids (input order)
        for kind in (
            kind for kind in range(NDISTANCES) if self.distance_options[kind].get()
        ):
            kind_name = distances_short_names[kind]
            out_table = (
                distance_table[["seqid (query 1)", "seqid (query 2)", kind_name]]
                .set_index(["seqid (query 1)", "seqid (query 2)"])
                .squeeze()
                .unstack()
            )
            out_table.columns.names = [None] * len(out_table.columns.names)
            out_table.index.names = [None] * len(out_table.index.names)
            self.output(
                f"{distances_names[kind]} between sequences",
                out_table,
            )
            del out_table

        self.show_progress("Seqid distance table 1")

        # The matrix of distances between seqids (alphabetical order)
        for kind in (
            kind for kind in range(NDISTANCES) if self.distance_options[kind].get()
        ):
            kind_name = distances_short_names[kind]
            out_table = (
                distance_table[["seqid (query 1)", "seqid (query 2)", kind_name]]
                .set_index(["seqid (query 1)", "seqid (query 2)"])
                .squeeze()
                .unstack()
                .sort_index()
                .sort_index(axis=1)
            )
            out_table.columns.names = [None] * len(out_table.columns.names)
            out_table.index.names = [None] * len(out_table.index.names)
            self.output(
                f"{distances_names[kind]} between sequences (Alphabetical order)",
                out_table,
            )
            del out_table

        self.show_progress("Seqid distance table 2")

        # clustering
        if self.perform_clustering.get():
            self.cluster_analysis(
                distance_table, species_table, os.path.basename(input_file)
            )
            self.show_progress("Cluster analysis")

        if self.species_analysis:
            # The matrix of distances between seqids (order by species)
            for kind in (
                kind for kind in range(NDISTANCES) if self.distance_options[kind].get()
            ):
                kind_name = distances_short_names[kind]
                out_table = (
                    pd.Series(
                        distance_table[kind_name].array,
                        index=pd.MultiIndex.from_frame(
                            distance_table[
                                [
                                    "species (query 1)",
                                    "species (query 2)",
                                    "seqid (query 1)",
                                    "seqid (query 2)",
                                ]
                            ]
                        ),
                    )
                    .squeeze()
                    .unstack(level=["species (query 2)", "seqid (query 2)"])
                    .sort_index()
                    .sort_index(axis=1)
                )
                out_table.columns.names = [None] * len(out_table.columns.names)
                out_table.columns.names = [None] * len(out_table.columns.names)
                self.output(
                    f"{distances_names[kind]} between sequences (Ordered by species)",
                    out_table,
                )
                del out_table
            self.show_progress("Seqid distance table 3")

            genus1 = (
                distance_table["species (query 1)"]
                .str.split(pat=r" |_", n=1, expand=True)
                .iloc[:, 0]
            )
            species1_index = distance_table.columns.get_loc("species (query 1)")
            distance_table.insert(
                loc=species1_index, column="genus (query 1)", value=genus1
            )

            genus2 = (
                distance_table["species (query 2)"]
                .str.split(pat=r" |_", n=1, expand=True)
                .iloc[:, 0]
            )
            species2_index = distance_table.columns.get_loc("species (query 2)")
            distance_table.insert(
                loc=species2_index, column="genus (query 2)", value=genus2
            )

            # The matrices of distances between species
            for kind in (
                kind for kind in range(NDISTANCES) if self.distance_options[kind].get()
            ):
                self.show_progress(distances_names[kind])
                kind_name = distances_short_names[kind]
                species_statistics = (
                    distance_table[
                        [
                            "species (query 1)",
                            "species (query 2)",
                            "genus (query 1)",
                            "genus (query 2)",
                            kind_name,
                        ]
                    ]
                    .groupby(["species (query 1)", "species (query 2)"])
                    .agg(
                        genus_1=pd.NamedAgg(column="genus (query 1)", aggfunc="first"),
                        genus_2=pd.NamedAgg(column="genus (query 2)", aggfunc="first"),
                        min=pd.NamedAgg(column=kind_name, aggfunc="min"),
                        mean=pd.NamedAgg(column=kind_name, aggfunc="mean"),
                        max=pd.NamedAgg(column=kind_name, aggfunc="max"),
                    )
                )
                species_statistics["minmax"] = (
                    species_statistics["min"]
                    .apply(format_float)
                    .str.cat(species_statistics["max"].apply(format_float), sep="-")
                )
                species_statistics["mean_minmax"] = (
                    species_statistics["mean"]
                    .apply(format_float)
                    .str.cat(species_statistics["minmax"] + ")", sep=" (")
                )

                square_species_statistics = species_statistics.unstack(
                    "species (query 2)"
                )

                self.output_square(
                    f"Mean {distances_names[kind]} between species",
                    square_species_statistics["mean"],
                )
                self.show_progress("Mean distance between species")

                self.output_square(
                    f"Minimum and maximum {distances_names[kind]} between species",
                    square_species_statistics["minmax"],
                )
                self.show_progress("Minimum and maximum distance between species")

                self.output_square(
                    f"Mean, minimum and maximum {distances_names[kind]} between species",
                    square_species_statistics["mean_minmax"],
                )
                self.show_progress("Mean, minimum and maximum distance between species")

                del square_species_statistics

                with open(
                    self.output_name(
                        f"Mean, minimum and maximum intra-species {distances_names[kind]}"
                    ),
                    mode="w",
                ) as outfile:
                    print(
                        f"Mean, minimum and maximum intra-species {distances_names[kind]}",
                        file=outfile,
                    )
                    species_statistics.reset_index(inplace=True)
                    species_statistics.loc[
                        species_statistics["species (query 1)"]
                        == species_statistics["species (query 2)"]
                    ].to_csv(
                        outfile,
                        sep="\t",
                        columns=["species (query 1)", "mean_minmax"],
                        header=False,
                        index=False,
                        line_terminator="\n",
                    )
                    outfile.write("\n")
                self.show_progress("Mean, minimum and maximum intra-species distance")

                closest_sequence_indices = (
                    distance_table.loc[
                        distance_table["species (query 1)"]
                        != distance_table["species (query 2)"],
                        ["species (query 1)", kind_name],
                    ]
                    .groupby("species (query 1)")
                    .idxmin()[kind_name]
                    .dropna()
                )
                closest_sequences = distance_table.loc[
                    closest_sequence_indices,
                    ["species (query 1)", kind_name, "seqid (query 2)"],
                ]

                with open(
                    self.output_name(
                        f"Closest sequence from different species with {distances_names[kind]}"
                    ),
                    mode="w",
                ) as outfile:
                    print(
                        f"Closest sequence from different species with {distances_names[kind]}",
                        file=outfile,
                    )
                    print(
                        "species\tdistance (closest sequence of different species)\tseqid (closest sequence of different species)",
                        file=outfile,
                    )
                    closest_sequences.to_csv(
                        outfile,
                        sep="\t",
                        float_format="%.4g",
                        header=False,
                        index=False,
                        line_terminator="\n",
                    )
                    outfile.write("\n")

                del closest_sequences
                del closest_sequence_indices

                self.show_progress("Closest sequences")

                genera_statistics = species_statistics.groupby(
                    ["genus_1", "genus_2"]
                ).agg({"min": "min", "mean": "mean", "max": "max"})
                genera_mean_minmax = (
                    genera_statistics["mean"].map(format_float)
                    + " ("
                    + genera_statistics["min"].map(format_float)
                    + "-"
                    + genera_statistics["max"].map(format_float)
                    + ")"
                )

                self.output_square(
                    f"Mean, minimum and maximum {distances_names[kind]} between genera",
                    genera_mean_minmax.unstack(),
                )
                self.show_progress("Mean, minimum and maximum distances between genera")

                genera_mean_minmax = genera_mean_minmax.reset_index(name="mean_minmax")
                intra_genus = genera_mean_minmax.loc[
                    genera_mean_minmax["genus_1"] == genera_mean_minmax["genus_2"],
                    ["genus_1", "mean_minmax"],
                ]

                with open(
                    self.output_name(
                        f"Mean, minimum and maximum intra-genus {distances_names[kind]}"
                    ),
                    mode="w",
                ) as outfile:
                    print(
                        f"Mean, minimum and maximum intra-genus {distances_names[kind]}",
                        file=outfile,
                    )
                    intra_genus.to_csv(
                        outfile,
                        sep="\t",
                        header=False,
                        index=False,
                        line_terminator="\n",
                    )
                    outfile.write("\n")
                self.show_progress("Mean, minimum and maximum intra-genus distance")

                del species_statistics
                del genera_statistics
                del genera_mean_minmax
                del intra_genus

            same_species = (
                distance_table["species (query 1)"]
                == distance_table["species (query 2)"]
            ).astype("int8")
            same_genus = (
                distance_table["genus (query 1)"] == distance_table["genus (query 2)"]
            ).astype("int8")

            if self.subspecies_analysis:
                same_subspecies = (
                    (
                        distance_table["subspecies (query 1)"]
                        == distance_table["subspecies (query 2)"]
                    )
                    >= same_species
                ).astype("int8")
                comparison_type = same_genus + same_species + same_subspecies + 1
            else:
                comparison_type = same_genus + same_species

            print(comparison_type)

            comparison_type = comparison_type.map(
                {
                    0: "inter-genus",
                    1: "inter-species",
                    2: "intra-species",
                    3: "intra-species (same intraspecific lineage)",
                    4: "intra-species (different intraspecific lineages)",
                }
            )

            # def comparison_type(same_species: bool, same_genus: bool) -> str:
            #     if same_genus:
            #         if same_species:
            #             return 'intra-species'
            #         else:
            #             return 'inter-species'
            #     else:
            #         return 'inter-genus'
            comparison_type_pos = int((len(distance_table.columns) - NDISTANCES) / 2)
            distance_table.insert(
                comparison_type_pos, "comparison_type", comparison_type
            )

            try:
                distance_table["species (query 1)"] = (
                    distance_table["species (query 1)"]
                    .str.split(pat=r" |_", n=1, expand=True)
                    .iloc[:, 1]
                )
            except IndexError as e:
                raise ValueError(
                    "The field organism or species must consist"
                    " of two separate expressions separated by a space or underscore"
                    " corresponding to genus and species name")
            distance_table["species (query 2)"] = (
                distance_table["species (query 2)"]
                .str.split(pat=r" |_", n=1, expand=True)
                .iloc[:, 1]
            )

        self.output("Summary statistics", distance_table, index=False)
        self.show_progress("Final table")

    def dereplicate(self, input_file: str) -> None:
        try:
            tables = self.input_format.load_chunks(input_file, chunk_size=1000)
        except ImportError:
            raise
        try:
            if self.dereplicate_settings.length_threshold.get():
                length_threshold: Optional[int] = int(
                    self.dereplicate_settings.length_threshold.get())
            else:
                length_threshold = 0
        except ValueError:
            logging.warning(
                "Can't parse length threshold for 'dereplicate'."
                " No filtering will be performed")
            length_threshold = None
        if self.dereplicate_settings.similarity.get():
            similarity_threshold = float(self.dereplicate_settings.similarity.get())
        else:
            similarity_threshold = 0.07
        filename, ext = os.path.splitext(input_file)
        dereplicated_file = filename + "_dereplicated" + ext
        excluded_replicates_file = filename + "_excluded_replicates" + ext
        for table in tables:
            self.dereplicate_table(table, length_threshold, similarity_threshold,
                                   dereplicated_file, excluded_replicates_file)

    def dereplicate_table(self, table: pd.DataFrame,
                          length_threshold: Optional[int], similarity_threshold: float,
                          dereplicated_file: str, excluded_replicates_file: str) -> None:
        if length_threshold:
            table = table.loc[table['sequence'].str.len() >= length_threshold]
        table.set_index("seqid", inplace=True)
        distance_table = make_alfpy_distance_table(table[['sequence']].copy())

        # preparing the table
        nodes = distance_table["seqid (query 1)"].unique()
        distance_table = distance_table[
            ["seqid (query 1)", "seqid (query 2)", distances_short_names[0]]
        ].copy()
        distance_table.columns = ["seqid1", "seqid2", "distance"]

        distance_threshold = similarity_threshold

        # calculating components
        connected_table = distance_table.loc[
            (distance_table["distance"] <= distance_threshold)
        ]
        graph = nx.from_pandas_edgelist(
            connected_table, source="seqid1", target="seqid2"
        )
        graph.add_nodes_from(nodes)
        components = nx.connected_components(graph)

        seqids_dereplicated: List[str] = []

        for component in components:
            chosen_seqid = table.loc[component, 'sequence'].str.len().idxmax()
            seqids_dereplicated.append(chosen_seqid)

        with open(dereplicated_file, mode='a', newline='') as outfile:
            table.loc[seqids_dereplicated].to_csv(
                outfile, header=(outfile.tell() == 0), sep='\t')

        table = table.drop(seqids_dereplicated)

        if self.dereplicate_settings.save_excluded_replicates:
            with open(excluded_replicates_file, mode='a', newline='') as outfile:
                table.to_csv(outfile, header=(outfile.tell() == 0), sep='\t')

    def cluster_analysis(
        self,
        distance_table: pd.DataFrame,
        species_table: Optional[pd.DataFrame],
        input_file: str,
    ) -> None:
        with open(self.output_name("Cluster analysis"), mode="w") as output_file:
            # extracting options
            distance_kind = distances_short_names[
                distances_names.index(self.cluster_distance.get())
            ]
            try:
                cluster_threshold = float(self.cluster_size.get())
            except Exception:
                logging.warning(
                    f"Invalid cluster threshold {self.cluster_size.get()}.\nUsing default: 0.3"
                )
                cluster_threshold = 0.3

            # preparing the table
            nodes = distance_table["seqid (query 1)"].unique()
            distance_table = distance_table[
                ["seqid (query 1)", "seqid (query 2)", distance_kind]
            ].copy()
            distance_table.columns = ["seqid1", "seqid2", "distance"]

            # calculating components
            connected_table = distance_table.loc[
                (distance_table["distance"] < cluster_threshold)
            ]
            graph = nx.from_pandas_edgelist(
                connected_table, source="seqid1", target="seqid2"
            )
            graph.add_nodes_from(nodes)
            components = nx.connected_components(graph)

            print(
                f"Samples were clustered at a threshold of {cluster_threshold:.3g} ({cluster_threshold*100:.2g}%) uncorrected p-distance",
                file=output_file,
            )

            # add cluster classification to the table
            cluster_of: Dict[str, int] = {}
            max_samples = 0
            min_samples = len(distance_table)
            for i, component in enumerate(components):
                print(f'Cluster{i+1}: {", ".join(component)}', file=output_file)
                min_samples = min(min_samples, len(component))
                max_samples = max(max_samples, len(component))
                for seqid in component:
                    cluster_of[seqid] = i
            num_clusters = i + 1
            distance_table["cluster1"] = distance_table["seqid1"].map(cluster_of)
            distance_table["cluster2"] = distance_table["seqid2"].map(cluster_of)
            print("\n", file=output_file)

            max_in_cluster_distances = (
                distance_table.loc[
                    distance_table["cluster1"] == distance_table["cluster2"]
                ][["cluster1", "distance"]]
                .groupby("cluster1")
                .max()["distance"]
            )

            print(
                "Maximum intra-sample distance within clusters (marked with # if above specified threshold):",
                file=output_file,
            )
            big_clusters = 0
            for cluster_i, distance in max_in_cluster_distances.items():
                if math.isnan(distance):
                    distance = 0
                if distance > cluster_threshold:
                    big_clusters += 1
                    print(f"Cluster{cluster_i+1}: {distance:.4g} #", file=output_file)
                else:
                    print(f"Cluster{cluster_i+1}: {distance:.4g}", file=output_file)

            output_file.write("\n")

            min_between_cluster_distance = (
                distance_table.loc[
                    distance_table["cluster1"] > distance_table["cluster2"]
                ][["cluster1", "cluster2", "distance"]]
                .groupby(["cluster1", "cluster2"])
                .min()["distance"]
                .unstack()
            )
            for i in range(num_clusters):
                min_between_cluster_distance.at[(i, i)] = 0
            min_between_cluster_distance.sort_index(axis=0, inplace=True)
            min_between_cluster_distance.index.name = ""
            min_between_cluster_distance.rename(
                index=(lambda i: f"Cluster{i+1}"),
                columns=(lambda i: f"Cluster{i+1}"),
                inplace=True,
            )

            print("Minimum distance between clusters:\n", file=output_file)

            min_between_cluster_distance.to_csv(
                output_file,
                sep="\t",
                line_terminator="\n",
                float_format="%.4g",
                index=True,
            )

            output_file.write("\n")
            print("Total number of clusters:", num_clusters, file=output_file)
            print(
                "Number of clusters violating threshold for intra-cluster distances:",
                big_clusters,
                file=output_file,
            )

            output_file.write("\n")
            print(
                f"A total of {num_clusters} clusters were found, containing between {min_samples} and {max_samples} samples.",
                file=output_file,
            )

            output_file.write("\n")

            if species_table is not None:
                species_table["cluster"] = (
                    species_table.index.to_series().map(cluster_of) + 1
                )
                species_table.drop_duplicates(inplace=True)

                print(
                    "Comparison of cluster assignment with species assignment in the input file:",
                    file=output_file,
                )
                exists_multicluster_species = False
                for (species, clusters) in species_table.groupby("species")["cluster"]:
                    if len(clusters) > 1:
                        print(
                            f"Sequences of {species} are included in {len(clusters)} clusters:",
                            ", ".join(clusters.astype(str)),
                            file=output_file,
                        )
                        exists_multicluster_species = True
                if exists_multicluster_species:
                    print(
                        "Sequences of all other species are included in only one cluster, respectively",
                        file=output_file,
                    )
                else:
                    print(
                        "Sequences of all species are included in only one cluster, respectively.",
                        file=output_file,
                    )

                output_file.write("\n")

                print(
                    "List of clusters containing sequences of more than one species (according to species assignment in the input file):",
                    file=output_file,
                )
                exists_multispecies_cluster = False
                for (cluster, species) in species_table.groupby("cluster")["species"]:
                    if len(species) > 1:
                        print(
                            f"Cluster {cluster} contains sequences of {len(species)} species:",
                            ", ".join(species),
                            file=output_file,
                        )
                        exists_multispecies_cluster = True
                if exists_multispecies_cluster:
                    print(
                        "All other clusters contain sequences of only a single species, respectively.",
                        file=output_file,
                    )
                else:
                    print(
                        "All clusters contain sequences of only one species, respectively.",
                        file=output_file,
                    )

                output_file.write("\n")

        time = datetime.datetime.now()
        with open(
            os.path.join(
                self.output_dir,
                "taxi2_cluster_" + time.strftime("%Y-%m-%dT%H%M%S") + ".spart",
            ),
            mode="w",
        ) as spart_file:
            print("begin spart;", file=spart_file)
            print("project_name = taxi2_clustering;", file=spart_file)
            print("Date = " + time.astimezone().isoformat(), file=spart_file)
            print(
                "N_spartitions = "
                + str(num_clusters)
                + " : "
                + spart_form(input_file)
                + ";",
                file=spart_file,
            )
            print("N_individuals = " + str(len(cluster_of)) + ";", file=spart_file)
            print("N_subsets = " + str(num_clusters) + ";", file=spart_file)
            print(
                f"[Generated by a simple {cluster_threshold*100:.2g}% threshold clustering in taxi2]",
                file=spart_file,
            )
            print(
                "[WARNING: The sample names below may have been changed to fit SPART specification (only alphanumeric characters and _ )]",
                file=spart_file,
            )
            print(
                f"[The following clusters included sequences differing by a distance above the threshold: {cluster_threshold:.3g}]",
                file=spart_file,
            )
            print("Individual_assignment = ", file=spart_file, end="")
            for specimen, cluster in cluster_of.items():
                print(
                    f"\n{spart_form(specimen)}: {cluster + 1}", file=spart_file, end=""
                )
            print(";\n", file=spart_file)
            print("end;", file=spart_file)

    def reference_comparison_process(
        self, input_file: str, reference_file: str
    ) -> None:
        self.start_time = time.monotonic()
        if self.input_format_name.get() in {"Genbank", "XLSX"}:
            raise ValueError(
                f"Comparison with reference database is not implemented for the {self.input_format_name.get()} format"
            )
        reference_table = self.input_format.load_table(reference_file)
        reference_table.set_index("seqid", inplace=True)
        if not self.already_aligned.get():
            reference_table["sequence"] = normalize_sequences(
                reference_table["sequence"]
            )

        with open(self.output_name("Closest reference sequences"), mode="w") as outfile:
            header = True
            sequences_num = 0
            for table in self.input_format.load_chunks(input_file, chunk_size=100):
                table.set_index("seqid", inplace=True)
                if not self.already_aligned.get():
                    table["sequence"] = normalize_sequences(table["sequence"])
                if not self.alignment_free.get():
                    distance_table = make_distance_table2(
                        table, reference_table, self.already_aligned.get()
                    )
                else:
                    distance_table = make_alfpy_distance_table2(table, reference_table)
                pdistance_name = distances_short_names[PDISTANCE]
                indices_closest = (
                    distance_table[["seqid (query 1)", pdistance_name]]
                    .groupby("seqid (query 1)")
                    .idxmin()[pdistance_name]
                    .squeeze()
                    .dropna()
                )
                closest_table = distance_table.loc[indices_closest].rename(
                    columns=(
                        lambda col: col.replace("query 2", "closest reference sequence")
                    )
                )
                closest_table = distance_table.loc[indices_closest].rename(
                    columns=(lambda col: col.replace("query 1", "query"))
                )
                closest_table.to_csv(
                    outfile,
                    sep="\t",
                    line_terminator="\n",
                    float_format="%.4g",
                    header=header,
                )
                header = False
                sequences_num += self.input_format.chunk_size
                outfile.flush()
                del closest_table
                del distance_table
                del table
                gc.collect()
                self.show_progress(f"{sequences_num} sequences processed")

    def decontaminate(
        self, input_file: str, reference_file: str
    ) -> None:
        self.start_time = time.monotonic()
        if self.input_format_name.get() in {"Genbank", "XLSX"}:
            raise ValueError(
                f"Decontamination is not implemented for the {self.input_format_name.get()} format"
            )
        reference_table = self.input_format.load_table(reference_file)
        reference_table.set_index("seqid", inplace=True)
        if not self.already_aligned.get():
            reference_table["sequence"] = normalize_sequences(
                reference_table["sequence"]
            )

        if self.decontaminate_settings.similarity.get():
            similarity_threshold = float(self.decontaminate_settings.similarity.get())
        else:
            similarity_threshold = 0.07

        filename, ext = os.path.splitext(input_file)
        decontaminated_file = filename + "_decontaminated" + ext
        contaminates_file = filename + "_contaminates" + ext
        summary_file = filename + "_decontaminate_summary" + ext

        header = True
        sequences_num = 0
        for table in self.input_format.load_chunks(input_file, chunk_size=1000):
            table.set_index("seqid", inplace=True)
            if not self.already_aligned.get():
                table["sequence"] = normalize_sequences(table["sequence"])
            distance_table = make_alfpy_distance_table2(table.copy(), reference_table)
            pdistance_name = distances_short_names[PDISTANCE]
            indices_closest = (
                distance_table[["seqid (query 1)", pdistance_name]]
                .groupby("seqid (query 1)")
                .idxmin()[pdistance_name]
                .squeeze()
                .dropna()
            )
            closest_table = distance_table.loc[indices_closest].rename(
                columns=(
                    lambda col: col.replace(
                        "query 2", "closest possible contaminant")
                )
            )
            closest_table = closest_table.rename(
                columns=(lambda col: col.replace("query 1", "query"))
            )
            closest_table["is_contaminant"] = ""
            closest_table.loc[closest_table[pdistance_name] <=
                              similarity_threshold, "is_contaminant"] = "contaminant"
            decontaminates_seqids = closest_table.loc[closest_table[pdistance_name]
                                                      > similarity_threshold, "seqid (query)"]
            with open(summary_file, mode="a") as outfile:
                closest_table.to_csv(
                    outfile,
                    sep="\t",
                    line_terminator="\n",
                    float_format="%.4g",
                    header=header,
                )
            with open(decontaminated_file, mode='a', newline='') as outfile:
                table.loc[decontaminates_seqids].to_csv(
                    outfile, header=(outfile.tell() == 0), sep='\t')

            table = table.drop(decontaminates_seqids)

            with open(contaminates_file, mode='a', newline='') as outfile:
                table.to_csv(outfile, header=(outfile.tell() == 0), sep='\t')

            header = False
            sequences_num += 1000
            del closest_table
            del distance_table
            del table
            gc.collect()
            self.show_progress(f"{sequences_num} sequences processed")

    def decontaminate2(
        self, input_file: str, outgroup_reference_file: str,
        ingroup_reference_file: str
    ) -> None:
        self.start_time = time.monotonic()
        if self.input_format_name.get() in {"Genbank", "XLSX"}:
            raise ValueError(
                f"Decontamination is not implemented for the {self.input_format_name.get()} format"
            )
        ingroup_table = self.input_format.load_table(ingroup_reference_file)
        ingroup_table.set_index("seqid", inplace=True)
        outgroup_table = self.input_format.load_table(outgroup_reference_file)
        outgroup_table.set_index("seqid", inplace=True)
        if not self.already_aligned.get():
            ingroup_table["sequence"] = normalize_sequences(
                ingroup_table["sequence"]
            )
            outgroup_table["sequence"] = normalize_sequences(
                outgroup_table["sequence"]
            )

        filename, ext = os.path.splitext(input_file)
        included_seq_file = filename + "_included_sequences" + ext
        excluded_seq_file = filename + "_excluded_sequences" + ext

        header = True
        sequences_num = 0
        for table in self.input_format.load_chunks(input_file, chunk_size=1000):
            table.set_index("seqid", inplace=True)
            if not self.already_aligned.get():
                table["sequence"] = normalize_sequences(table["sequence"])
            if self.decontaminate2_settings.alignment_free.get():
                ingroup_distance_table = make_alfpy_distance_table2(
                    table.copy(), ingroup_table)
                outgroup_distance_table = make_alfpy_distance_table2(
                    table.copy(), outgroup_table)
            else:
                ingroup_distance_table = make_distance_table2(
                    table.copy(), ingroup_table, self.already_aligned.get())
                outgroup_distance_table = make_distance_table2(
                    table.copy(), outgroup_table, self.already_aligned.get())
            pdistance_name = distances_short_names[PDISTANCE]
            ingroup_closest_distance = (
                ingroup_distance_table[["seqid (query 1)", pdistance_name]]
                .groupby("seqid (query 1)")
                .min()
                .squeeze()
                .dropna()
            )
            outgroup_closest_distance = (
                outgroup_distance_table[["seqid (query 1)", pdistance_name]]
                .groupby("seqid (query 1)")
                .min()
                .squeeze()
                .dropna()
            )
            included_sequences = table.loc[
                ingroup_closest_distance < outgroup_closest_distance]
            excluded_sequences = table.loc[
                ingroup_closest_distance >= outgroup_closest_distance]
            with open(included_seq_file, mode="a") as outfile:
                included_sequences.to_csv(
                    outfile,
                    sep="\t",
                    line_terminator="\n",
                    float_format="%.4g",
                    header=header,
                )

            with open(excluded_seq_file, mode="a") as outfile:
                excluded_sequences.to_csv(
                    outfile,
                    sep="\t",
                    line_terminator="\n",
                    float_format="%.4g",
                    header=header,
                )

            header = False
            sequences_num += 1000
            del included_sequences
            del excluded_sequences
            del ingroup_closest_distance
            del outgroup_closest_distance
            del ingroup_distance_table
            del outgroup_distance_table
            del table
            gc.collect()
            self.show_progress(f"{sequences_num} sequences processed")

    def simple_sequence_statistics(self, table: pd.DataFrame) -> None:
        table_no_dashes = table.copy()
        table_no_dashes["sequence"] = table_no_dashes["sequence"].str.replace("-", "")

        total_stats = table_no_dashes["sequence"].agg(sequence_statistics)
        total_stats_with_gaps = table["sequence"].agg(sequence_statistics_with_gaps)
        with open(
            os.path.join(self.output_dir, "Sequence summary statistics.txt"), mode="w"
        ) as outfile:
            for stat_name, stat_value in total_stats.items():
                print(stat_name, f"{stat_value:.4g}", sep=": ", file=outfile)
            for stat_name, stat_value in total_stats_with_gaps.items():
                print(stat_name, f"{stat_value:.4g}", sep=": ", file=outfile)

            print(file=outfile)
            if self.species_analysis:
                species_stats = table_no_dashes.groupby("species")["sequence"].agg(
                    list(sequence_statistics.values())
                )
                species_stats.columns = list(sequence_statistics)
                species_stats[list(sequence_statistics_with_gaps)] = table.groupby(
                    "species"
                )["sequence"].agg(list(sequence_statistics_with_gaps.values()))
                species_stats.to_csv(
                    outfile, sep="\t", line_terminator="\n", float_format="%.4g"
                )


def format_float(x: float) -> str:
    return f"{x:.4g}"


def spart_form(s: str) -> str:
    return "_".join(re.compile("[A-Za-z0-9_]+").findall(s))


def select_distance(distance_table: pd.DataFrame, kind: int) -> pd.DataFrame:
    """
    Transforms table of arrays of distances into the table of distances

    kind should be one of the distance selection constants
    """
    if kind >= NDISTANCES:
        raise ValueError(f"{kind} doesn't corresponds to a valid distance")
    return distance_table.applymap(lambda arr: arr[kind])


def seqid_distance_table(distance_table: pd.DataFrame) -> pd.DataFrame:
    """
    Changes the index of the table to only seqid
    """
    result = distance_table.copy()
    result.index = distance_table.index.get_level_values("seqid")
    result.columns = distance_table.columns.get_level_values("seqid")
    return result


def species_distance_table(distance_table: pd.DataFrame) -> pd.DataFrame:
    """
    Changes the index of the table to seqid and species
    """
    result = distance_table.copy()
    to_drop = [
        level for level in result.index.names if level not in {"seqid", "species"}
    ]
    result.index = result.index.droplevel(to_drop)
    result.columns = result.columns.droplevel(to_drop)
    return result


def normalize_sequences(sequences: pd.Series) -> pd.Series:
    return sequences.str.upper().str.replace("?", "N").str.replace("-", "")


def alignment_file_name(output_file: str) -> str:
    output_file_base, output_file_ext = os.path.splitext(output_file)
    return output_file_base + "_alignments" + output_file_ext


def print_alignments(sequences: pd.Series, alignment_file: TextIO) -> None:
    sequences = sequences.copy()
    sequences.index = sequences.index.get_level_values("seqid")
    aligner = make_aligner()
    for (seqid_target, target) in sequences.items():
        for (seqid_query, query) in sequences.items():
            print(f"{seqid_target} <-> {seqid_query}", file=alignment_file)
            show_alignment(aligner, target, query, alignment_file)


def table_closest(distance_table: pd.DataFrame) -> pd.DataFrame:
    pdistance_name = distances_short_names[PDISTANCE]
    indices_closest = (
        distance_table[["seqid (query 1)", pdistance_name]]
        .groupby("seqid (query 1)")
        .idxmin()[pdistance_name]
        .squeeze()
        .dropna()
    )
    closest_table = distance_table.loc[indices_closest].rename(
        columns=(
            lambda col: col.replace("query 2", "most similar sequence in the dataset")
        )
    )
    return closest_table


def series_append(series_tuple: pd.Series, series_elem: pd.Series) -> pd.Series:
    """
    Combines a Series of tuples with a Series of scalars by appending componentwise
    """
    return series_tuple.combine(series_elem, lambda tuple, elem: tuple + (elem,))


def show_min_max(mean_min_max: Tuple[float, float, float]) -> str:
    if np.isnan(mean_min_max[0]):
        return ""
    else:
        return f"{mean_min_max[1]:.4g}-{mean_min_max[2]:.4g}"


def show_mean_min_max(mean_min_max: Tuple[float, float, float]) -> str:
    if np.isnan(mean_min_max[0]):
        return ""
    else:
        return f"{mean_min_max[0]:.4g} ({mean_min_max[1]:.4g}-{mean_min_max[2]:.4g})"


def find_closest_from_another(table: pd.DataFrame) -> pd.Series:
    """
    Returns as Series of tuples:
    species | (seqid_of_closest, species_of_closest, seqid_of_self)
    """
    table = table.copy()
    for lbl in table.index.levels[1]:
        table.loc[(slice(None), lbl), (slice(None), lbl)] = np.nan
    return table.stack(level=0).idxmin()


def select_genus(species: str) -> str:
    genus, _ = re.split(r"[ _]", species, maxsplit=1)
    return genus


def convert_percent(table: pd.DataFrame) -> None:
    """
    Multiplies all numbers in table by 100
    """
    table[table.select_dtypes(include='number').columns] *= 100
    num_regex = re.compile(r'(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?')
    for column in table.columns:
        try:
            table[column] = table[column].str.replace(
                num_regex, lambda m: f"{(100 * float(m.group())):.4g}", regex=True)
        except AttributeError:
            pass
