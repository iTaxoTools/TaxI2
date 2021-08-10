import os
import sys
import shutil
from typing import Any, Callable, Iterator
import logging

import tkinter as tk
import tkinter.ttk as ttk
import tkinter.font as tkfont
import tkinter.messagebox as tkmessagebox
import tkinter.filedialog as tkfiledialog

from library.programstate import *
from library.gui_utils import display_errors_and_warnings
from library.plot_taxi import Plot

resource_path = getattr(sys, "_MEIPASS", sys.path[0])


class TkLogger(logging.Handler):
    """docstring for TkWarnLogger."""

    def __init__(self, level=logging.NOTSET):
        logging.Handler.__init__(self, level)
        self.addFilter(
            lambda record: record.levelno == logging.WARNING
            or record.levelno == logging.ERROR
        )

    def emit(self, record: logging.LogRecord) -> None:
        if record.levelno == logging.WARNING:
            tkmessagebox.showwarning("Warning", record.getMessage())
        else:
            tkmessagebox.showerror("Error", record.getMessage())
        print(record.pathname, record.lineno, sep=": ")
        print(record.stack_info)
        if record.levelno == logging.WARNING:
            print("Warning:", record.getMessage(), "\n")
        else:
            print("Error:", record.getMessage(), "\n")


class TaxiGUI(ttk.Frame):
    def __init__(self, *args: Any, preview_dir, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)

        self.images = {}
        self.load_images(
            {
                "txt_icon": "file-text.png",
                "graph_icon": "file-graph.png",
                "log_icon": "file-log.png",
                "open_button": "open.png",
                "save_button": "save.png",
                "save_all_button": "save_all.png",
                "run_button": "run.png",
                "clear_button": "clear.png",
            }
        )
        self.preview_dir = preview_dir
        self.programstate = ProgramState(self, self.preview_dir)

        # make directory for graph previews
        os.mkdir(os.path.join(self.preview_dir, "graph_previews"))

        self.panes = ttk.Panedwindow(self, orient="horizontal")
        self.panes.grid(row=8, column=0, sticky="nsew")
        self.create_top_frame()
        self.create_parameters_frame()
        self.create_filelist_frame()
        self.create_preview_frame()

        ttk.Separator(self, orient="horizontal").grid(row=1, column=0, sticky="we")

        self.input_file = tk.StringVar()
        ttk.Label(self, text="Input file").grid(row=2, column=0, sticky="w")
        ttk.Entry(self, textvariable=self.input_file).grid(row=3, column=0, sticky="we")

        # spacing
        ttk.Label(self, font=tkfont.Font(size=5)).grid(row=4, column=0)

        self.reference_file = tk.StringVar()
        ttk.Label(self, text="Reference").grid(row=5, column=0, sticky="w")
        ttk.Entry(self, textvariable=self.reference_file).grid(
            row=6, column=0, sticky="we"
        )

        ttk.Label(self, font=tkfont.Font(size=5)).grid(row=7, column=0)

        self.rowconfigure(8, weight=1)
        self.columnconfigure(0, weight=1)
        self.grid(row=0, column=0, sticky="nsew")
        logger = logging.getLogger()
        logger.addHandler(TkLogger())
        logger.setLevel(logging.WARNING)

    def load_images(self, image_dict: Dict[str, str]) -> None:
        for key, file in image_dict.items():
            self.images[key] = tk.PhotoImage(
                file=os.path.join(resource_path, "data", file)
            )

    def create_top_frame(self) -> None:
        top_frame = ttk.Frame(self, relief="sunken", padding=4)
        top_frame.rowconfigure(0, weight=1)
        top_frame.grid(row=0, column=0, sticky="nsew")

        current_column = 0

        ttk.Label(top_frame, text="TaxI3", font=tkfont.Font(size=20), padding=5).grid(
            row=0, column=current_column
        )
        current_column += 1
        ttk.Label(top_frame, text="Taxonomic identifications\nfrom DNA barcodes").grid(
            row=0, column=current_column
        )
        current_column += 1
        ttk.Separator(top_frame, orient="vertical").grid(
            row=0, column=current_column, sticky="nsew")
        current_column += 1

        ttk.Radiobutton(
            top_frame,
            text="Compare sequences\nagainst reference\ndatabase",
            variable=self.programstate.mode,
            value=ProgramState.COMPARE_REFERENCE,
        ).grid(row=0, column=current_column, sticky="nsew")
        current_column += 1
        ttk.Radiobutton(
            top_frame,
            text="All-against-all\nsequence comparison\nwith genetic distance\nanalysis and clustering",
            variable=self.programstate.mode,
            value=ProgramState.COMPARE_ALL,
        ).grid(row=0, column=current_column, sticky="nsew")
        current_column += 1
        ttk.Radiobutton(
            top_frame,
            text="DEREP",
            variable=self.programstate.mode,
            value=ProgramState.DEREPLICATE,
        ).grid(row=0, column=current_column, sticky="nsew")
        current_column += 1
        ttk.Radiobutton(
            top_frame,
            text="DECONT",
            variable=self.programstate.mode,
            value=ProgramState.DECONTAMINATE,
        ).grid(row=0, column=current_column, sticky="nsew")
        current_column += 1

        for image_key, text, command in (
            (
                "open_button",
                "open reference\nsequence database",
                self.open_reference_command,
            ),
            ("open_button", "open input file\n(query sequences)", self.open_command),
            ("save_button", "save", self.save_command("selected")),
            ("save_all_button", "save_all", self.save_command("all")),
            ("run_button", "run", self.run_command),
            ("clear_button", "clear", self.clear_command),
        ):
            ttk.Button(
                top_frame,
                text=text,
                image=self.images[image_key],
                compound="top",
                style="Toolbutton",
                padding=(10, 0),
                command=command,
            ).grid(row=0, column=current_column, sticky="w")
            current_column += 1

        ttk.Separator(top_frame, orient="vertical").grid(
            row=0, column=current_column, sticky="nsew"
        )
        current_column += 1
        self.images["logo"] = tk.PhotoImage(
            file=os.path.join(
                resource_path, "data", "iTaxoTools Digital linneaeus MICROLOGO.png"
            )
        )
        ttk.Label(top_frame, image=self.images["logo"]).grid(
            row=0, column=current_column, sticky="nse"
        )
        current_column += 1

        top_frame.columnconfigure(current_column - 3, weight=1)

    def open_command(self) -> None:
        path = tkfiledialog.askopenfilename()
        if not path:
            return
        self.input_file.set(os.path.abspath(path))

    def open_reference_command(self) -> None:
        path = tkfiledialog.askopenfilename()
        if not path:
            return
        self.reference_file.set(os.path.abspath(path))

    def save_command(self, which: str) -> Callable[[], None]:
        """
        which should be either "all" or "selected"
        """

        def command():
            save_folder = tkfiledialog.askdirectory()
            if not save_folder:
                return
            for file in self.outfilenames(which):
                full_filename = os.path.join(self.preview_dir, file)
                shutil.copy(full_filename, save_folder)

        return command

    def show_progress(self, message: str) -> None:
        """
        Adds message to preview textbox
        """
        self.preview.insert("end", message)
        self.update()

    def only_pairwise_distance(
        self,
    ) -> bool:
        chosen_distances = list(
            map(tk.BooleanVar.get, self.programstate.distance_options)
        )
        return chosen_distances[0] and chosen_distances.count(True) == 1

    def show_incompatible_options_error(
        self,
    ) -> bool:
        errors = [
            (
                self.programstate.alignment_free.get()
                and not self.only_pairwise_distance(),
                "Only pairwise uncorrelated distance is supported for alignment-free distance calculation",
            ),
            (
                self.programstate.alignment_free.get()
                and self.programstate.print_alignments.get(),
                "The Print alignments option cannot be selected with Alignment-free distance calculation",
            ),
        ]
        for condition, msg in errors:
            if condition:
                logging.log(logging.ERROR, msg)
                return True
        else:
            return False

    def show_incompatible_options_warnings(
        self,
    ) -> None:
        warnings = [
            (
                self.programstate.mode.get() == ProgramState.COMPARE_REFERENCE
                and self.programstate.perform_clustering.get(),
                'Clustering is not performed in the "Compare against reference" mode',
            ),
            (
                self.programstate.mode.get() == ProgramState.COMPARE_REFERENCE
                and self.programstate.print_alignments.get(),
                'Printing alignments is not implemented for "Compare against reference" mode',
            ),
            (
                self.programstate.mode.get() != ProgramState.COMPARE_REFERENCE
                and self.programstate.mode.get() != ProgramState.DECONTAMINATE
                and self.reference_file.get(),
                'A reference database is not needed in the selected mode and the selected reference database file will be ignored.',
            ),
        ]
        for condition, msg in warnings:
            if condition:
                logging.log(logging.WARNING, msg)

    def run_command(self) -> None:
        self.clear_command()
        self.update()
        with display_errors_and_warnings(debug=True):
            input_file = self.input_file.get()
            if self.show_incompatible_options_error():
                return
            self.show_incompatible_options_warnings()
            if self.programstate.mode.get() == ProgramState.COMPARE_REFERENCE:
                self.programstate.reference_comparison_process(
                    input_file, self.reference_file.get()
                )
            elif self.programstate.mode.get() == ProgramState.COMPARE_ALL:
                output_dir = self.preview_dir
                self.programstate.process(input_file)
                plot_input = os.path.join(
                    self.preview_dir, ProgramState.SUMMARY_STATISTICS_NAME
                )
                distance_name = [
                    distance
                    for distance, is_chosen in zip(
                        distances_short_names, self.programstate.distance_options
                    )
                    if is_chosen.get()
                ]
                if self.programstate.species_analysis:
                    self.show_progress("Starting plotting\n")
                    Plot(plot_input, output_dir, distance_name)
                    self.show_progress("Plotting complete\n")
            elif self.programstate.mode.get() == ProgramState.DEREPLICATE:
                self.programstate.dereplicate(input_file)
            elif self.programstate.mode.get() == ProgramState.DECONTAMINATE:
                self.programstate.decontaminate(input_file, self.reference_file.get())
            self.fill_file_list()
            tkmessagebox.showinfo("Done", "Calculation complete.")

    def clear_command(self) -> None:
        for file in os.scandir(self.preview_dir):
            if file.is_file():
                os.unlink(file)
        for file in os.scandir(os.path.join(self.preview_dir, "graph_previews")):
            if file.is_file():
                os.unlink(file)
        self.filelist.delete(*self.filelist.get_children())
        self.preview.delete("1.0", "end")
        self.preview_frame.configure(text="Preview")

    def outfilenames(self, which: str) -> Iterator[str]:
        if which == "all":
            index_list = self.filelist.get_children()
        elif which == "selected":
            index_list = self.filelist.selection()
        else:
            raise ValueError(f"Don't know how to save {which}")
        for index in index_list:
            yield self.filelist.item(index, option="text")

    def create_parameters_frame(self) -> None:
        parameters_frame_scroll = ttk.LabelFrame(self, text="Parameters")
        parameters_frame_scroll.rowconfigure(0, weight=1)
        parameters_frame_scroll.columnconfigure(0, weight=1)
        parameters_container = tk.Canvas(parameters_frame_scroll)
        parameters_frame = ttk.Frame(parameters_container)
        self.panes.add(parameters_frame_scroll, weight=0)
        parameters_frame.rowconfigure(12, weight=1)
        parameters_frame.columnconfigure(0, weight=1)
        parameters_frame.grid(row=0, column=0, sticky="nsew")
        parameters_container.grid(row=0, column=0, sticky="nsew")
        parameters_scroll = ttk.Scrollbar(
            parameters_frame_scroll, orient='vertical', command=parameters_container.yview)
        parameters_container.configure(yscrollcommand=parameters_scroll.set)
        parameters_scroll.grid(row=0, column=1, sticky="nsew")

        ttk.Label(parameters_frame, text="Input file format").grid(
            row=0, column=0, sticky="w"
        )

        format_combobox = ttk.Combobox(
            parameters_frame,
            textvariable=self.programstate.input_format_name,
            state="readonly",
            values=list(ProgramState.formats.keys()),
        )
        format_combobox.current(0)
        format_combobox.grid(row=1, column=0, sticky="w")

        distances_frm = ttk.LabelFrame(
            parameters_frame, text="Distances to calculate", padding=5, relief="sunken"
        )
        distances_frm.grid(row=2, column=0, sticky="we")

        for kind in range(NDISTANCES):
            ttk.Checkbutton(
                distances_frm,
                variable=self.programstate.distance_options[kind],
                text=distances_names[kind],
            ).grid(row=kind, column=0, sticky="w")

        ttk.Checkbutton(
            parameters_frame,
            variable=self.programstate.already_aligned,
            text="Already aligned",
        ).grid(row=3, column=0, sticky="w")
        ttk.Checkbutton(
            parameters_frame,
            variable=self.programstate.alignment_free,
            text="Alignment-free distance calculation",
        ).grid(row=4, column=0, sticky="w")
        ttk.Checkbutton(
            parameters_frame,
            variable=self.programstate.print_alignments,
            text="Print alignments",
        ).grid(row=5, column=0, sticky="w")

        ttk.Checkbutton(
            parameters_frame,
            variable=self.programstate.intra_species_lineages,
            text="Include 'intraspecific lineage' category",
        ).grid(row=6, column=0, sticky="w")

        ttk.Checkbutton(
            parameters_frame,
            variable=self.programstate.perform_clustering,
            text="Perform clustering",
        ).grid(row=7, column=0, sticky="w")
        ttk.Label(parameters_frame, text="Clustering by:").grid(
            row=8, column=0, sticky="w"
        )

        ttk.Combobox(
            parameters_frame,
            textvariable=self.programstate.cluster_distance,
            state="readonly",
            values=list(distances_names),
            width=30,
        ).grid(row=9, column=0, sticky="w")

        cluster_size_frame = ttk.Frame(parameters_frame)
        cluster_size_frame.grid(row=10, column=0, sticky="w")

        ttk.Label(
            cluster_size_frame, text="with distance threshold \n(between 0 and 1)"
        ).grid(row=0, column=0, sticky="w")
        ttk.Entry(cluster_size_frame, textvariable=self.programstate.cluster_size).grid(
            row=0, column=1, sticky="w"
        )

        # dereplicate_frame = ttk.LabelFrame(
        #     parameters_frame, text="Dereplicate parameters")
        # dereplicate_frame.grid(row=11, column=0, sticky="w")

        # similarity_frame = ttk.Frame(dereplicate_frame)
        # similarity_frame.grid(row=0, column=0, sticky="w")
        # ttk.Label(similarity_frame, text="Distance similarity threshold").grid(
        #     row=0, column=0, sticky="w")
        # ttk.Combobox(similarity_frame, values=[
        #              "0.07", "0.10", "0.25", "0.31"],
        #              textvariable=self.programstate.dereplicate_settings.similarity
        #              ).grid(row=0, column=1, sticky="w")

        # length_threshold_frame = ttk.Frame(dereplicate_frame)
        # length_threshold_frame.grid(row=1, column=0, sticky="w")
        # ttk.Label(length_threshold_frame, text="Remove sequences shorter than: ").grid(
        #     row=0, column=0, sticky="w")
        # ttk.Entry(length_threshold_frame,
        #           textvariable=self.programstate.dereplicate_settings.length_threshold).grid(
        #     row=0, column=1, sticky="w")

        # ttk.Radiobutton(dereplicate_frame, variable=self.programstate.dereplicate_settings.keep_most_complete,
        #                 text="Keep sequences with smallest amount of missing data").grid(row=2, column=0, sticky="w")

        # ttk.Radiobutton(dereplicate_frame, variable=self.programstate.dereplicate_settings.save_excluded_replicates,
        #                 text="Save excluded replicates to separate file").grid(row=3, column=0, sticky="w")

    def create_filelist_frame(self) -> None:
        filelist_frame = ttk.Labelframe(self, text="Files")
        filelist_frame.rowconfigure(0, weight=1)
        filelist_frame.columnconfigure(0, weight=1)
        self.panes.add(filelist_frame, weight=0)

        self.filelist = ttk.Treeview(
            filelist_frame, height=15, selectmode="extended", show="tree"
        )
        self.filelist.grid(row=0, column=0, sticky="nsew")

        filelist_scroll = ttk.Scrollbar(
            filelist_frame, orient="vertical", command=self.filelist.yview
        )
        self.filelist.configure(yscrollcommand=filelist_scroll.set)
        filelist_scroll.grid(row=0, column=1, sticky="nsew")

        filelist_scroll_x = ttk.Scrollbar(
            filelist_frame, orient="horizontal", command=self.filelist.xview
        )
        self.filelist.configure(xscrollcommand=filelist_scroll_x.set)
        filelist_scroll_x.grid(row=1, column=0, sticky="nsew")

        self.filelist.bind("<<TreeviewSelect>>", self.preview_selected)

    def icon_for_file(self, filename) -> tk.PhotoImage:
        TXT_EXTS = {".txt", ".tab", ".tsv", ".csv", ".spart"}
        _, ext = os.path.splitext(filename)
        if ext in TXT_EXTS:
            return self.images["txt_icon"]
        elif ext == ".log":
            return self.images["log_icon"]
        else:
            return self.images["graph_icon"]

    def fill_file_list(self) -> None:
        def by_ext(name):
            name, ext = os.path.splitext(name)
            return (ext, name)

        files = sorted(
            (file for file in os.listdir(self.preview_dir) if file != "graph_previews"),
            key=by_ext,
        )

        for filename in files:
            name = os.path.basename(filename)
            img = self.icon_for_file(name)
            self.filelist.insert(parent="", index="end", text=name, image=img)

    def create_preview_frame(self) -> None:
        self.preview_frame = ttk.LabelFrame(self, text="Preview")
        self.preview_frame.rowconfigure(0, weight=1)
        self.preview_frame.columnconfigure(0, weight=1)
        self.panes.add(self.preview_frame, weight=1)

        self.preview = tk.Text(self.preview_frame, height=15, width=30, wrap="none")
        self.preview.grid(row=0, column=0, sticky="nsew")

        yscroll = ttk.Scrollbar(
            self.preview_frame, orient="vertical", command=self.preview.yview
        )
        self.preview.config(yscrollcommand=yscroll.set)
        yscroll.grid(row=0, column=1, sticky="nsew")

        xscroll = ttk.Scrollbar(
            self.preview_frame, orient="horizontal", command=self.preview.xview
        )
        self.preview.config(xscrollcommand=xscroll.set)
        xscroll.grid(row=1, column=0, sticky="nsew")

    def preview_selected(self, _) -> None:
        self.preview.delete("1.0", "end")
        if not self.filelist.selection():
            return
        selected_index = self.filelist.selection()[-1]
        self.preview_frame.configure(
            text=f'Preview - {self.filelist.item(selected_index, option="text")}'
        )
        file_to_preview = self.filelist.item(selected_index, option="text")
        full_file_to_preview = os.path.join(self.preview_dir, file_to_preview)
        TXT_EXTS = {".txt", ".tab", ".tsv", ".csv", ".log", ".spart"}
        IMG_EXTS = {".gif", ".png", ".pbm", ".pgm", ".ppm", ".pnm"}
        _, ext = os.path.splitext(file_to_preview)
        if ext in TXT_EXTS:
            self.preview_txt(full_file_to_preview)
        elif ext in IMG_EXTS:
            self.preview_img(full_file_to_preview)
        elif ext == ".pdf":
            self.preview_pdf(file_to_preview)
        else:
            self.no_preview(full_file_to_preview)

    def preview_txt(self, filename) -> None:
        with open(filename) as file:
            self.preview.insert("1.0", file.read())

    def preview_img(self, filename) -> None:
        self.images["current"] = tk.PhotoImage(file=filename)
        self.preview.image_create("1.0", image=self.images["current"])

    def preview_pdf(self, filename) -> None:
        name, _ = os.path.splitext(filename)
        image_path = os.path.join(self.preview_dir, "graph_previews", name + ".png")
        self.images["current"] = tk.PhotoImage(file=image_path)
        self.preview.insert("1.0", "Approximate preview.\n")
        self.preview.image_create("2.0", image=self.images["current"])

    def no_preview(self, _) -> None:
        self.preview.insert("1.0", "Preview is not possible")


def test_look() -> None:
    root = tk.Tk()
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)
    gui = TaxiGUI(root, preview_dir="/tmp/out_dir")
    gui.fill_file_list()
    root.mainloop()
