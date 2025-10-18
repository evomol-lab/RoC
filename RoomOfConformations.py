import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext
from PIL import Image, ImageTk
import threading
import queue
import traceback
import sys
import os

# --- PyInstaller Hidden Import Fix ---
# This block forces PyInstaller to find hidden compiled modules.
try:
    import MDAnalysis.lib.formats.cython_util
    import MDAnalysis.lib._transformations
    import PIL._tkinter_finder
except ImportError:
    pass
# --- End of Fix ---


# --- Backend Scientific Libraries ---
import numpy as np
import prody as pr
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import plotly.express as px
import pandas as pd # THE FIX IS HERE: aliasing pandas as pd
import matplotlib
matplotlib.use('Agg') # Use a non-interactive backend for Matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# This helper function is CRITICAL for PyInstaller.
def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)


class ResultsWindow(tk.Toplevel):
    """A new window to display the Matplotlib plot images."""
    def __init__(self, master):
        super().__init__(master)
        self.title("Analysis Results")
        self.geometry("1000x500")

        # Create a main frame
        main_frame = ttk.Frame(self, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)
        main_frame.columnconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(0, weight=1)
        
        # --- RMSD Plot ---
        self.rmsd_fig = Figure(figsize=(5, 4), dpi=100)
        self.rmsd_ax = self.rmsd_fig.add_subplot(111)
        self.rmsd_canvas = FigureCanvasTkAgg(self.rmsd_fig, master=main_frame)
        self.rmsd_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew", padx=5)

        # --- RMSF Plot ---
        self.rmsf_fig = Figure(figsize=(5, 4), dpi=100)
        self.rmsf_ax = self.rmsf_fig.add_subplot(111)
        self.rmsf_canvas = FigureCanvasTkAgg(self.rmsf_fig, master=main_frame)
        self.rmsf_canvas.get_tk_widget().grid(row=0, column=1, sticky="nsew", padx=5)
        
        self.withdraw()
        self.protocol("WM_DELETE_WINDOW", self.withdraw)

    def draw_plots(self, rmsd_data, rmsf_data, rmsf_resids):
        """Draws the RMSD and RMSF data onto the Matplotlib canvases."""
        # --- Draw RMSD ---
        self.rmsd_ax.clear()
        self.rmsd_ax.plot(rmsd_data[:, 0], rmsd_data[:, 2])
        self.rmsd_ax.set_title("RMSD vs. Original Structure")
        self.rmsd_ax.set_xlabel("Model / Frame")
        self.rmsd_ax.set_ylabel("RMSD (Å)")
        self.rmsd_fig.tight_layout()
        self.rmsd_canvas.draw()

        # --- Draw RMSF ---
        self.rmsf_ax.clear()
        self.rmsf_ax.plot(rmsf_resids, rmsf_data)
        self.rmsf_ax.set_title("RMSF per Residue")
        self.rmsf_ax.set_xlabel("Residue ID")
        self.rmsf_ax.set_ylabel("RMSF (Å)")
        self.rmsf_fig.tight_layout()
        self.rmsf_canvas.draw()


class RoomOfConformationsApp:
    def __init__(self, master):
        self.master = master
        master.title("The Room of Conformations")
        master.geometry("700x650")

        self.input_file_path = tk.StringVar()
        self.output_file_path = tk.StringVar()
        self.num_confs_var = tk.StringVar(value="50")
        self.num_modes_var = tk.StringVar(value="15")
        self.log_queue = queue.Queue()
        
        self.save_html_var = tk.BooleanVar(value=True)
        self.save_csv_var = tk.BooleanVar(value=True)
        
        self.results_window = ResultsWindow(self.master)
        self.create_widgets()
        self.process_log_queue()

    def create_widgets(self):
        main_frame = ttk.Frame(self.master, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)
        try:
            logo_path = resource_path("RoC-Logo.png")
            original_image = Image.open(logo_path)
            resized_image = original_image.resize((100, 100), Image.Resampling.LANCZOS)
            self.logo_image = ImageTk.PhotoImage(resized_image)
            logo_label = ttk.Label(main_frame, image=self.logo_image)
            logo_label.pack(pady=10)
        except Exception as e:
            self.log(f"Warning: Could not load logo. Error: {e}")
        io_frame = ttk.LabelFrame(main_frame, text="1. Input & Output Files", padding=10)
        io_frame.pack(fill=tk.X, expand=True, pady=5)
        io_frame.columnconfigure(1, weight=1)
        ttk.Label(io_frame, text="Input File:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        ttk.Entry(io_frame, textvariable=self.input_file_path, state="readonly").grid(row=0, column=1, sticky="ew", padx=5, pady=5)
        ttk.Button(io_frame, text="Browse...", command=self.browse_input).grid(row=0, column=2, sticky="e", padx=5, pady=5)
        ttk.Label(io_frame, text="Save File As:").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        ttk.Entry(io_frame, textvariable=self.output_file_path, state="readonly").grid(row=1, column=1, sticky="ew", padx=5, pady=5)
        ttk.Button(io_frame, text="Browse...", command=self.browse_output).grid(row=1, column=2, sticky="e", padx=5, pady=5)
        params_frame = ttk.LabelFrame(main_frame, text="2. Parameters", padding=10)
        params_frame.pack(fill=tk.X, expand=True, pady=5)
        ttk.Label(params_frame, text="Number of Conformations:").pack(side=tk.LEFT, padx=5)
        ttk.Entry(params_frame, textvariable=self.num_confs_var, width=10).pack(side=tk.LEFT, padx=5)
        ttk.Label(params_frame, text="Number of Modes:").pack(side=tk.LEFT, padx=5)
        ttk.Entry(params_frame, textvariable=self.num_modes_var, width=10).pack(side=tk.LEFT, padx=5)
        options_frame = ttk.LabelFrame(main_frame, text="3. Output Options", padding=10)
        options_frame.pack(fill=tk.X, expand=True, pady=5)
        ttk.Checkbutton(options_frame, text="Save interactive plots (.html)", variable=self.save_html_var).pack(side=tk.LEFT, padx=10)
        ttk.Checkbutton(options_frame, text="Save analysis data (.csv)", variable=self.save_csv_var).pack(side=tk.LEFT, padx=10)
        run_frame = ttk.LabelFrame(main_frame, text="4. Run", padding=10)
        run_frame.pack(fill=tk.X, expand=True, pady=5)
        self.run_button = ttk.Button(run_frame, text="Run Ensemble Generation", command=self.start_run_thread, state="disabled")
        self.run_button.pack(pady=5, fill=tk.X, expand=True)
        log_frame = ttk.LabelFrame(main_frame, text="Log", padding=10)
        log_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        self.log_text = scrolledtext.ScrolledText(log_frame, wrap=tk.WORD, height=10)
        self.log_text.pack(fill=tk.BOTH, expand=True)
        self.log_text.configure(state='disabled')

    def browse_input(self):
        file_path = filedialog.askopenfilename(title="Select a PDB or mmCIF file", filetypes=(("PDB files", "*.pdb"), ("mmCIF files", "*.cif"), ("All files", "*.*")))
        if file_path: self.input_file_path.set(file_path); self.check_ready_to_run()

    def browse_output(self):
        file_path = filedialog.asksaveasfilename(title="Save Ensemble As", filetypes=(("PDB files", "*.pdb"),), defaultextension=".pdb")
        if file_path: self.output_file_path.set(file_path); self.check_ready_to_run()
            
    def check_ready_to_run(self):
        if self.input_file_path.get() and self.output_file_path.get(): self.run_button.config(state="normal")
        else: self.run_button.config(state="disabled")

    def log(self, message):
        self.log_queue.put(message)

    def process_log_queue(self):
        try:
            while True:
                msg = self.log_queue.get_nowait()
                if isinstance(msg, tuple) and msg[0] == "PLOT_DATA":
                    self.results_window.draw_plots(msg[1], msg[2], msg[3])
                elif msg == "JOB_COMPLETE":
                    self.run_button.config(state="normal")
                    self.log("✅ Analysis complete. Showing results window.")
                    self.results_window.deiconify()
                elif msg == "JOB_FAILED":
                    self.run_button.config(state="normal")
                    self.log("🚨 ERROR: Analysis failed. Check log for details.")
                else:
                    self.log_text.configure(state='normal')
                    self.log_text.insert(tk.END, str(msg) + '\n')
                    self.log_text.configure(state='disabled')
                    self.log_text.see(tk.END)
        except queue.Empty:
            pass
        finally:
            self.master.after(100, self.process_log_queue)

    def start_run_thread(self):
        self.run_button.config(state="disabled")
        self.log_text.configure(state='normal'); self.log_text.delete(1.0, tk.END); self.log_text.configure(state='disabled')
        self.log("Starting job...")
        self.results_window.withdraw()
        try:
            in_file, out_file, num_confs, num_mds, save_html, save_csv = (
                self.input_file_path.get(), self.output_file_path.get(),
                int(self.num_confs_var.get()), int(self.num_modes_var.get()),
                self.save_html_var.get(), self.save_csv_var.get()
            )
            self.analysis_thread = threading.Thread(
                target=self.run_backend,
                args=(in_file, out_file, num_confs, num_mds, save_html, save_csv),
                daemon=True
            )
            self.analysis_thread.start()
        except ValueError:
            self.log("🚨 ERROR: Parameters must be integers."); self.run_button.config(state="normal")
        except Exception as e:
            self.log(f"🚨 ERROR: {e}"); self.run_button.config(state="normal")

    def run_backend(self, in_file, out_file, num_confs, num_mds, save_html, save_csv):
        try:
            self.log("─" * 80); self.log("GENERATING THE ALL-ATOM ENSEMBLE ⚙️")
            self.log("--> STEP 1/4: Loading structure and building ANM...")
            protein_prody = pr.parsePDB(in_file)
            calphas_prody = protein_prody.select('calpha')
            if calphas_prody is None: raise ValueError("ProDy could not select C-alpha atoms.")
            base_name = os.path.splitext(os.path.basename(in_file))[0]
            anm = pr.ANM(f'{base_name} ANM'); anm.buildHessian(calphas_prody); anm.calcModes(n_modes=num_mds)
            self.log("--> STEP 2/4: Generating C-alpha ensemble...")
            calpha_ensemble = pr.sampleModes(anm, atoms=calphas_prody, n_confs=num_confs)
            self.log("--> STEP 3/4: Loading structure with MDAnalysis...")
            u = mda.Universe(in_file); all_atoms_mda = u.select_atoms('all'); ref_calphas_mda = u.select_atoms('name CA')
            self.log(f"--> STEP 4/4: Reconstructing and writing {num_confs} models...")
            with mda.Writer(out_file, multiframe=True) as pdb_writer:
                for i in range(calpha_ensemble.numConfs()):
                    if (i+1) % 50 == 0: self.log(f"    ... writing model {i+1} of {num_confs}")
                    target_calpha_coords = calpha_ensemble[i].getCoords()
                    ref_calphas_mda.positions = target_calpha_coords
                    mda.analysis.align.alignto(all_atoms_mda, ref_calphas_mda, select='name CA', weights='mass')
                    pdb_writer.write(all_atoms_mda)
            self.log(f"\n🎉 Ensemble generation complete! File saved as '{out_file}'."); self.log("─" * 80)

            self.log("\nANALYZING AND PLOTTING RESULTS 📈")
            ens_universe = mda.Universe(out_file); calphas = ens_universe.select_atoms("name CA")
            ref_structure = ens_universe; ref_structure.trajectory[0]

            self.log("Calculating RMSD...")
            R = rms.RMSD(ens_universe, ref_structure, select="name CA", ref_frame=0).run()
            rmsd_results = R.results.rmsd
            if save_html:
                px.line(x=rmsd_results[:, 0], y=rmsd_results[:, 2], labels={'x': 'Model / Frame', 'y': 'RMSD (Å)'}, title='RMSD vs. Original Structure').write_html(os.path.join(os.path.dirname(out_file), "rmsd_plot.html"))
                self.log("Interactive RMSD plot (.html) saved.")
            if save_csv:
                pd.DataFrame(rmsd_results, columns=["Frame", "Time (ps)", "RMSD (A)"]).to_csv(os.path.join(os.path.dirname(out_file), "rmsd_data.csv"), index=False)
                self.log("RMSD data (.csv) saved.")

            self.log("Calculating RMSF...")
            align.AlignTraj(ens_universe, ref_structure, select='name CA', in_memory=True, ref_frame=0).run()
            R_fluct = rms.RMSF(calphas).run()
            rmsf_results = R_fluct.results.rmsf
            if save_html:
                px.line(x=calphas.resids, y=rmsf_results, labels={'x': 'Residue ID', 'y': 'RMSF (Å)'}, title='RMSF per Residue').write_html(os.path.join(os.path.dirname(out_file), "rmsf_plot.html"))
                self.log("Interactive RMSF plot (.html) saved.")
            if save_csv:
                pd.DataFrame({"Residue_ID": calphas.resids, "RMSF (A)": rmsf_results}).to_csv(os.path.join(os.path.dirname(out_file), "rmsf_data.csv"), index=False)
                self.log("RMSF data (.csv) saved.")
            
            # Send the raw data to the GUI for Matplotlib plotting
            self.log_queue.put(("PLOT_DATA", rmsd_results, rmsf_results, calphas.resids))
            self.log_queue.put("JOB_COMPLETE")

        except Exception as e:
            self.log("\n" + "="*80); self.log("🚨 AN ERROR OCCURRED! The process stopped unexpectedly. 🚨")
            self.log(f"ERROR: {e}"); self.log(traceback.format_exc()); self.log("="*80)
            self.log_queue.put("JOB_FAILED")

if __name__ == "__main__":
    root = tk.Tk()
    app = RoomOfConformationsApp(root)
    root.mainloop()
