import tkinter as tk
from  tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter import messagebox
from subprocess import call
from functools import partial


class MainMenu(tk.Frame):
    count = 0
    genome_file = None

    def __init__(self, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        compression_settings_opener = tk.Button(self, text="Compression Settings", command=self.new_window)
        compression_settings_opener.place(x=200,y=200)

        var_genome = StringVar()
        label_genome = Label(self, textvariable=var_genome, relief=RAISED, padx = 4,pady=4)
        var_genome.set("!No genome loaded!")
        label_genome.place(x=120,y=20)

        label_genome_info = Label(self,text="Genome path", relief=RAISED,padx=4,pady=4)
        label_genome_info.place(x=20,y=20)

        genome_loader = tk.Button(self, text="Load genome", command=partial(self.load_genome, var_genome))
        genome_loader.place(x=200,y=300)
        program_opener = tk.Button(self, text="Run program", command=self.run_program)
        program_opener.place(x=200,y=400)

    def load_genome(self,var_genome):
        root.withdraw()
        self.genome_file = askopenfilename()
        var_genome.set(self.genome_file)
        root.deiconify()

    def new_window(self):
        self.count += 1
        id = "New window #%s" % self.count
        window = tk.Toplevel(self)
        label = tk.Label(window, text=id)
        label.pack(side="top", fill="both", padx=10, pady=10)
        b = tk.Button(window, text="Confirm", command=window.destroy)
        b.pack(side="top")

        window.lift()
        window.focus_force()
        window.grab_set()

    def run_program(self):
        if self.genome_file is None:
            messagebox.showerror("No file for genome set", "Please load genome")
        else:
            call(["../myapp", self.genome_file, "to", "spa"])

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Main Menu")
    root.geometry("700x500")
    view = MainMenu(root)
    view.pack(side="top", fill="both", expand=True)
    root.mainloop()