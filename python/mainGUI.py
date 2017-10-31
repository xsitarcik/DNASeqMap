import tkinter as tk
from subprocess import call


class MainMenu(tk.Frame):
    count = 0
    def __init__(self, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        b = tk.Button(self, text="Compression Settings", command=self.new_window)
        b.pack(side="top")
        b1 = tk.Button(self, text="Run program", command=self.run_program)
        b1.pack(side="top")

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
        call(["SequenceAligner/myapp", "args", "to", "spa"])

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Main Menu")
    root.geometry("250x250")
    view = MainMenu(root)
    view.pack(side="top", fill="both", expand=True)
    root.mainloop()