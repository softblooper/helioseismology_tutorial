# using Tkinter's Optionmenu() as a combobox
    # Python3
import tkinter as tk
def select():
    sf = "value is %s" % var.get()
    root.title(sf)
root = tk.Tk()
# use width x height + x_offset + y_offset (no spaces!)
root.geometry("%dx%d+%d+%d" % (330, 80, 200, 150))
root.title("tk.Optionmenu as combobox")
var = tk.StringVar(root)
# initial value
var.set('Tutorial 1')
choices = ['Tutorial 1', 'Tutorial 2', 'Tutorial 3', 'Tutorial 4','Tutorial 5', 'Tutorial 6', 'Tutorial 7']
option = tk.OptionMenu(root, var, *choices)
option.pack(side='left', padx=10, pady=10)
button = tk.Button(root, text="Select Tutorial", command=select)
button.pack(side='left', padx=20, pady=10)
root.mainloop()
