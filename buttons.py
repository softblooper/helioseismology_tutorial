import sys
import os
import tkinter
from tkinter import *

root = Tk() #Blank window
theLabel = Label(root, text="Select a Tutorial") # Creates labels
theLabel.pack()
topFrame = Frame(root)
topFrame.pack()
bottomFrame = Frame(root)
bottomFrame.pack(side=BOTTOM)

def CallBack1():
    os.system("01_read_data.py")
button1 = Button(topFrame, text="Tutorial 1", fg="black", command = CallBack1)
button1.pack(side=LEFT)

def CallBack2():
    os.system("compare_frames_c2.py")
button2 = Button(topFrame, text="Tutorial 2", fg="black", command = CallBack2)
button2.pack(side=LEFT)

def CallBack3():
    os.system("03_acoustic_signal_c2.py")
button3 = Button(topFrame, text="Tutorial 3", fg="black", command = CallBack3)
button3.pack(side=LEFT)

def CallBack4():
    os.system("04_time_slices.py")
button4 = Button(topFrame, text="Tutorial 4", fg="black", command = CallBack4)
button4.pack(side=LEFT)

def CallBack5():
    os.system("05_time_averages.py")
button5 = Button(topFrame, text="Tutorial 5", fg="black", command = CallBack5)
button5.pack(side=LEFT)

def CallBack6():
    os.system("06_residuals.py")
button6 = Button(topFrame, text="Tutorial 6", fg="black",command = CallBack6)
button6.pack(side=LEFT)

def CallBack7():
    os.system("07_power_spectra.py")
button7 = Button(topFrame, text="Tutorial 7", fg="black", command = CallBack7)
button7.pack(side=LEFT)


root.mainloop() #Keeps the window running
