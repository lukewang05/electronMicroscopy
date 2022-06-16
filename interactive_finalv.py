# Generates strain maps in the material based on the diffraction data acquired using
# Precession Electron Diffraction (PED)
# May also be used for 4D-STEM data
# The first part of the algorithm filters the PED data
# The second part calculates the distance in diffraction patterns
# The third part generates strain maps


from os import remove, path
import plotly.express as px
from concurrent.futures import ProcessPoolExecutor
from csv import writer
from ctypes import c_double
from hyperspy.api import load
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import cm
from matplotlib import colors
import matplotlib.pyplot as plt
from multiprocessing import Array
from numpy import sqrt, array, ndenumerate, arange, percentile, linspace, zeros
from numpy.ctypeslib import as_array
from pandas import DataFrame
from PIL import Image, ImageTk
from pixstem.api import PixelatedSTEM
from seaborn import heatmap
import tkinter as tk
from tkinter import filedialog
from skimage.feature import match_template

file = None
distances = None
single_values = None
curr_func = None


def set_curr_func(func_name, current_file, s_values):
    global curr_func
    curr_func = str(func_name)
    entry.delete(0, tk.END)
    if curr_func == "load_file":
        entry.bind("<Return>", get_entry(curr_func))
        label1['text'] = label1['text'] \
            + "Please enter the path of the input file in the text box provided then press Enter.\n"
    elif curr_func == "to_csv":
        if current_file is None:
            label1['text'] = label1['text'] + "Please load a file before saving data.\n"
        elif s_values is None:
            label1['text'] = label1['text'] + "Please analyze the file before saving data.\n"
        else:
            entry.bind("<Return>", get_entry(curr_func))
            label1['text'] = label1['text'] \
                + "Please enter the path of the file you want to save to in the text box provided then press Enter.\n"
    elif curr_func == "analysis":
        if current_file is None:
            label1['text'] = label1['text'] + "Please load a file before starting analysis.\n"
        else:
            entry.bind("<Return>", get_entry(curr_func))
            label1['text'] = label1['text'] \
                + "Please enter the number of rows and columns you would like to analyze, as integers, separated by" \
                  " spaces. Press Enter when ready.\n"


def get_entry(current_function):
    # global curr_func  # replaced global variable with passed through variable
    if current_function == "load_file":
        entry.unbind("<Return>")
        load_file(entry.get())
    elif current_function == "analysis":
        entry.unbind("<Return>")
        start_analysis(entry.get())
    elif current_function == "to_csv":
        entry.unbind("<Return>")
        to_csv(entry.get())


def load_file():
    global file
    input_file = filedialog.askopenfilename()
    label1['text'] = label1['text'] + "Loading file...\n"
    root.update()
    try:
        file = load(input_file)
        label1['text'] = label1['text'] + "File loaded.\n"
    except:
        label1['text'] = label1['text'] + "Error loading. Please check path and try again.\n"
    entry.delete(0, tk.END)
    # entry.unbind("<Return>")


def distance(x1, y1, x2, y2):
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2))


def intensity(values):
    s = PixelatedSTEM(file.inav[values[0], values[1]])
    im_array = array(s)
    single_values[values[1]][values[0]] = im_array[values[3]][values[2]]


def find_center(im, peak):
    center = (352, 382)
    minimum = 144
    for (x, y) in ndenumerate(peak):
        for (a, b) in y:
            length = len(im)
            d = distance(350, 380, b, a)
            # d = distance(length/2, length/2, b, a)
            if (int(a) & int(b)) > 340 and d < minimum and int(a) < 390 and int(b) < 370:
                minimum = d
                center = (b, a)
    return center


def multiprocessing_func(values):
    global single_values, distances
    s = PixelatedSTEM(file.inav[values[0], values[1]])

    original = array(s)
    ##################################################################################################################
    # # FILTERS

    # sigma_est = mean(estimate_sigma(original, ))
    # # patch_size for 580 - 1, patch_size for 144 = 3
    # nlm = denoise_nl_means(original, h=1.15*sigma_est, fast_mode=True, patch_size=1, patch_distance=6, )
    # gaussian = gaussian_filter(original, 1.15*sigma_est)
    # wien = wiener(original, 5, 3)

    # original = array(nlm)
    #################################################################################################################
    # # PIXSTEM

    # s = s.rotate_diffraction(0,show_progressbar=False)
    # st = s.template_match_disk(disk_r=5, lazy_result=False, show_progressbar=False)
    # peak_array = st.find_peaks(lazy_result=False, show_progressbar=False)
    # peak_array_com = s.peak_position_refinement_com(peak_array, lazy_result=False, show_progressbar=False)
    # s_rem = s.subtract_diffraction_background(lazy_result=False, show_progressbar=False)
    # peak_array_rem_com = s_rem.peak_position_refinement_com(peak_array_com, lazy_result=False, show_progressbar=False)
    ##################################################################################################################
    # MY METHOD

    # defines template and templates matches
    # spot for 580 - [265:320, 265:320]
    # spot for 144 - [65:80, 65:80]
    template = original[265:320, 265:320]
    result = match_template(original, template, pad_input=True)
    # only takes points greater than the threshold r-value
    temp_list = []
    for i in range(len(result)):
        for j in range(len(result[i])):
            if result[i][j] > 0.87:  # change correlation value
                temp_list.append((i, j))
    # removes duplicate spots that are too close to each other
    i = 0
    list1 = []  # need better name original name was "l"
    while i < len(temp_list):
        j = 0
        temp = []
        point = temp_list[i]
        while j < len(temp_list):
            if distance(point[0], point[1], temp_list[j][0], temp_list[j][1]) < 15:  # change minimum center distance
                temp.append(temp_list[j])
                temp_list.pop(j)
            else:
                j = j + 1
        # max = 0  # variable not used
        pnt = temp[0]
        for j in range(len(temp)):
            if result[pnt[0]][pnt[1]] < result[temp[j][0]][temp[j][1]]:
                # max = result[temp[j][0]][temp[j][1]]  # function not used
                pnt = temp[j]
        list1.append(pnt)
    peak_array_rem_com = [[], list1]
    ##################################################################################################################
    center = find_center(original, peak_array_rem_com)
    # finds the specific spot and adding that distance to the array
    # pos_distance = 0  # variable not used
    closest_point = center
    idx = 0
    length = len(original)
    for (x, y) in ndenumerate(peak_array_rem_com):
        minimum = 999999
        for (a, b) in y:
            if 2 < b < length - 2 and 2 < a < length - 2:
                di = distance(center[0], center[1], b, a)
                distances[values[1]][values[0]][idx] = round(di, 3)
                idx += 1
            dis = distance(values[2], values[3], b, a)
            if dis < minimum and dis < length / 10:
                minimum = dis
                closest_point = (b, a)
    pos_distance = distance(closest_point[0], closest_point[1], center[0], center[1])
    single_values[values[1]][values[0]] = round(pos_distance, 4)
    print(values[0], values[1], closest_point, pos_distance, center)


def start_analysis(values=None):
    global file, curr_func
    point_xy = None
    analysis_method = ""
    if point_xy is None:
        def assign_method(method):
            nonlocal analysis_method
            analysis_method = method

        def mouse_coords(event):
            nonlocal analysis_method
            length = len(array(PixelatedSTEM(file.inav[25, 25])))
            mouse_xy = (int(event.x * length / 400), int(event.y * length / 400))  # get the mouse position from event
            label1['text'] = label1['text'] + str(mouse_xy[0]) + " " + str(mouse_xy[1]) + "\n"
            label1['text'] = label1['text'] + "Starting analysis...\n"
            r.update()
            analysis(mouse_xy, values, analysis_method)
            remove("temp.png")
            c2.unbind('<Button-1>')
            r.destroy()
            label1['text'] = label1['text'] + "Analysis complete.\n"

        s = PixelatedSTEM(file.inav[25, 25])
        s.save("temp.png")
        img = Image.open("temp.png")
        img = img.resize((400, 400), Image.ANTIALIAS)
        img.save('temp.png')

        r = tk.Toplevel(root)

        c = tk.Canvas(r, height=720, width=1080)
        c.pack()
        f = tk.Frame(r, bg='#333333')
        f.place(relwidth=1, relheight=1)
        log = tk.Message(f, bg='#999999', font=('Calibri', 15), anchor='nw', justify='left', highlightthickness=0, bd=0,
                         width=1000)
        log.place(relx=0.05, rely=0.7, relwidth=0.9, relheight=0.2)

        b1 = tk.Button(f, text='Intensity Mapping', bg='#620000', font=('Calibri', 15), highlightthickness=0, bd=0,
                       activebackground='#800000', activeforeground='#ffffff',
                       command=lambda: assign_method("intensity"), pady=0.02, fg='#ffffff')

        b1.place(relx=0.2, rely=0.6, relwidth=0.2, relheight=0.05)

        b2 = tk.Button(f, text='Strain Mapping', bg='#620000', font=('Calibri', 15), highlightthickness=0, bd=0,
                       activebackground='#800000', activeforeground='#ffffff', command=lambda: assign_method("strain"),
                       pady=0.02, fg='#ffffff')
        b2.place(relx=0.6, rely=0.6, relwidth=0.2, relheight=0.05)
        c2 = tk.Canvas(r, width=400, height=400)
        c2.place(relx=0.3)
        img = ImageTk.PhotoImage(Image.open("temp.png"))
        c2.create_image(0, 0, anchor='nw', image=img)
        c2.bind('<Button-1>', mouse_coords)
        log['text'] = log['text'] + "Please click on the method of analysis and then the point you would like to " \
            "analyze from the diffraction pattern above.\n"
        r.mainloop()
        if path.exists("temp.png"):
            remove("temp.png")


def analysis(point_xy, values, analysis_method=""):
    global file, single_values, distances
    t = values.split(" ")
    col = int(t[1])
    row = int(t[0])

    list1 = []  # need better variable name
    for r in range(row):
        for c in range(col):
            list1.append((r, c, point_xy[0], point_xy[1]))

    shared_array_base = Array(c_double, row * col)
    single_values = as_array(shared_array_base.get_obj())
    single_values = single_values.reshape(col, row)

    shared_array = Array(c_double, row * col * 50)
    distances = as_array(shared_array.get_obj())
    distances = distances.reshape(col, row, 50)

    with ProcessPoolExecutor() as executor:
        if analysis_method == "strain":
            executor.map(multiprocessing_func, list1)
        else:
            executor.map(intensity, list1)
    entry.delete(0, tk.END)
    f = open("Distances", "w")
    w = writer(f)
    for i in distances:
        w.writerow(i)
    f.close()
    label1['text'] = label1['text'] + "File saved.\n"
    entry.delete(0, tk.END)
    # entry.unbind("<Return>")


def to_csv(filename=None):
    global single_values
    f = open(filename, "w")
    w = writer(f)
    for i in single_values:
        w.writerow(i)
    f.close()
    label1['text'] = label1['text'] + "File saved.\n"
    entry.delete(0, tk.END)
    # entry.unbind("<Return>")


def heat_map_maker(minimum, maximum, parity=0):
    global single_values, distances

    if parity == 0:
        data = single_values.copy()
        df = DataFrame(data, columns=arange(len(data[0])), index=arange(len(data)))
        _, a = plt.subplots(figsize=(6, 5.5))
        chart1 = heatmap(df, cmap=cm.get_cmap("rainbow"), ax=a, vmin=minimum, vmax=maximum, square=True)
        return chart1.get_figure()
    else:
        data = zeros((len(single_values), len(single_values[0])), dtype=float)
        for i in range(len(distances)):
            for j in range(len(distances[i])):
                dist_sum = 0
                num = 0
                for k in distances[i][j]:
                    if minimum < k < maximum:
                        dist_sum += k
                        num += 1
                if num > 0:
                    data[i][j] = round(dist_sum / num, 1)

        df = DataFrame(data, columns=arange(len(data[0])), index=arange(len(data)))
        _, a = plt.subplots(figsize=(6, 5.5))
        gray = cm.get_cmap('gray', 512)
        new_colors = gray(linspace(0.15, 0.85, 2048))
        white = array([255 / 256, 255 / 256, 255 / 256, 1])
        new_colors[:1, :] = white
        new_colors[2047:, :] = white
        new_colormap = colors.ListedColormap(new_colors)
        chart = heatmap(df, cmap=new_colormap, vmin=minimum, vmax=maximum, square=True)
        return chart.get_figure()


def bar_chart():
    global distances
    if file is None:
        label1['text'] = label1['text'] + "Please load a file before creating a bar chart.\n"
    elif distances is None:
        label1['text'] = label1['text'] + "Please analyze the file before creating a bar chart.\n"
    else:
        label1['text'] = label1['text'] + "Creating bar chart. This might take several minutes depending " \
                                          "on the size of data.\n"
        root.update()
        dist = single_values.flatten()

        # fig, a = plt.subplots(figsize=(6,5.5))
        plt.xlabel('Distance from center peek', fontsize=10)
        plt.ylabel('Counts', fontsize=10)
        plt.title('Distance Counts', fontsize=10)
        # plt.bar(y_pos, counts, align='center', alpha=0.95) # creates the bar plot
        plt.hist(dist, bins=500)

        def scope_heat_map(chart):
            values = e.get().split(" ")
            minimum = float(values[0])
            maximum = float(values[1])
            f = heat_map_maker(minimum, maximum, 1)
            chart = FigureCanvasTkAgg(f, bar_chart_window)
            chart.draw()
            chart.get_tk_widget().place(relx=0.51, rely=0.2)

        bar_chart_window = tk.Toplevel(root)
        bar_chart_window.geometry('1920x1080')
        chart_type = FigureCanvasTkAgg(plt.gcf(), bar_chart_window)
        chart_type.draw()
        chart_type.get_tk_widget().place(relx=0.0, rely=0.2, relwidth=0.5)
        m = tk.Message(bar_chart_window, font=('Calibri', 15), highlightthickness=0, bd=0, width=1000, justify='center')
        m['text'] = "Enter the minimum value and the maximum value (exclusive) separated by a space. " \
                    "Press Enter to create the heatmap with these specifications"
        m.place(relx=0.25, rely=0.05)
        e = tk.Entry(bar_chart_window, font=('Calibri', 15))
        e.place(relx=0.44, rely=0.1)
        e.bind("<Return>", scope_heat_map(chart_type))


def outlier(data):
    data = data.flatten()
    q1 = percentile(data, 25)
    q3 = percentile(data, 75)
    iqr = q3 - q1
    minimum = q1 - (1.5 * iqr)
    maximum = q3 + (1.5 * iqr)
    return minimum, maximum


def heat_map():
    import hyperspy.api as hs
    global single_values
    if file is None:
        label1['text'] = label1['text'] + "Please load a file before creating a heat map.\n"
    elif distances is None:
        label1['text'] = label1['text'] + "Please analyze the file before creating a heat map.\n"
    else:
        data = single_values.copy()
        df = DataFrame(data, columns=arange(len(data[0])), index=arange(len(data)))
        print(df)
        print(single_values)
        fig = px.imshow(df, color_continuous_scale=["blue", "green", "red"])
        fig.show()

        def image_gallery():
            global file
            values = e.get().split(" ")
            x0 = int(values[0])
            y0 = int(values[1])
            x1 = int(values[2])
            y1 = int(values[3])
            # x = 0  # variable not used
            # y = 0
            index_x = x1

            for x in range(x1 - x0 + 1):
                index_y = y1
                for y in range(y1 - y0 + 1):
                    s = PixelatedSTEM(hs.signals.Signal2D(file.inav[index_x, index_y]))
                    st = s.template_match_ring(r_inner=1, r_outer=6, lazy_result=True, show_progressbar=False)
                    peak_array = st.find_peaks(method='dog', min_sigma=0.8, max_sigma=15, sigma_ratio=1.9,
                                               threshold=0.42, overlap=0.5, lazy_result=False, show_progressbar=True)
                    s.add_peak_array_as_markers(peak_array)
                    # plt.plot(s)
                    s.plot()
                    ax = s.plot.signal_plot.ax  # was originally "ax = s._plot.signal_plot.ax"
                    ax.set_xlabel("pixel(" + str(index_x) + "_" + str(index_y) + ")")
                    # plt.title("pixel(" + str(index_x) + "_" + str(index_y) + ")")
                    # plt.show()
                    index_y -= 1
                    y += 1
                index_x -= 1
                x += 1
            plt.show()

        bar_chart_window = tk.Toplevel(root)
        bar_chart_window.geometry('1280x720')
        m = tk.Message(bar_chart_window, font=('Calibri', 15), highlightthickness=0, bd=0, width=600, justify='left')
        m['text'] = "A new window should open displaying the heatmap created, if you would like to view specific " \
            "diffraction patterns, " "enter the starting x and the y value and the ending x and y value " \
            "separated by a space. " "Press Enter to display these diffraction patterns"
        m.place(relx=0.25, rely=0.05)
        e = tk.Entry(bar_chart_window, font=('Calibri', 15))
        e.place(relx=0.44, rely=0.3)
        e.bind("<Return>", image_gallery)

        bar_chart_window.mainloop()


if __name__ == "__main__":

    HEIGHT = 1080
    WIDTH = 1920

    root = tk.Tk()

    canvas = tk.Canvas(root, height=HEIGHT, width=WIDTH)
    canvas.pack()
    frame = tk.Frame(root, bg='#450000')
    frame.place(relwidth=1, relheight=1)

    # Menu Label
    label = tk.Label(frame, text='Menu', bg='#450000', font=('Times New Roman', 50), fg='#ffffff')
    label.place(relx=0.40, rely=0.05, relwidth=0.2, relheight=0.05)

    # Text Output box
    label1 = tk.Message(frame, bg='#ffffff', font=('Calibri', 15), anchor='nw', justify='left', highlightthickness=0,
                        bd=0, width=1500)
    label1.place(relx=0.1, rely=0.5, relwidth=0.8, relheight=0.35)

    # Entry box
    entry = tk.Entry(frame, font=('Calibri', 15))
    entry.place(relx=0.1, rely=0.9, relwidth=0.8, relheight=0.05)



    # Creating buttons for UI. The buttons change the current function to the button name.
    button = tk.Button(frame, text='Load File', bg='#620000', font=('Calibri', 30), highlightthickness=0, bd=0,
                       activebackground='#800000', activeforeground='#ffffff', 
                       command=lambda: load_file(), pady=0.02, fg='#ffffff')

    button.place(relx=0.42, rely=0.15, relwidth=0.16, relheight=0.05)

    button1 = tk.Button(frame, text='Start Analysis', bg='#620000', font=('Calibri', 30), highlightthickness=0, bd=0,
                        activebackground='#800000', activeforeground='#ffffff',
                        command=lambda: set_curr_func("analysis", file, single_values), pady=0.02, fg='#ffffff')
    button1.place(relx=0.39, rely=0.22, relwidth=0.22, relheight=0.05)

    button2 = tk.Button(frame, text='Create Bar Chart', bg='#620000', font=('Calibri', 30), highlightthickness=0, bd=0,
                        activebackground='#800000', activeforeground='#ffffff', command=lambda: bar_chart(), pady=0.02,
                        fg='#ffffff')
    button2.place(relx=0.375, rely=0.29, relwidth=0.25, relheight=0.05)

    button3 = tk.Button(frame, text='Create Heat Map', bg='#620000', font=('Calibri', 30), highlightthickness=0, bd=0,
                        activebackground='#800000', activeforeground='#ffffff', command=lambda: heat_map(), pady=0.02,
                        fg='#ffffff')
    button3.place(relx=0.38, rely=0.36, relwidth=0.24, relheight=0.05)

    button4 = tk.Button(frame, text='Transfer Data to .csv', bg='#620000', font=('Calibri', 30), highlightthickness=0,
                        bd=0, activebackground='#800000', activeforeground='#ffffff',
                        command=lambda: set_curr_func("to_csv", file, single_values), pady=0.02, fg='#ffffff')
    button4.place(relx=0.34, rely=0.43, relwidth=0.32, relheight=0.05)

    root.mainloop()
    if path.exists("temp.png"):
        remove("temp.png")
