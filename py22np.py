import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j1
from scipy.optimize import curve_fit
import matplotlib

def file_extraction(filename):
    """
    Function that takes in input file and then does the following:
    - Converts the metadata into a dictionary for later use in the code
    :param filename:
    :return:
    """
    metadata = dict()

    """
    Code sequence below copied from Computing 2 resource videos on Minerva
    """

    try:
        with open(filename, "r") as intensity_data:
            for line in intensity_data:
                line = line.strip()
                if line == "&END":
                    break
                if "=" in line:  # o this will be some metadata
                    parts = line.split("=")
                    if len(parts) != 2:
                        raise RuntimeError(f"Unable to parse {line} as name:value")
                    metadata[parts[0]] = parts[1]

            intensity_data.close()
    except:
        raise TypeError("File is not readable/has no metadata")

    return metadata


def data_extraction(filename):
    """

    :param filename:
    :return:
    """
    global meta, wavelength, k, units, symbol
    meta = file_extraction(filename)
    units = str(meta["Horizontal, Vertical units"])
    symbol = {"m": 0, "cm": -2, "mm": -3}
    wavelength = (float(meta["Wavelength (nm)"])) * (1e-9)
    k = (2 * np.pi) / wavelength

    word = "&END"
    with open(filename, 'r') as fp:
        # read all lines in a list
        lines = fp.readlines()
        for line in lines:
            # check if string present on a current line
            if line.find(word) != -1:
                position = lines.index(line)

    data = np.genfromtxt(filename, delimiter="", skip_header=(position + 1))

    length = len(data)

    # Defining a global 2D array containing intensity values from data file
    Intensity = []

    # Defining a temporary set to store intensity values

    # For loop to iterate through k-th row of 'data'
    for k in range(1, length):
        # Defining a temporary list to store intensity values
        Temp_Intensity = []
        for l in range(1, length):
            Temp_Intensity.append(data[l][k])

        # Append Temp_Intensity values to the Intensity list to make 2D array
        Intensity.append(Temp_Intensity)

    # Checking if x/y values are given as sin(angles)
    if meta["As Projected"] == "False":

        # Making an empty list to store y-values
        y_data = []

        # Iterating through data set to add y-values to the y-value list
        for i in range(1, length):
            y_data.append(data[0][i])

        # Making an empty list to store x-values
        x_data = []

        # Iterating through data set to add y-values to the x-value list
        for j in range(1, length):
            x_data.append(data[j][0])

    # Checking if the x/y values are given as positions on the screen
    elif meta["As Projected"] == "True":

        L = float(meta["Distance to Screen (m)"])

        if L == 0:
            raise ZeroDivisionError("Length cannot be 0 as division by 0 is undefined")

        # Making an empty list to store y-values
        y_data = []

        # Iterating through data set to add y-values to the y-value list
        for i in range(1, length):
            y_data.append(data[0][i])

        # Making an empty list to store x-values
        x_data = []

        # Iterating through data set to add y-values to the x-value list
        for j in range(1, length):
            x_data.append(data[j][0])

        # Making temporary x-data list
        x_temp = []

        # Adding every element in x_data to x_temp divided by length L in correct units [m]
        for el in x_data:
            x_temp.append((el * 10 ** (symbol[units])) / L)

        # Clearing x_data list
        x_data.clear()

        # Redefining x_data list
        x_data = x_temp

        # Making temporary y-data list
        y_temp = []

        # Adding every element in y_data to y_temp divided by length L in correct units [m]
        for el in y_data:
            y_temp.append((el * 10 ** (symbol[units])) / L)

        # Clearing x_data list
        y_data.clear()

        # Redefining x_data list
        y_data = y_temp

    else:
        raise TypeError("Metadata contains no information on projection")

    # Defining an empty set for the line sin(phi) = 0
    phi = []

    # Defining an empty set for the line sin(theta) = 0
    theta = []

    # Check if 0 is in the x_data set
    if 0 in y_data:
        # Storing index of 0 in x_data to position variable
        positiony = y_data.index(0)
        for m in range(1, length):
            # Appending the values in the position-th column of the data array
            phi.append(data[m][positiony])
    else:
        # Check if the length of the x/y set is even
        if length % 2 == 0:
            positiony = int(length / 2)
            for m in range(1, length):
                phi.append(data[m][positiony])
        else:
            positiony = int((length + 1) / 2)
            for m in range(1, length):
                phi.append(data[m][positiony])

    # Check if 0 is in the x_data set
    if 0 in x_data:
        # Storing index of 0 in x_data to position variable
        positionx = x_data.index(0)
        for m in range(1, length):
            # Appending the values in the position-th column of the data array
            theta.append(data[positionx][m])
    else:
        # Check if the length of the x/y set is even
        if length % 2 == 0:
            positionx = int(length / 2)
            for m in range(1, length):
                theta.append(data[positionx][m])
        else:
            positionx = int((length + 1) / 2)
            for m in range(1, length):
                theta.append(data[positionx][m])

    return x_data, y_data, Intensity, phi, theta


def k_converter(x_val):
    """
    Function that takes in data-set and returns the same set multiplied by the wave number k
    :param x_val:
    :return:
    """

    # Calling file_extraction function to retrieve the value of wavelength from metadata
    meta = file_extraction(filename)

    wavelength = float(meta["Wavelength (nm)"]) * 1e-9
    k = 2 * np.pi / wavelength

    # Converting given x_val array into a numpy array
    x_val = np.asarray(x_val)

    return k * x_val


def square(x, i0, width):
    """
    Function that creates an array / list of the intensity values that one would expect were the aperture a square shape
    :param width: Float
    :param i0: Float
    :param x: List
    :return: y-values
    """

    # Calling k_converter function to multiply elements in x by the wave-number k
    kx = k_converter(x)

    return i0 * (np.sinc(kx * (width / 2)) ** 2)


def diamond(x, io, width):
    """
    Function that creates an array / list of the intensity values that one would expect were the aperture a diamond shape
    :param x: List
    :param io: Float
    :param width: Float
    :return: y-values
    """

    # Calling k_converter function to multiply elements in x by the wave-number k
    kx = k_converter(x)

    return io * (np.sinc(kx * width / (2 * np.sqrt(2)))) ** 4


def circle(x, io, D):
    """
    Function that creates an array / list of the intensity values that one would expect were the aperture a circle
    :param x: List
    :param io: Float
    :param D: Float
    :return: y-values
    """

    # Calling k_converter function to multiply elements in x by the wave-number k
    kx = k_converter(x)

    # Changing the kx array value to the intensity estimate to prevent division by 0
    kx[kx == 0] = io

    return io * ((2 * j1(np.pi * kx * D)) / (np.pi * kx * D)) ** 2


def chi_squared(fitted_data, dataSet):
    """
    Function that takes two lists, namely the intensity values from the data file slices, and the y-values from the
    fitting functions and returns a value of the chi-squared parameter.
    :param fitted_data: List
    :param dataSet: List
    :return: Total
    """

    # Defining variable total to which chi-squared values for each data point are added
    total = 0

    for i in range(len(dataSet)):
        total += ((dataSet[i] - fitted_data[i]) ** 2)

    total = total / (len(dataSet) - 2 - 1)

    return total


def plotting_contour(x, y, z):
    """
    This function receives the x, y and z data from the data file and accordingly plots a contour/heat map showing
    the intensity values as functions of x and y values
    :param y: List
    :param x: List
    :param z: array
    :return: Contour plot
    """

    variable = plt.pcolor(z, norm=matplotlib.colors.LogNorm())

    # Plotting the main contour graph for the intensity data
    plt.pcolormesh(x, y, z)
    if meta["As Projected"] == "False":
        plt.xlabel(r"$sin(\theta)/ rad$")
        plt.ylabel(r"$sin(\phi)/ rad$")
        plt.title("Diffraction pattern data for py22np")
        plt.xlim([min(x), max(x)])
        plt.ylim([min(y), max(y)])
        plt.colorbar(variable)
        plt.show()
    elif meta["As Projected"] == "True":
        plt.xlabel("Distance, x/ mm")
        plt.ylabel("Distance, y/ mm")
        plt.title("Diffraction pattern data for py22np")
        plt.xlim([min(x), max(x)])
        plt.ylim([min(y), max(y)])
        plt.colorbar(variable)
        plt.show()
    else:
        raise TypeError("Metadata 'As Projected' setting not as expected")


def parameters(x, y):
    """
    Function that takes x and y values, i.e., distance and intensity, and returns the fitting parameters for
    I0 and the width of the aperture
    :param x: List :param y:
    :return: Fitting parameters
    """

    # Curve fitting for parameters of the square function
    poptsq, pcovsq = curve_fit(square, x, y, (max(y), 25 * (10 ** -6)))
    perrsq = np.sqrt(np.diag(pcovsq))

    # Curve fitting for parameters of the diamond function
    poptd, pcovd = curve_fit(diamond, x, y, (max(y), 25 * (10 ** -6)))
    perrd = np.sqrt(np.diag(pcovd))

    # Curve fitting for parameters of the circle function
    poptc, pcovc = curve_fit(circle, x, y, (max(y), 25 * (10 ** -6)))
    perrc = np.sqrt(np.diag(pcovc))

    # Adding all elements of the fitting parameters and their uncertainties to a list that can be accessed late
    fit_values = [poptsq, perrsq, poptd, perrd, poptc, perrc]

    return fit_values


def intensity_plots(x, y, z, phi, theta):
    """
    Function that receives arguments from the function titled data_extraction and from that does the following:
    - creates a new list that contains the diagonal entries of the intensity slice
    - creates 3 dictionaries to store the chi-squared value for each fit on each axis
    - determines the smallest chi-squared and hence finds the shape
    - plots two intensity cross-section graphs each of which indicate the smallest chi-squared value and the set of
    raw data plus the fits of the three other fitting functions
    :param x:
    :param y:
    :param z:
    :param phi:
    :param theta:
    :return:
    """

    # Defining an empty list to hold diagonal intensity values
    pos_diag_values = []

    # Adding diagonal elements of the Intensity array to pos_diag_values
    for i in range(len(x)):
        pos_diag_values.append(z[len(x) - 1 - i][i])

    # Calling parameters function to return an array of fitted parameters for fitting
    diag1 = parameters(x, phi)
    vertical = parameters(y, theta)
    horizontal = parameters(x, phi)

    # Defining dictionary to compare chi-squared values for diagonal slice
    positive_diag_chi_squared = {"square": chi_squared(square(x, diag1[0][0], diag1[0][1]), pos_diag_values),
                                 "diamond": chi_squared(diamond(x, diag1[2][0], diag1[2][1]), pos_diag_values),
                                 "circle": chi_squared(circle(x, diag1[4][0], diag1[4][1]), pos_diag_values)}

    # Defining dictionary to compare chi-squared values for vertical slice
    vert = {"square": chi_squared(square(y, vertical[0][0], vertical[0][1]), theta),
            "diamond": chi_squared(diamond(y, vertical[2][0], vertical[2][1]), theta),
            "circle": chi_squared(circle(y, vertical[4][0], vertical[4][1]), theta)}

    # Defining dictionary to compare chi-squared values for horizontal slice
    hori = {"square": chi_squared(square(x, horizontal[0][0], horizontal[0][1]), phi),
            "diamond": chi_squared(diamond(x, horizontal[2][0], horizontal[2][1]), phi),
            "circle": chi_squared(circle(x, horizontal[4][0], horizontal[4][1]), phi)}

    # Defined variables to hold the keys which correspond to the smallest chi-squared values for each slice
    pos_diag_shape = min(positive_diag_chi_squared, key=positive_diag_chi_squared.get)
    vert_shape = min(vert, key=vert.get)
    hori_shape = min(hori, key=hori.get)

    # Checking if the smallest chi-squared value of the diagonal slice is smaller than the vertical and horizontal
    if (positive_diag_chi_squared[pos_diag_shape] < vert[vert_shape] and
            positive_diag_chi_squared[pos_diag_shape] < hori[hori_shape]):
        shape = pos_diag_shape
        # print(shape)
    # Checking if the smallest chi-squared value of the diagonal slice is greater than the vertical and horizontal

    elif (vert[vert_shape] < (positive_diag_chi_squared[pos_diag_shape]) and
          hori[hori_shape] < positive_diag_chi_squared[pos_diag_shape]):
        shape = vert_shape
        # print(shape)
    else:
        raise Exception("Chi-squared is inconclusive")

    # Intensity plot from condition of sin(phi) = 0
    plt.figure()
    plt.title(r"Slice of $sin(\phi) = 0$")
    plt.plot(x, phi, "r.", label="Data")
    plt.xlabel(r"$sin(\theta)/ rad$")
    plt.ylabel("Intensity/ $Wm^{-2}$")

    # Checking the shape variable and plotting the fitted data lines together with an asterisk to indicate
    # which fit had the lowest chi-squared value

    if shape == "square":
        plt.plot(x, square(x, vertical[0][0], vertical[0][1]), label=r"Square$^*$")
        plt.plot(x, diamond(x, vertical[2][0], vertical[2][1]), label="Diamond")
        plt.plot(x, circle(x, vertical[4][0], vertical[4][1]), label="Circle")

        intensity_x = {"square": 0, "I0": vertical[0][0], "sigma": vertical[1][0]}
        width_x = {"square": 0, "width": vertical[0][1], "sigma": vertical[1][1]}

    elif shape == "diamond":
        plt.plot(x, square(x, diag1[0][0], diag1[0][1]), label="Square")
        plt.plot(x, diamond(x, diag1[2][0], diag1[2][1]), label=r"Diamond$^*$")
        plt.plot(x, circle(x, diag1[4][0], diag1[4][1]), label="Circle")

        intensity_x = {"diamond": 0, "I0": diag1[2][0], "sigma": diag1[3][0]}
        width_x = {"diamond": 0, "width": diag1[2][1], "sigma": diag1[3][1]}

    elif shape == "circle":
        plt.plot(x, square(x, vertical[0][0], vertical[0][1]), label="Square")
        plt.plot(x, diamond(x, vertical[2][0], vertical[2][1]), label="Diamond")
        plt.plot(x, circle(x, vertical[4][0], vertical[4][1]), label=r"Circle$^*$")

        intensity_x = {"circle": 0, "I0": vertical[4][0], "sigma": vertical[5][0]}
        width_x = {"circle": 0, "width": vertical[4][1], "sigma": vertical[5][1]}

    else:
        raise Exception("Cannot determine best fit of aperture")

    plt.legend()
    plt.show()

    # Intensity plot from condition of sin(theta) = 0
    plt.figure()
    plt.title(r"Slice of $sin(\theta) = 0$")
    plt.plot(y, theta, "r.", label="Data")
    plt.xlabel(r"$sin(\phi)/ rad$")
    plt.ylabel("Intensity/ $Wm^{-2}$")

    # Checking the shape variable and plotting the fitted data lines together with an asterisk to indicate
    # which fit had the lowest chi-squared value

    if shape == "square":
        plt.plot(y, square(y, horizontal[0][0], horizontal[0][1]), label=r"Square$^*$")
        plt.plot(y, diamond(y, horizontal[2][0], horizontal[2][1]), label="Diamond")
        plt.plot(y, circle(y, horizontal[4][0], horizontal[4][1]), label="Circle")

        intensity_y = {"square": 0, "I0": horizontal[0][0], "sigma": horizontal[1][0]}
        width_y = {"square": 0, "width": horizontal[0][1], "sigma": horizontal[1][1]}

    elif shape == "diamond":
        plt.plot(y, square(y, diag1[0][0], diag1[0][1]), label="Square")
        plt.plot(y, diamond(y, diag1[2][0], diag1[2][1]), label=r"Diamond$^*$")
        plt.plot(y, circle(y, diag1[4][0], diag1[4][1]), label="Circle")

        intensity_y = {"diamond": 0, "I0": diag1[2][0], "sigma": diag1[3][0]}
        width_y = {"diamond": 0, "width": diag1[2][1], "sigma": diag1[3][1]}

    elif shape == "circle":
        plt.plot(y, square(y, horizontal[0][0], horizontal[0][1]), label="Square")
        plt.plot(y, diamond(y, horizontal[2][0], horizontal[2][1]), label="Diamond")
        plt.plot(y, circle(y, horizontal[4][0], horizontal[4][1]), label=r"Circle$^*$")

        intensity_y = {"circle": 0, "I0": horizontal[4][0], "sigma": horizontal[5][0]}
        width_y = {"circle": 0, "width": horizontal[4][1], "sigma": horizontal[5][1]}

    else:
        raise Exception("Cannot determine best fit of aperture")

    plt.legend()
    plt.show()

    return shape, intensity_x, width_x, intensity_y, width_y


def Results(shape, I0x, wx, I0y, wy):
    """
    Function that takes arguments of the shape - calculated from the chi-squared function - as well as the values and
    uncertainties of the horizontal and vertical intensities and widths
    :param shape:
    :param I0x:
    :param I0y:
    :param wx:
    :param wy:
    :return:
    """

    # Defining variables to hold the widths and associated uncertainties from the dictionaries created in
    # the previous function
    w1, w2 = np.abs(wx["width"]), np.abs(wy["width"])
    sigma1, sigma2 = wx["sigma"], wy["sigma"]

    I0A1, I0A2 = I0x["I0"], I0y["I0"]
    I0A1err, I0A2err = I0x["sigma"], I0y["sigma"]

    # Checking if the chosen aperture shape is square
    if shape == "square":
        # Checking if the width estimates are equal or not
        if w1 > w2:
            # Checking if the widths lie within at most 3 standard deviations of each other
            if w2 + 3 * sigma2 > w1:
                # Assigning the shape variable 'square'
                shape = "square"
            else:
                shape = "rectangle"
        elif w2 > w1:
            # Checking if the widths lie within at most 3 standard deviations of each other
            if w1 + 3 * sigma1 > w2:
                shape = "square"
            else:
                shape = "rectangle"
        elif w1 == 0 or w2 == 0:
            raise AttributeError("Inconclusive curve fitting returned value of 0 for one dimension")

    Area = w1 * w2
    Area_sigma = Area * np.sqrt((sigma1 / w1) ** 2 + (sigma2 / w2) ** 2)

    if I0A1 > I0A2:
        # Checking if the intensities lie within at most 3 standard deviations of each other
        if I0A2 + 3 * I0A2err > I0A1:
            # Assigning the variable E_density to the average of the two intensity values
            I0_combined = (I0A1 + I0A2) / 2
            I0_combined_sigma = np.sqrt(I0A1err ** 2 + I0A2err ** 2)
            E_density = I0_combined / Area
            E_d_sigma = E_density * np.sqrt((I0_combined_sigma / I0_combined) ** 2 + (Area_sigma / Area) ** 2)
        else:
            raise Exception("Intensities are not within 3 standard deviations of each other")
    elif I0A1 < I0A2:
        # Checking if the widths lie within at most 3 standard deviations of each other
        if I0A1 + 3 * I0A1err > I0A2:
            I0_combined = (I0A1 + I0A2) / 2
            I0_combined_sigma = np.sqrt(I0A1err ** 2 + I0A2err ** 2)
            E_density = I0_combined / Area
            E_d_sigma = E_density * np.sqrt((I0_combined_sigma / I0_combined) ** 2 + (Area_sigma / Area) ** 2)
        else:
            raise Exception("Intensities are not within 3 standard deviations of each other")
    elif I0A1 == 0 or I0A2 == 0:
        I0_combined = (I0A1 + I0A2) / 2
        I0_combined_sigma = np.sqrt(I0A1err ** 2 + I0A2err ** 2)
        E_density = I0_combined / Area
        E_d_sigma = E_density * np.sqrt((I0_combined_sigma / I0_combined) ** 2 + (Area_sigma / Area) ** 2)
        raise AttributeError("Inconclusive curve fitting returned value of 0 for one dimension")
    else:
        I0_combined = (I0A1 + I0A2) / 2
        I0_combined_sigma = np.sqrt(I0A1err ** 2 + I0A2err ** 2)
        E_density = I0_combined / Area
        E_d_sigma = E_density * np.sqrt((I0_combined_sigma / I0_combined) ** 2 + (Area_sigma / Area) ** 2)

    return shape, w1, sigma1, w2, sigma2, E_density, E_d_sigma

def ProcessData(filename):
    """
    Final code function that calls the functions below, storing all but the plotting contour function
    as a variable as this function does not return any numerical values, only plots/figures
    :param filename: 
    :return: 
    """

    # Function that returns the x,y,z data sets as well as the slice sets of sin(phi) and sin(theta) = 0
    values = data_extraction(filename)

    # Function that takes in the x, y data sets as well as the intensity data set
    plotting_contour(values[0], values[1], values[2])

    var_1 = intensity_plots(values[0], values[1], values[2], values[3], values[4])

    final = Results(var_1[0], var_1[1], var_1[2], var_1[3], var_1[4])

    results = {
        "shape": final[0],
        # one of "squares", "rectangle", "diamond", "circle" - must always be present.
        "dim_1": final[1],  # a floating point number of the first dimension expressed in microns
        "dim_1_err": final[2],  # The uncertainty in the above, also expressed in microns
        "dim_2": final[3],
        # For a rectangle, the second dimension, for other shapes, the same as dim_1
        "dim_2_err": final[4],  # The uncertainty in the above, also expressed in microns
        "I0/area": final[5],  # The fitted overall intensity value/area of the aperture.
        "I0/area_err": final[6],  # The uncertainty in the above.
    }

    return results


if __name__ == "__main__":
    # Copde placed here will not be executed when the autograder imports your function so it is a good place
    # to put code that you want to run whilst you are testing whether your ProcessData() function works correctly.
    # You can leave this code in your solution when you submit it (so long as it does not contain any SyntaxErrors!)
    # The default below should be a good first start
    from pprint import pprint  # pprint aka 'pretty-print' makes printed dictionaries look better,

    filename = "assessment_data_py22np.dat"

    test_results = ProcessData(filename)
    pprint(test_results)
