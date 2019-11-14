##########################################################
# Bioinformatische Datenanalyse WS 18/19 - Übungsblatt 4 #
# Aufgabe 1.b/c) Markov-Kette 0. Ordnung                 #
# Autoren: Jonas Heinzel (931167), Tanja Nündel (931179) #
# Abgabe: 07.01.2019                                     #
##########################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Erstellt ein Dictionary aus Wahrscheinlichkeitswerten
def dictgenerator(keys, values):
    probabilities = {}
    if len(keys) == len(values):
        for index, key in enumerate(keys):
            probabilities[key] = values[index]
    return probabilities


# Kombiniert alle Einträge des Dictionaries miteinander und multipliziert deren Schlüssel
# Beispiel: {A:2, B:1} => {AA:4, AB:2, BA:2, BB:1}
def dictsquare(dictionary):
    probabilities = {}
    for first in dictionary:
        for second in dictionary:
            probabilities["".join((first, second))] = dictionary[first] * dictionary[second]
    return probabilities


# Bekommt ein festes Dictionary (mit Schlüsseln und zugehörigen Wahrscheinlichkeits-Werten) und eine gewünschte
# Sequenzlänge. Generiert daraus ein Array entsprechender Länge
def generatebasearray(dict, length):
    allkeys = []
    allvalues = []
    for k, v in iter(dict.items()):
        allkeys.append(k)
        allvalues.append(v)
    # Hier: Erzeugung der Bernoulli-Sequenz mit Übergabe der dict-Schlüssel, der gewünschten Länge und der dict-Werte
    return np.random.choice(allkeys, length, True, allvalues)


# Für die Achsenbeschriftung
def plotformat(x, pos):
    'The two args are the value and tick position'
    return '%1.1f' % (x)


# Funktion zum Plotten der Wahrscheinlichkeitsverteilung
def hairyplotter(array, elements):
    formatter = FuncFormatter(plotformat)
    x = np.arange(len(elements))
    fig, ax = plt.subplots()
    ax.yaxis.set_major_formatter(formatter)
    plt.bar(x, array)
    alphabet = elements
    axisdescription = []
    for times in range(0, int(round(len(elements) / len(alphabet)))):
        for chars in alphabet:
            axisdescription.append(chars)
    plt.xticks(x, axisdescription)  # ('A', 'C', 'G', 'T')
    plt.show()


# Funktion, welche die jeweiligen Basen eines Arrays zählt, die ebenso im zweiten, übergebenen Array zu finden sind.
# (sozusagen Zählen nach vorgegebenem Alphabet) - nutzt zur key-value-Zuordnung ein Dict (countdict)
# enthält Ausgabefunktion der absoluten Basenhäufigkeiten und das dazugehörige Array für den Plotter
def basecounter(array, combinations):
    countdict = {}
    basearray = []
    for checkbase in combinations:
        countdict[checkbase] = 0
    for index, key in enumerate(array):
        if key in combinations:
            countdict[key] += 1
    for key, value in iter(countdict.items()):
        print("\tBasen " + str(key) + ": " + str(value))
        basearray.append(value)
    return basearray


# Erstellt ein dictionary aus zwei Listen mit je Basen und deren Einzelwahrscheinlichkeiten
# Kombiniert per dictsquare alle Einträge und multipliziert deren Schlüssel
def getallcombinations(keys, values):
    dict = dictgenerator(keys, values)
    return dictsquare(dict)


# Bekommt ein dictionary übergeben und gibt ein Array zurück, das alle Schlüssel des dicts enthält.
def getallkeys(dict):
    allkeys = []
    for k, v in iter(dict.items()):
        allkeys.append(k)
    return allkeys


# Allgemeiner Programmablauf
# Die einfachen Kombis und
def markov0():

    # Wahrscheinlichkeitsberechnung der Einzelbasen
    base = dictgenerator(("A", "C", "G", "T"), (0.4, 0.1, 0.1, 0.4))
    print("\nEinzelwahrscheinlichkeiten der Basen: " + "\n\t" + str(base))

    # Generierung des Arrays für die Bernoulli-Sequenz mit einer Länge von 1000 Basen
    sequencearray = generatebasearray(base, 1000)
    print("\nGenerierte Sequenz: " + "\n\t" + str(sequencearray))

    # Wahrscheinlichkeitsberechnung der 2-mere (für eine separat erzeugte Sequenz mit einer Länge von 1000 Basen)
    bases = getallcombinations(("A", "C", "G", "T"), (0.4, 0.1, 0.1, 0.4))
    combinationsarray = generatebasearray(bases, 1000)
    print("\nErrechnete Wahrscheinlichkeiten der möglichen 2-mere : " + "\n\t" + str(bases))
    print("\t(Format ist intern bedingt durch Float-Ungenauigkeiten)")
    print("\nErzeugte 2-mere (seperate Sequenz): " + "\n\t" + str(combinationsarray))

    # Häufigkeiten
    #print("\ntestarray: " + str(testarray))
    print("\nAbsolute Häufigkeiten der Basen und Basenkombinationen der generierten 2-mere-Sequenz: ")
    basearray = basecounter(combinationsarray, getallkeys(bases))


    # Aufruf der Plotfunktionen
    hairyplotter(basearray, getallkeys(bases))
    hairyplotter(basecounter(generatebasearray(base, 1000), getallkeys(base)), getallkeys(base))


# ausführen der markov0() zum Starten der Erzeugung einer Markov-Kette 0. Ordnung
markov0()
