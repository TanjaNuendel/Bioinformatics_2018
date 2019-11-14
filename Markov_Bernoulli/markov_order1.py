##########################################################
# Bioinformatische Datenanalyse WS 18/19 - Übungsblatt 4 #
# Aufgabe 2.b/c) Markov-Kette 1. Ordnung                     #
# Autoren: Jonas Heinzel (931167), Tanja Nündel (931179) #
# Abgabe: 07.01.2019                                     #
##########################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


# Dictionary mit Arrays, welche die Wahrscheinlichkeiten des zukünftigen Symbols anhand des Vorgängers aufschlüsseln
# siehe dazu Tabelle aus Aufg. 2.a)
# Prinzip: Startwert A, Übergang zur Base [A, C, G, T] mit den Werten [0.6, 0.1, 0.1, 0.2]
state_dict = {'A': np.array([['A', 'C', 'G', 'T'], [0.6, 0.1, 0.1, 0.2]]),
              'C': np.array([['A', 'C', 'G', 'T'], [0.1, 0.5, 0.3, 0.1]]),
              'G': np.array([['A', 'C', 'G', 'T'], [0.05, 0.2, 0.7, 0.05]]),
              'T': np.array([['A', 'C', 'G', 'T'], [0.4, 0.05, 0.05, 0.5]])
              }


# Klasse für das Markov-Objekt und seine Zustände
class Markov(object):

    def __init__(self, state_dict):
        self.state_dict = state_dict
        self.state = list(self.state_dict.keys())[0]

    # Überprüfung, in welchem Zustand die Kette aktuell ist - kann optional eingebaut/genutzt werden
    def check_state(self):
        print("\nAktueller Basenzustand: %s" % self.state)

    # Setzen und Ausgabe des Startzustands
    def set_state(self, state):
        self.state = state
        print("\nErster Basenzustand wird gesetzt auf: %s\t" % self.state)

    # Setzen des nächstes Zustandes, basierend auf übergebenen Wahrscheinlichkeiten (state_dict)
    def next_state(self):
        A = self.state_dict[self.state]
        # Hier wird der nächste Status aus dem Dictionary randomisiert gebildet
        # mit a = key und p = value aus dem übergebenen dict
        # mit [0] 1. Arraydimension (Base), [1] 2. Arraydimension (Wahrscheinlichkeit)
        self.state = np.random.choice(a=list(A[0]), p=list(A[1]))
        return self.state


# Kombiniert alle Einträge des dictionaries miteinander und verbindet sie zu Paaren
# Hierbei wird aufgeschlüsselt in die Unterschiede bspw. der Folgen AT und TA
# Beispiel: {A:2, B:1} => {AA:4, AB:2, BA:2, BB:1}
def generate2merdict(state_dict):
    gen2merdict= {}
    print("\nÜbergangswahrscheinlichkeiten der Markovkette 1. Ordnung: ")
    for key, val in iter(state_dict.items()):
        secKeyarray = val[0]
        valarray = val[1]
        print("\tFür Base " + str(key) + "\t-\tÜbergangswerte zu A, C, G, T: " + str(valarray))
        for index, sk in enumerate(secKeyarray):
            gen2merdict["".join((key, sk))] = valarray[index]
    return gen2merdict


# Generierungsfunktion für die Sequenzen - erhält die gewünschte Länge, die in der markov1() festgelegt wird
def generatebasearray(length):
    basearray = []

    # Instanziierung des Markov-Objekts
    diagram_a = Markov(state_dict)

    # Startstatusfestlegung über die np.random.choice der next_state() des Markov-Objekts
    startstate = diagram_a.next_state()
    diagram_a.set_state(startstate)

    # Berechnung des next State über die übergebene Länge
    # Übergabe des jeweiligen Next State an ein Sequenzarray (exklusive Startbase)
    for x in range(length):
        diagram_a.next_state()
        basearray.append(diagram_a.state)
    # print("Basearray aus der generate(): " + str(basearray))
    return basearray


# Zusammenstellung aller Key in einem Array
def getallkeys(dictionary):
    allkeys = []
    for k, v in iter(dictionary.items()):
        allkeys.append(k)
    return allkeys


# Zusammenstellung aller Values in einem Array
def getallvalues(dictionary):
    allvalues = []
    for k,v in iter(dictionary.items()):
        allvalues.append(v)
    return allvalues


# Für die Achsenbeschriftung
def plotformat(x, pos):
    'The two args are the value and tick position'
    return '%1.1f' % (x)


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
    print("\nBasenhäufigkeiten:")
    for key, value in iter(countdict.items()):
        print("\tBasen " + str(key) + ": " + str(value))
        basearray.append(value)
    return basearray


# Generiert ein Array der 2mere
def generatetwomeresarray(onemeresarray):
    twomeresarray = []
    i = 0
    while i<=len(onemeresarray):
        if(i<=len(onemeresarray)-1):
            twomeresarray.append(onemeresarray[i] + onemeresarray[i+1])
        i += 2
    return twomeresarray


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


# Allgemeiner Programmablauf
def markov1():

    # generiert das Array mit vorgegebener Länge
    sequencearray = generatebasearray(length=1000)
    print("\nGenerierte Bernoulli-Sequenz: " + "\n\t" + str(sequencearray))

    # generiert das Zustandsübergangs-Dictionary und zieht die keys in ein Array
    gen2merdict = generate2merdict(state_dict)
    allkeysarray = getallkeys(gen2merdict)
    print("\nWahrscheinlichkeiten der 2-mere: " + "\n\t" + str(gen2merdict))
    print("\nAlle möglichen 2-mere: " + "\n\t" + str(allkeysarray))

    # Absolute Häufigkeitsberechnung der einfachen Basen
    simplekeys = ["A", "C", "G", "T"]
    base1count = basecounter(sequencearray, simplekeys)

    # generiert ein 2-mer-Array aus der Bernoulli-Sequenz und berechnet die absolute Häufigkeit dieser
    twomeresarray = generatetwomeresarray(sequencearray)
    base2count = basecounter(twomeresarray, allkeysarray)

    # Aufruf der Plotfunktionen für die absoluten Häufigkeiten der einfachen Basen und 2 mere
    hairyplotter(base1count, simplekeys)
    hairyplotter(base2count, allkeysarray)


# ausführen der markov1() zum Starten der Erzeugung einer Markov-Kette 1. Ordnung
markov1()