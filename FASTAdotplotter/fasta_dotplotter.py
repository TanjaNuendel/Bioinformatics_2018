# Bioinformatische Datenanalyse WS 18/19 - Übungsblatt 1
# Autoren: Tanja Nündel (931179), Jonas Heinzel (931167)
# Abgabe: 15.10.2018

import numpy as np
from Bio import SeqIO
from matplotlib import pyplot as plt

bio_arr1 = []
bio_arr2 = []


# Festlegen der Fenstergröße und Übereinstimmungszahl:

window = 2
threshold = 1


# Einlesen der Sequenzdaten aus der ersten .fasta-Datei in ein Array mit Offset 1
# Bitte die relativen Pfade der .fasta -Dateien beachten!

with open("../data/gi7305369.fasta", "rU") as handle:
    for seq_record in SeqIO.parse(handle, "fasta"):
        templiste = list(seq_record.seq)
        bio_arr1 = np.array(templiste)
        print("Laenge Array 1: " + str(len(bio_arr1)))

        print("Sequenz-Array 1: " + str(bio_arr1))

        # Einteilung des Arrays in Fenster erfolgt bereits während des Parsens,
        # um diese im weiteren Verlauf einfacher vergleichen zu können

        bio_arr1_long = []
        for item in range(len(bio_arr1)-window):
            if bio_arr1[item]:
                for element in range(window):

                    # die .join() -Methode erzeugt einen String aus mehreren Einzelstrings, wonach eine zusammen-
                    # hängende Sequenz in der gewünschten Fensterlänge aus einzelnen Aminosäuren erstellt wird

                    substring = ''.join(str(x) for x in bio_arr1[item:item+int(window)])
                bio_arr1_long.append(substring)
            else:
                break


# Einlesen der Sequenzdaten aus der zweiten .fasta-Datei in ein Array, ebenfalls mit Offset 1

with open("../data/gi12643549.fasta", "rU") as handle:
    for seq_record in SeqIO.parse(handle, "fasta"):
        templiste = list(seq_record.seq)
        bio_arr2 = np.array(templiste)
        print("Laenge Array 2: " + str(len(bio_arr2)))

        print("Sequenz-Array 2: " + str(bio_arr2))

        # Einteilung des Arrays in Fenster erfolgt bereits während des Parsens,
        # um diese im weiteren Verlauf einfacher vergleichen zu können

        bio_arr2_long = []
        for item in range(len(bio_arr2)-window):
            if bio_arr2[item]:
                for element in range(window):
                    # .join() -Methode siehe oben.
                    substring = ''.join(str(x) for x in bio_arr2[item:item+int(window)])
                bio_arr2_long.append(substring)
            else:
                break


# Funktion match zur Ermittlung von Matches zwischen den beiden übergebenen Sequenzen
# Übergeben werden ein 2-dimensionales Array, welches am Index 0 das Array der ersten Sequenz und am Index 1 das
# Array der zweiten Sequenz enthält. Zurückgegeben wird ein 2D-Array, welches die Koordinaten der Matches enthält.

def match(givenarray, thrsh):

    matcharray1 = []
    matcharray2 = []

    # 1./2. for-Schleife iterieren über alle Elemente
    # matcharray.append fügt bei einem Match den Zielarrays mit dem Index des Matches die entsprechende Koordinate hinzu

    for indx, x in enumerate(givenarray[0]):
        for indy, y in enumerate(givenarray[1]):
            if x == y:
                matcharray1.append(indx)
                matcharray2.append(indy)

            # bei mismatch werden die Aminosäuren im aktuellen Fenster einzeln auf Übereinstimmung verglichen
            # und in der Summe gegen den Threshold-Wert aufgerechnet.
            else:
                thrsh_counter = 0
                for a, b in zip(str(x), str(y)):
                    if a == b:
                        thrsh_counter += 1
                if thrsh_counter / window >= thrsh:
                    matcharray1.append(indx)
                    matcharray2.append(indy)
    return [matcharray1, matcharray2]


# Funktion hairyplotter nimmt das Match-Array, das Merge-Array mit den beiden Sequenzen, die Fenstergröße und die
# Übereinstimmungszahl entgegen. Das Match-Array wird für die Darstellung der Matchings, das Merge-Array für die
# Länge der Achsen und Fenstergröße sowie Übereinstimmungszahl lediglich für die Darstellung im Plot benötigt.

def hairyplotter(matcharray, superarray, window, threshold):

    endx = len(superarray[0])
    endy = len(superarray[1])
    fig = plt.figure()
    ax= fig.add_subplot(111)
    ax.plot(matcharray[0], matcharray[1], 'k.')
    ax.set_aspect(endx/endy)
    plt.axis([0, endx, 0, endy])
    plt.gca().invert_yaxis()
    plt.gca().xaxis.set_ticks_position('top')
    plt.xlabel('Maus PAX-6', fontsize=16)
    plt.gca().xaxis.set_label_position('top')
    plt.ylabel('Drosophila  ' + '$\mathit{eyeless}$', fontsize=16)
    plt.title('window = %i' %window, loc='left')
    plt.title('threshold = %.1f' % threshold, loc='right')
    plt.show()


# Aufrufe zur Überprüfung der Arrays und der Funktionalität in der run-Konsole:

print("Erstes sequenziertes Array: " + str(bio_arr1_long))
print("Zweites sequenziertes Array: " + str(bio_arr2_long))

mergearray = [bio_arr1_long, bio_arr2_long]
print("Merge-Array: " + str(mergearray))

matcharray = match(mergearray, threshold)
print("Match-Array " + str(matcharray))


# Aufruf der Plotfunktion

hairyplotter(matcharray, mergearray, window, threshold)
