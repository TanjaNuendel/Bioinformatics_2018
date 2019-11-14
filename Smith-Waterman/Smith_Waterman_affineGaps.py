##########################################################
# Bioinformatische Datenanalyse WS 18/19 - Übungsblatt 3 #
# Aufgabe 3) Affine Gapkosten                            #
# Autoren: Tanja Nündel (931179), Jonas Heinzel (931167) #
# Abgabe: 13.01.2019                                     #
##########################################################

import numpy as np
import os.path as op
from Bio import SeqIO


# Festlegen der Werte fürs hier angewandte Scoring-Modell
matVar = 1      # Variable für das Match-Scoring
misVar = -1     # Variable für das Mismatch-Scoring
gapOpVar = -5   # Variable für den "Gap opening penalty"
gapExVar = -1   # Variable für den "Gap extension penalty"
pt = {'match': matVar, 'mismatch': misVar, 'gapEx': gapExVar, 'gapOp': gapOpVar}

global MIN
MIN = -float("inf")


# Funktion checkfasta zur Überprüfung ob der übergebene Pfad eine gültige .fasta -Datei enthält
def checkfasta(path):
    if op.isfile(path) and (path[-6:] == ".fasta"):
        print("\t## Reading .fasta file successful.")
        return True
    else:
        print("\tERROR: no valid file format or path.")
        return False


# Definition der Funktion readfasta(), welche zum Einlesen der .fasta -Dateien einen Pfad als String übergeben bekommt
# Sie gibt ein zweidimensionales Array zurück [[Sequenz1], [Sequenz2], ... , [SequenzN]], wenn der übergebene String
# eine gültige .fasta -Datei enthält.
def readfasta(path):
    zielarray = []
    sequenz = 0
    if checkfasta(path):
        with open(path, "r") as handle:
            # Iteration über alle Sequenzen in der eingelesenen fasta-Datei. Diese werden ins Zielarray geschrieben
            for seq_record in SeqIO.parse(handle, "fasta"):
                templiste = list(seq_record.seq)
                seq_array = np.array(templiste)
                if zielarray == []:
                    zielarray = [seq_array]
                    print("\tSequenz %i eingelesen" % len(zielarray))
                else:
                    zielarray.append(seq_array)
                    print("\tSequenz %i eingelesen" % len(zielarray))
                sequenz += 1
    else:
        # Errorfall, wenn der Dateityp nicht korrekt ist. Die Funktion gibt in diesem Fall ein leeres Array zurück.
        print("\t## ERROR: no valid data type given")
    return zielarray


# Definition der Funktion match(), welche zwei übergebene Werte a und b auf Gleichheit (match), Sequenzlücken (gaps),
# bzw Ungleichheit (mismatch) überprüft und den zugehörigen Scoring-Wert hierfür zurückgibt.
def match(a, b):  # a = Spalten (vertikal), b = Zeilen (horizontal)
    if a == b:
        return pt['match']
    else:
        return pt['mismatch']


# Funktion, welche die Matrix für das Backtracking initialisiert und mit leeren Strings füllt
def init_backtracking_matrix(i,j):
    if i >= 0 and j >= 0:
        return "   "


# Funktion, welche die Matrix für Deletions initialisiert und mit entsprechenden Werten füllt
def init_deletion_matrix(i, j):
    if i > 0 and j == 0:
        return MIN
    elif j > 0 and i == 0:
        return MIN
    else:
        return 0


# Funktion, welche die Matrix für Insertions initialisiert und mit entsprechenden Werten füllt
def init_insertion_matrix(i, j):
    if j > 0 and i == 0:
        return MIN
    elif i > 0 and j == 0:
        return MIN
    else:
        return 0


# Funktion, welche die Scoring-Matrix initialisiert / mit Nullen füllt
def init_matrix(i, j):
    return 0


# Definition der Funktion smithwaterman()
def smithwatermanaffine(seq1, seq2, count):
    len1 = len(seq1)    # Referenz-Sequenz (In der Matrix vertikal)
    len2 = len(seq2)    # Read-Sequenz (In der Matrix horizontal)
    readscore = []

    # Initialisierung der 4 Matrizen
    backtrackScore = [[init_backtracking_matrix(i, j) for i in range(0, len2+1)] for j in range(0, len1+1)]  #BS
    scoreArr = [[init_matrix(i, j) for i in range(0, len2 + 1)] for j in range(0, len1 + 1)]  # M
    delArr = [[init_deletion_matrix(i, j) for i in range(0, len2 + 1)] for j in range(0, len1 + 1)]  # y
    insArr = [[init_insertion_matrix(i, j) for i in range(0, len2 + 1)] for j in range(0, len1 + 1)]  # x

    print('\nLambda-Referenz: \t\tLänge: ' + str(len1))
    print('\t' + str(seq1))
    print('\nRead ' + str(count) + '\t\t\t\t\tLänge: ' + str(len2))
    print('\t' + str(seq2))

    # Zeilenweises Füllen der 3 Matrizen mit den Regeln nach Gotoh (lokale Alignments mit affinen Gapkosten)
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):

            delPreTop = delArr[i-1][j]
            insPreLeft = insArr[i][j-1]
            scorePreLeft = scoreArr[i][j-1]
            scorePreTop = scoreArr[i-1][j]
            scorePreDiag = scoreArr[i-1][j-1]

            delArr[i][j] = max(scorePreTop + pt['gapOp'] + pt['gapEx'],  #
                               delPreTop + pt['gapEx'])    # (P)

            insArr[i][j] = max(scorePreLeft + pt['gapOp'] + pt['gapEx'],  #
                               insPreLeft + pt['gapEx'])    # (Q)

            scoreArr[i][j] = max(scorePreDiag + match(seq1[i-1], seq2[j-1]),  # match x(i), y(j)
                                 delArr[i][j],
                                 insArr[i][j],
                                 0)

            # Füllen der backtrackScore-Matrix, welche wir für das Backtracking verwenden. (wie Pfeile in einer
            # Scoring-Matrix die auf das Element zeigen, aus welchem der aktuelle Wert errechnet wurde)
            if scoreArr[i][j] == scorePreDiag + match(seq1[i-1], seq2[j-1]):
                backtrackScore[i][j] = "mat"
            elif scoreArr[i][j] == scorePreLeft + pt['gapOp'] + pt['gapEx']:
                backtrackScore[i][j] = "ins"
            elif scoreArr[i][j] == scorePreTop + pt['gapOp'] + pt['gapEx']:
                backtrackScore[i][j] = "del"

    # Test-Plots für die Matrizen

    # print("\nScore-Matrix: (D)")
    # print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in scoreArr]))

    # print("\nDeletion-Matrix: (P)")
    # print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in delArr]))

    # print("\nInsertion-Matrix: (Q)")
    # print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in insArr]))

    # print("\nBacktracking-Matrix (Score): ")
    # print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in backtrackScore]))

    # Ermittlung der Maximalwerte der Scoring-Matrix (Alle Alignment-Openings) um die lokalen Alignments zu finden

    maxVal = np.amax(scoreArr)
    maximumValues = []

    print('\n\tMaximalwert der Score-Matrix: ' + str(maxVal))

    for i in range(0, len1+1):
        for j in range(0, len2+1):
            if scoreArr[i][j] == maxVal:
                maxi = i
                maxj = j
                appender = [maxi, maxj]
                maximumValues.append(appender)
                print('\n\t##\tMaximum %d found at: ' % maxVal)
                print('\t\t\ti: ' + str(i) + '\tj: ' + str(j))

    # Ab hier: Backtracking - Ermittlung des Alignments von "unten rechts" nach "oben links" durch die Sequenz-Matrix.

    for maxValueArray in maximumValues:
        countDeletion = 0
        countInsertion = 0
        align1 = ''
        align2 = ''
        i = maxValueArray[0]
        j = maxValueArray[1]

        # Fallunterscheidung für das Backtracking in der Sequenzmatrix, je nach Wert in der Backtrack - Matrix.
        # Abhängig von match / insertion / deletion werden hier die Alignment-Strings zusammengebaut.
        while i > 0 and j > 0:
            # Sprung diagonal -> Match oder Mismatch
            if backtrackScore[i][j] == "mat":
                align1 += seq1[i-1]
                align2 += seq2[j-1]
                i -= 1
                j -= 1
            # Sprung nach oben -> Insertion im Read
            elif backtrackScore[i][j] == "del":
                align1 += seq1[i-1]
                align2 += '-'
                countInsertion += 1
                i -= 1
            # Sprung nach links -> Deletion im Read
            elif backtrackScore[i][j] == "ins":
                align1 += '-'
                align2 += seq2[j-1]
                countDeletion += 1
                j -= 1
            else:
                i, j = 0, 0

        # Umdrehen der Alignment-Strings, da sie in umgekehrter Richtung zusammengesetzt wurden.
        align1 = align1[::-1]
        align2 = align2[::-1]
        sequencelength = len(align1)  # Länge der Referenz-Sequenz
        sym = ''
        seq_score = 0
        matches = 0

        for i in range(sequencelength):
            a1 = align1[i]
            a2 = align2[i]
            if a1 == a2:
                sym += "|"
                matches += 1
                seq_score += match(a1, a2)
            else:
                seq_score += match(a1, a2)
                sym += ' '

        # Ausgabe, wenn es ein sinnvolles Alignment gibt.
        if sequencelength > 0 and matches > 0:
            prozent = matches / sequencelength * 100
            print("\n\tReferenz:  " + str(align1))
            print("\tSymmetrie: " + str(sym))
            print("\tRead-Seq:  " + str(align2) + "\n")
            print('\tSequence Length: \t' + str(sequencelength))
            print('\tÜbereinstimmung: \t%2.1f Prozent' % prozent)
            print('\tScoring: \t\t\t%d\t\t(Match: %d ; Mismatch: %d ; Gap-Op: %d ; Gap-Ex: %d)' %
                  (seq_score, matVar, misVar, gapOpVar, gapExVar))
            print('\tMatches: \t\t\t%d' % matches)
            print('\tMisMatches: \t\t%d' % (sequencelength - matches))
            print('\tInsertions: \t\t%d' % countInsertion)
            print('\tDeletions: \t\t\t%d' % countDeletion)
        else:
            print('\n## Kein Alignment gefunden. ##\n')
        print('--------------------------')
        readscore.append(seq_score)
        # readscore gibt die Scorings der einzelnen Alignments in einem Array zurück [[score 1],[score 2],...,[score n]]
    return readscore


# Definition der Funktion aligner(), welche eine referenz-Sequenz und Reads jeweils in einem Array entgegennimmt und
# diese über die Funktion needleman() auf Globale Alignments überprüft.
def aligner(ref, reads):
    reference = []
    scorings = []
    bestAlignment = []
    for lambdaref in ref:
        reference = list(lambdaref)
    for count, sequence in enumerate(reads):
        singleread = list(sequence)
        scorings.append(smithwatermanaffine(reference, singleread, count+1))

    print('\n\t## Scorings insgesamt:')
    for index, reads in enumerate(scorings):
        if reads:
            print('\n\t## Read ' + str(index+1) + ' hat Scorings: ')
            for index2, scorings in enumerate(reads):
                print('\t\t## Alignment ' + str(index2 + 1) + ' hat Scoring: ' + str(scorings))
            bestAlignment.append(np.amax(reads))

    if bestAlignment:
        print('\n\t## Bester Read: ' + str(bestAlignment.index(np.amax(bestAlignment))+1))
        print('\t\t## Bestes Scoring ' + str(np.amax(bestAlignment)))
    else:
        print('\n\t## Keine Alignments gefunden')


# Einlesen der gegebenen .fasta Sequenz-Datensätze. Zu Testzwecken wurden jedoch kürzere Sequenzen verwendet!
# Bitte auf die Dateinamen achten!

lambdaSeq = readfasta("venv/data/lambda_ref_test.fasta")
reads = readfasta("venv/data/reads.fasta")
aligner(lambdaSeq, reads)


# Testsequenzen, da inbesondere die gegebene Referenzsequenz mit 48.000 Basen sehr lange zur Iteration benötigt.

# seq1 = [["G", "G", "G", "C", "G", "G", "C", "G", "A", "C", "C", "T", "C", "G", "C", "G", "G", "G", "T", "T"]]
# seq2 = [["T", "G", "A", "A", "C", "T", "T", "G", "C", "T", "C", "T", "C", "T", "G", "A", "T", "T"],
#         ["G", "C", "G", "G", "C", "G", "T", "A", "G", "C", "T", "G", "C", "G", "G", "G", "T", "T", "A"]]

# aligner(seq1, seq2)