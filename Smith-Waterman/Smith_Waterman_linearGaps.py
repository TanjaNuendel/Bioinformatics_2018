##########################################################
# Bioinformatische Datenanalyse WS 18/19 - Übungsblatt 3 #
# Aufgabe 2) Locales Sequenz-Alignment                   #
# Autoren: Tanja Nündel (931179), Jonas Heinzel (931167) #
# Abgabe: 10.12.2018                                     #
##########################################################

import numpy as np
import os.path as op
from Bio import SeqIO


# Festlegen der Werte fürs hier angewandte Scoring-Modell
matVar = 1      # Variable für das Match-Scoring
misVar = -1     # Variable für das Mismatch-Scoring
gapVar = -2     # Variable für den "Gap penalty"
pt = {'match': matVar, 'mismatch': misVar, 'gap': gapVar}


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
def match(a, b):
    if a == b:
        return pt['match']
    elif a == '-' or b == '-':
        return pt['gap']
    else:
        return pt['mismatch']


# Definition der Funktion needleman(), welche zwei Arrays übergeben bekommt, welche wiederum die Sequenzen enthalten,
# welche mit Hilfe des Needleman-Wunsch-Algorithmus auf ein Globales Alignment hin überprüft werden sollen.
def smithwaterman(seq1, seq2, count):
    len1 = len(seq1)    # Referenz-Sequenz (In der Matrix vertikal)
    len2 = len(seq2)    # Read-Sequenz (In der Matrix horizontal)
    readscore = []

    print('\nLambda-Referenz: \t\tLänge: ' + str(len1))
    print('\t' + str(seq1))
    print('\nRead ' + str(count) + '\t\t\t\t\tLänge: ' + str(len2))
    print('\t' + str(seq2))

    # Erstellen der Sequenzmatrix, zunächst gefüllt mit Nullen, mit zusätzlicher Null-Spalte und -Zeile.
    score = np.zeros((len1 + 1, len2 + 1), dtype=int)

    # Spaltenweises Füllen der Sequenzmatrix mit dem maximalen Scoring für jede Position
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            diag = score[i-1][j-1] + match(seq1[i-1], seq2[j-1])
            delete = score[i-1][j] + pt['gap']
            insert = score[i][j-1] + pt['gap']
            score[i][j] = max(diag, delete, insert, 0)

    # Ermittlung der Maximalwerte der Scoring-Matrix (Alle Alignment-Openings) um die lokalen Alignments zu finden

    # print("\nScore-Matrix: (D)")
    # print('\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in score]))

    maxVal = np.amax(score)
    maximumValues = []

    print('\n\tMaximalwert der Sequenzmatrix: ' + str(maxVal))

    for i in range(0, len1+1):
        for j in range(0, len2+1):
            if score[i][j] == maxVal:
                maxi = i
                maxj = j
                appender = [maxi, maxj]
                maximumValues.append(appender)
                print('\n\t##\tMaximum %d found at: ' % maxVal)
                print('\t\t\ti: ' + str(i) + '\tj: ' + str(j))


    # Ab hier: Backtracking - Ermittlung des Alignments von jedem Max aus in Richtung "oben-links" in der Sequenzmatrix

    for maxValueArray in maximumValues:
        countDeletion = 0
        countInsertion = 0
        align1 = ''
        align2 = ''
        i = maxValueArray[0]
        j = maxValueArray[1]

        # Fallunterscheidung für das Backtracking in der Sequenzmatrix, je nach Score, ob beim Backtracking diagonal,
        # nach oben, oder nach links gegangen wird (und das für die Fälle: (1) innerhalb der Matrix (i und j > 0),
        # (2) am linken Rand der Matrix (nur i>0) , (3) am oberen Rand der Matrix (nur j>0).
        # Zusammenbau der Alignments während des Backtrackings.

        # Backtracking:
        while i > 0 and j > 0:
            score_current = score[i][j]
            score_diag = score[i-1][j-1]
            score_left = score[i][j-1]
            score_up = score[i-1][j]

            # Sprung diagonal
            if score_current == score_diag + match(seq1[i-1], seq2[j-1]):
                align1 += seq1[i-1]
                align2 += seq2[j-1]
                i, j = i-1, j-1
            # Sprung nach oben
            elif score_current == score_up + pt['gap']:
                align1 += seq1[i-1]
                align2 += '-'
                countDeletion += 1
                i -= 1
            # Sprung nach links
            elif score_current == score_left + pt['gap']:
                align1 += '-'
                align2 += seq2[j-1]
                countInsertion += 1
                j -= 1
            else:
                i, j = 0, 0

        # Umdrehen der Alignment-Strings, da sie in umgekehrter Richtung zusammengesetzt wurden.
        align1 = align1[::-1]
        align2 = align2[::-1]
        sequencelength = len(align1)
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
            print('\tScoring: \t\t\t%d\t\t(Match: %d ; Mismatch: %d ; Gap: %d)' % (seq_score, matVar, misVar, gapVar))
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
        scorings.append(smithwaterman(reference, singleread, count+1))
    print('\n\t## Scorings insgesamt:')
    for index, reads  in enumerate(scorings):
        print('\n\t## Read ' + str(index+1) + ' hat Scorings: ')
        for index2, scorings in enumerate(reads):
            print('\t\t## Alignment ' + str(index2 + 1) + ' hat Scoring: ' + str(scorings))
        bestAlignment.append(np.amax(reads))
    print('\n\t## Bester Read: ' + str(bestAlignment.index(np.amax(bestAlignment))+1))
    print('\t\t## Bestes Scoring ' + str(np.amax(bestAlignment)))


# Einlesen der gegebenen .fasta Sequenz-Datensätze. Zu Testzwecken wurden jedoch kürzere Sequenzen verwendet!
# Bitte auf die Dateinamen achten!

lambdaSeq = readfasta("venv/data/lambda_ref_test.fasta")
reads = readfasta("venv/data/reads.fasta")
aligner(lambdaSeq, reads)


# Testsequenzen, da inbesondere die gegebene Referenzsequenz mit 48.000 Basen sehr lange zur Iteration benötigt.

# seq1 = [["G", "G", "G", "C", "G", "G", "C", "G", "A", "C", "C", "T", "C", "G", "C", "G", "G", "G", "T", "T"]]
# seq2 = [["T", "G", "A", "A", "C", "T", "T", "G", "C", "T", "C", "T", "C", "T", "G", "A", "T", "T"]]
#         ["G", "C", "G", "G", "C", "G", "T", "A", "G", "C", "T", "G", "C", "G", "G", "G", "T", "T", "A"]]

# Test der Funktionalität mit Hilfe der Testsequenzen:

# aligner(seq1, seq2)