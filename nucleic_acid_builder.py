#Andy Jiang and Asim Alenaizan, Sherrill Group, Georgia Institute of Technology

import sys

from openbabel import openbabel as ob

import numpy as np
from numpy import linalg
import math

import base_step_final 

#This function allows the calculation of the direction vectors (x, y, and z) of a given molecule (OBMol object)
def getOriginAndVectors():

    origin = np.array([0, 0, 0])
    x_vector = np.array([1, 0, 0])
    y_vector = np.array([0, 1, 0])
    z_vector = np.array([0, 0, 1])

    return (origin, x_vector, y_vector, z_vector)

#Changes the coordinates of a molecule based on its direction vectors and origin
def getNewMolecule(mol, newCoords):
    """
    mol = the nucleotide
    newCoords = the coordinates of the nucleotide after the transformation
    """
    oldCoords = getOriginAndVectors()

    for atom in ob.OBMolAtomIter(mol):

        atomCoords = np.array([atom.GetX(), atom.GetY(), atom.GetZ()], dtype=float) - oldCoords[0]

        oldTransMatrixInv = np.array([oldCoords[1], oldCoords[2], oldCoords[3]], dtype=float)
        newTransMatrix = np.transpose(np.array([newCoords[1], newCoords[2], newCoords[3]], dtype=float))

        transMatrix = np.matmul(newTransMatrix, oldTransMatrixInv)
        atomCoords = np.matmul(transMatrix, atomCoords) + newCoords[0]

        atom.SetVector(atomCoords[0], atomCoords[1], atomCoords[2])

    return mol

#Returns the number of rings in a molecule, useful for checking if a nucleobase is a pyrimidine or a purine
def GetNumRings(mol):
    count = 0

    for ring in ob.OBMolRingIter(mol):
        count += 1

    return count

#Prints the name of each atom in a molecule, useful for debugging purposes
def printAtomNames(mol):
    for res in ob.OBResidueIter(mol):
        for atom in ob.OBResidueAtomIter(res):
            name = res.GetAtomID(atom)
            print(name)

#This function returns a nucleobase after the translations and rotations have been applied, also sets the chain of the nucleobase
#i = the order of the nucleobase in the stack
#newMol = the nucleobase to be added
def build(i, mol, newCoords):

    mol = getNewMolecule(mol, newCoords)

    mol.GetResidue(0).SetChain(chr(ord('A') + i - 1))

    mol.SetChainsPerceived()

    return mol

#Returns the name of a nucleobase file based on its letter (A, U, T, G, C)
def letterToNucleotide(letter):
    if letter == 'A':
        return "adenine.pdb"
    elif letter == 'U':
        return "uracil.pdb"
    elif letter == 'T':
        return "thymine.pdb"
    elif letter == 'G':
        return 'guanine.pdb'
    elif letter == 'C':
        return "cytosine.pdb"

#This function stacks nucleotides based on inputs from a file (stackFile is the file name)
def stack(stackFile):

    conv = ob.OBConversion()

    ladder = ob.OBMol()

    stackFile = open(stackFile, "r")
    stackLines = stackFile.readlines()
    stackFile.close()

    #ladder is the stack of nucleobases that will be returned by this function
    ladder = ob.OBMol()

    baseStepVectors = []

    for i in range(1, len(stackLines)):
        el = stackLines[i].split()
        mol = ob.OBMol()
        conv.ReadFile(mol, letterToNucleotide(el[0]))

        Dx = float(el[1])
        Dy = float(el[2])
        Dz = float(el[3])
        Rx = math.radians(float(el[4]))
        Ry = math.radians(float(el[5]))
        Rz = math.radians(float(el[6]))

        oldCoords = getOriginAndVectors()

        if i == 1:
            baseStepVectors.append(oldCoords)

        newCoords = base_step_final.calcNewBSVectors(baseStepVectors[-1][0], baseStepVectors[-1][1], baseStepVectors[-1][2], baseStepVectors[-1][3], Dx, Dy, Dz, Rx, Ry, Rz)
        baseStepVectors.append(newCoords)

        mol = build(i, mol, newCoords)

        ladder += mol

        ladder.SetChainsPerceived()

    conv.WriteFile(ladder, "ladder.pdb")
