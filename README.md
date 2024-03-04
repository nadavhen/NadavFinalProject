# Nadav Final Project

in my lab we have countred a problem were we need to design [PCR primers](https://www.youtube.com/watch?v=NODrmBHHni8&ab_channel=Henrik%27sLab) but the tools on the market these days are expensive so the planning is done manually and inefficiently.

when planned manualy you need to select the sequence yourself and then analyze it using various online tools to get GC content, secondary structure and [Primer dimer](https://kilobaser.com/the-pain-of-primer-dimer/) which will all have effect on the PCR reaction efficiency.

I want to build a tool that will get a sequence and produce primers for that sequence.

In order to do so I am using these libraries and modules:

- [Biopython](https://biopython.org/)
- [ViennaRNA Package 2](https://www.tbi.univie.ac.at/RNA/documentation.html#)
- other standard packges like numpy and pandas

The application features a graphical user interface for inputting DNA sequences for analysis. It supports direct text entry or uploading of FASTA files, with validation to ensure input accuracy. The interface dynamically adjusts based on the chosen input method, streamlining the user experience. This tool is designed to facilitate DNA sequence analysis by providing a user-friendly platform for data input and validation.
