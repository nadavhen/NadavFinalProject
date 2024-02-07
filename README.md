# Nadav Final Project
in my lab we have countred a problem were we need to design [PCR primers](https://www.youtube.com/watch?v=NODrmBHHni8&ab_channel=Henrik%27sLab) but the tools on the market these days are expensive and so the planning is done manually and inefficiently.
when planned manualy you need to select the sequence yourself and then analyze it using various online tools to get GC content, secondary structure and [Primer dimer](https://kilobaser.com/the-pain-of-primer-dimer/) which will all have effect on the PCR reaction efficiency.
I want to build a tool that will get a sequence and produce primers for that sequence.

In order to do so i will be using these libraries and modules:

<ul>
    <li>Biopython<l/i>
    <li>foldseq</li>
    <li>ViennaRNA Package 2</li>
<ul>



the input for the program:

<ol>
    <li>the sequence in raw text or fasta format</li>
    <li>from which nucleutides to start and end (optional by default it will create primers for the entire sequence)</li>
    <li>anealing temperature (optional difault will be 55C)</li>
    <li>minimal free energy of secondary structure (optional difault will be set to -10)</li>
</ol>

the out put would be an excel file containing the primers with all their parameters and secondary structure

