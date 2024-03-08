# Nadav Final Project

in my lab we have countred a problem were we need to design [PCR primers](https://www.youtube.com/watch?v=NODrmBHHni8&ab_channel=Henrik%27sLab) but the tools on the market these days are expensive so the planning is done manually and inefficiently.

when planned manualy you need to select the sequence yourself and then analyze it using various online tools to get GC content, secondary structure and [Primer dimer](https://kilobaser.com/the-pain-of-primer-dimer/) which will all have effect on the PCR reaction efficiency.

I want to build a tool that will get a sequence and produce primers for that sequence.

In order to do so I am using these libraries and modules:

- [Biopython](https://biopython.org/)
- [ViennaRNA Package 2](https://www.tbi.univie.ac.at/RNA/documentation.html#)
- other standard packges like numpy and pandas

The application can be used to design primers for any DNA sequence, providing a comprehensive solution for primer design and analysis. The graphical user interface simplifies the process of entering sequence data and viewing the results, enhancing the user experience and productivity. The tool is designed to be intuitive and user-friendly, making it accessible to researchers and professionals involved in molecular biology and genetics.

The application is designed to be modular and extensible, allowing for the addition of new features and functionalities in the future. It can be integrated with other tools and platforms to enhance its capabilities and provide a more comprehensive solution for primer design and analysis. The codebase is well-documented and structured, making it easy to understand and maintain, and the application can be further optimized for performance and scalability.


Step 1: Install all required packages whcih include: [Biopython](https://biopython.org/), - [ViennaRNA Package 2](https://www.tbi.univie.ac.at/RNA/documentation.html#), [tkinter](https://docs.python.org/3/library/tkinter.html), [pandas](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html).

```python
pip install biopython
pip install viennarna
pip install tkinter
pip install pandas
```

Step 2: download python files: Main.py, Primerdesigner.py, GUI.py and place them in the same directory.

Step 3: Run the Main.py file using the command line or an IDE. After running the file the next window will open up:
<div align="center">
    <img src="Pictures/start_menu.png" alt="Start Menu" width="400">
</div>
Step 4: Enter the DNA sequence or upload a FASTA file to begin the analysis.

Step 5: Customize primer design parameters and review the results for further use.




Step 4: Enter the DNA sequence or upload a FASTA file to begin the analysis.

Step 5: Customize primer design parameters and review the results for further use.

