# 3D-BMPP refactor

This is a refactor of the  3D Beta-barrel Membrane Protein Predictor code by Wei Tian. 
The original source code is available at : http://sts.bioe.uic.edu/3dbmpp/

### Full Citation
```
    @article {Tian201716817,
        author = {Tian, Wei and Lin, Meishan and Tang, Ke and Liang, Jie and Naveed, Hammad},
        title = {High-resolution structure prediction of β-barrel membrane proteins},
        year = {2018},
        doi = {10.1073/pnas.1716817115},
        publisher = {National Academy of Sciences},
        issn = {0027-8424},
        URL = {http://www.pnas.org/content/early/2018/01/25/1716817115},
        eprint = {http://www.pnas.org/content/early/2018/01/25/1716817115.full.pdf},
        journal = {Proceedings of the National Academy of Sciences}
    }
```

## Requirements
### SCRWL
Sidechain prediction relies on [SCWRL](http://dunbrack.fccc.edu/scwrl4/). 
To download scwrl, you must apply for a license on their website. 
We then recommend unpacking it into a folder name `scwrl`in the top level folder of this repo. 
You can also change the path to an absolute path on line 38 of `defs.py`.

You also need to chage line 3 in `scrwl/Scwrl4.ini` to `FilePath = [absolute path to scwrl]/bbDepRotLib.bin`

### Python
Python requirements: `numpy, pandas, biopython`

`pip install -r requirements.txt` to install required python packages.

The code was developed using `python 2.7` but should be compatible with python 3.

## Compiling
Part of the application is written in C and needs to be compiled. In a terminal from the main folder:

```
cd bb_register/pred_combine/
make
```

If you're getting an error while trying to run the program that says something like 
`Can't find file construction_inputs/results/regs/*.regs` you probably forgot to do this step.

## Usage

To predict beta-barrel structure, you must provide a fasta file with the protein sequence 
as well as the identified transmembrane strands. 
See the example at the end of this readme and `example_fasta.txt` for reference.
You can also provide a fasta file with multiple sequences and they will be processed individually. 

**The transmembrane strands are NOT predicted by our tool.**
We recommend using [PRED-TMBB2](http://www.compgen.org/tools/PRED-TMBB2) to predict the transmembrane strands beforehand.

These files will be processed in order to create the appropriate files for the predictor in the `inputs/` directory.
For example, if we predict structure for the example fast file which contains a protein named `A9WGN5_CHLAA`
the `inputs/` directory will end up looking like this:
```
|-inputs/
    |- A9WGN5_CHLAA/
        |- A9WGN5_CHLAA.res
        |- A9WGN5_CHLAA.strands
        |- A9WGN5_CHLAA.tmregs
        |- A9WGN5_CHLAA.tmstrands
        
```
Feel free to clean up the inputs directory between runs.

### Running the code 

From a shell terminal call: `python main.py -f [fastafile]`

to predict the beta barrel structure.  The newly created pdb files can be found in `output/` and can be inspected using [PyMOL](https://pymol.org/2/). 

#### Other parameters:
`--level` or `-l` : \[1-5\] or -1 
 *This changes the set of hyperparameters used during register prediction.

A detailed description of these parameters can be found on page 2 of the appendix in `docs/`.
If the parameter is set to -1, a predicted structure is generated with parameters for each level.
The level of parameters used is identified in output files by a `_l0[#]` suffix.

When the parameter is not defined, we infer the level based on the number of transmembrane strands (N) in the provided input file. 
* For N<16, 2 barrels are predicted using both levels 1 and 2
* For 16 <= N < 20, 2 barrels are predicted using levels 3 and 4

Quick reference :

* 1 : Small BMPs (N < 16) without inplugs or outclamps
* 2 : Small BMPs (N < 16) with inplugs or outclamps
* 3 : Medium oligomeric BMPs (16 <= N < 20)
* 4 : Medium monomeric BMPs (16 <= N < 20)
* 5 : Large BMPs (N >= 20)


### Example fasta file
```
>tr|A9WGN5|A9WGN5_CHLAA
MTMINRSRLTIFALLLTGILGSIIAIWSWSANAQTASLTVSPTVARQNTTVTLYGSGFVP
GEKVSIWITYPDYTVYGVTVLTIDERGQFSHPYLPDFLGATFTPTGRYTYTARGWQSGRE
AYASIDVDIAPAPGTTAGVQLTVDRAVQTQGNTFTFSGSGYKPGERVALWLRYPNNAVAD
LGVQIADGQGRIGLAIDSNGVPVGRYALTARGLQSGGNGIVEFEVQVGDALRPRGTAGLE
VGPGSSQQRSAVSLRGTGFLPGEVITIWATRPDYSTEWLGDVTAAADGSFTTELYLSEQN
PAGRYAFSAYGNRSERRAVAEYTLLPGR
>tr|A9WGN5|A9WGN5_CHLAA
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIMMMMMMMMMOOOOOOMMMMMMMMMI
IIMMMMMMMMOOOOOOOOOOOOOOOOOOMMMMMMMMMIIMMMMMMMMMMMOOOOOOOOOO
MMMMMMMMMIIIIIIIMMMMMMMMMMMOOOOOOOOOOOOOOOOOOOMMMMMMMIIIMMMM
MMMMMOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOMMMMMMMMMIIIMMMMMMMMMM
MOOOOOOOOOOOOMMMMMMMIIMMMMMMMMMOOOOOOMMMMMMMMMIIMMMMMMMMOOOO
OOOOOOOOOOOOOOOOOMMMMMMMMMMM
```

The two lines starting with `>tr|[short]|[full]` indicate a shortened version of the protein name and its full name. 
The first is followed by a line containing the amino acid sequence. The second identifying line is followed by the 
transmembrane strand prediction, where each amino acid is predicted as being either 
inner (`I`), outer (`O`), or transmembrane(`M`). 


### Adding loops with MODELLER
After building the beta barrels, the program will print the command to run Modeller for adding the
loops to the transmembrane strands already modelled. By default, 100 models are generated but this
can be changed with the -n parameter. We recommend looking at the generated transmembrane strands if
you have used more than one "level" parameter, and choosing the one that seems the most realistic beta-barrel
for modelling the loops.
`python ./LoopModeller.py -b A9WGN5_CHLAA_ext_l03 -l 3 -n 100`

